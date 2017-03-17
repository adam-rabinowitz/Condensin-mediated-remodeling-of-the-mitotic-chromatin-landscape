import gzip
import multiprocessing
import subprocess
import os
import random
import itertools

def FastqGeneralIterator(handle): 
    ''' Function taken from Bio.SeqIO.Quality.IO '''
    handle_readline = handle.readline
    while True: 
        line = handle_readline() 
        if not line: 
            return
        if line[0] == "@": 
            break 
        if isinstance(line[0], int): 
            raise ValueError("Is this handle in binary mode not text mode?")
    while line: 
        if line[0] != "@": 
            raise ValueError( 
                "Records in Fastq files should start with '@' character") 
        title_line = line[1:].rstrip()
        seq_string = handle_readline().rstrip() 
        while True: 
            line = handle_readline() 
            if not line: 
                raise ValueError("End of file without quality information.") 
            if line[0] == "+": 
                second_title = line[1:].rstrip() 
                if second_title and second_title != title_line: 
                    raise ValueError("Sequence and quality captions differ.") 
                break 
            seq_string += line.rstrip()  # removes trailing newlines 
        seq_len = len(seq_string) 
        quality_string = handle_readline().rstrip() 
        while True: 
            line = handle_readline() 
            if not line: 
                break
            if line[0] == "@": 
                if len(quality_string) >= seq_len: 
                    break 
                quality_string += line.rstrip() 
        if seq_len != len(quality_string): 
            raise ValueError("Lengths of sequence and quality values differs " 
                " for %s (%i and %i)." 
                % (title_line, seq_len, len(quality_string))) 
        yield (title_line, seq_string, quality_string) 
    raise StopIteration 

class parseFastq(object):
    
    ''' Class functions as an iterator to extract reads from single or paired
    FASTQ files. Classes uses multiprocessing to speed up the extraction of
    the FASTQ entries and the FastqGeneralIterator from Bio.SeqIO to parse
    indivdual reads.
    '''
    
    def __init__(
        self, fastq1, fastq2, shell = True
    ):
        ''' Function to initialise readFastq object. Checks FASTQ files
        exist and creates lists to store read and write processes.
        
        Args:
            fastq1 (str)- Full path to read1 FASTQ file or list of paths.
            fastq2 (str)- Full path to read2 FASTQ file or list of paths.
            shell (bool)- Whether to use shell to read gzip files.
        
        '''
        # Store fastq files
        fastq_list = []
        # Process fastq1 string argument
        if isinstance(fastq1, str):
            if isinstance(fastq2, str):
                fastq_list.append((fastq1, fastq2))
            else:
                raise TypeError('fastq2 must be same type as fastq1')
        # Process fastq 1
        if isinstance(fastq1, list):
            if not isinstance(fastq2, list):
                raise TypeError('fastq2 must be None or same type as fastq1')
            if len(fastq1) != len(fastq2):
                raise ValueError('fastq1 and fastq2 have different lengths')
            for f1, f2 in zip(fastq1, fastq2):
                fastq_list.append((f1, f2))
        # Raise error
        else:
            raise TypeError('fastq1 must be string or list')
        # Check if fastq files exist
        for fastq_tuple in fastq_list:
            for fastq in fastq_tuple:
                if not os.path.isfile(fastq):
                    raise IOError('File {} could not be found'.format(fastq))
        # Store data
        self.fastq_list = fastq_list
        self.shell = shell
        self.read_processes = []
    
    def __read_handle_create(self, fastq):
        ''' Function to create filehandle and subprocess for reading of
        FASTQ files.
        
        Args:
            fastq (str)- Full path to fastq file.
        
        Returns:
            fh - File handle for fastq file.
            sp - zcat subprocess reading the FASTQ file if self.shell=True
                and files ends with '.gz', else None.
        
        '''
        # Create fastq handle and subprocess
        if self.shell and fastq.endswith('.gz'):
            sp = subprocess.Popen(['zcat', fastq], stdout = subprocess.PIPE,
                bufsize = 1)
            fh = sp.stdout
        elif fastq.endswith('.gz'):
            sp = None
            fh = gzip.open(fastq)
        else:
            sp = None
            fh = open(fastq)
        # Return handle and subprocess
        return(fh, sp)
    
    def __send(self, data, conn):
        ''' Function to send data down the end of multiprocessing pipe.
        First checks that termination signal, None, has not been received
        first.
        
        Args:
            data - Data to send down pipe.
            conn - End of pipe.
            
        Raises:
            StopIteration: If termination signal of None is received.
            ValueError: If anythind other than None is received.
       
        '''
        # Check pipe for incoming signal
        if conn.poll():
            # Extract signal and process
            recv = conn.recv()
            if recv is None:
                raise StopIteration('Termination signal received')
            else:
                raise ValueError('Unknown signal received')
        # Add read to connection
        else:
            conn.send(data)
    
    def __read_processes_recv(self):
        ''' Function to receive data from running processes.
        
        Returns:
            data - A list of data received from each process.
        
        Raises:
            IOError - If not all processes return data.
            EOFError - If all processes have no data
        
        '''
        # Extract signal and process
        data = []
        try:
            data.append(self.read_processes[0][1].recv())
        except EOFError:
            try:
                self.read_processes[1][1].recv()
                raise IOError('Differing number of data points')
            except EOFError:
                raise EOFError
        try:
            data.append(self.read_processes[1][1].recv())
        except EOFError:
            raise IOError('Differing number of data points')
        return(data)
    
    def __read_process_stop(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. The list 
        self.process_list is then emptied.
        '''
        # Loop through processes and terminate
        for process, conn in self.read_processes:
            # Add termination signal to pipes
            try:
                conn.send(None)
            except IOError:
                pass
            # Join process and close pipes
            process.join()
            conn.close()
        # Empty process list
        self.read_processes = []
    
    def __interleave_trim_read_process(
            self, trim, number, conn
        ):
        ''' Function to trim fastq reads until after supplied sequence.
        
        Args:
            fastq (str)- Full path to the FASTQ file to read.
            trim (str)- Sequence after which sequence is trimmed.
            pend - End of multiprocessing pipe down which data will be sent.
        
        Sends:
            name (str)- Read name
            number (int)- Read number
        
        '''
        # Process pipes
        conn1, conn2 = conn
        conn1.close()
        # Extract numbe for string find adjustment
        adj = len(trim)
        # Loop through FASTQ files
        stop = False
        for pair in self.fastq_list:
            if stop:
                continue
            fastq = pair[number]
            # Loop through fastq files
            fh, sp = self.__read_handle_create(fastq)
            for read in FastqGeneralIterator(fh):
                # Extract read data
                header, sequence, quality = read
                name, description = header.split(None, 1)
                readno, other = description.split(':', 1)
                # Find trim sequence
                trimLoc = sequence.find(trim)
                if trimLoc != -1:
                    trimLoc += adj
                    sequence = sequence[:trimLoc]
                    quality = quality[:trimLoc]
                # Send data or break iteration
                read = (
                    '{}:{} {}'.format(name, readno, description),
                    sequence,
                    quality
                )
                try:
                    self.__send((name, readno, trimLoc, read), conn2)
                except StopIteration:
                    stop = True
                    break
            # Clean up
            if sp:
                sp.terminate()
            fh.close()
        # Close connection
        conn2.close()
    
    def interleave_trim_reads(
        self, trim, outFastq, minLength = 20
    ):
        ''' Function interleaves paired fastq files into a single fastq
        file.
        
        Args:
            label (bool)- Add ':1' label to read1 name and ':2' label to
                read2 name.
            check_pairs (bool)- Check read pairing prior to interleaving.
        
        Returns:
            count (int)- Number of paired reads processed
        
        '''
        # Loop through fastq files:
        for number in (0, 1):
            # Create pipe and process
            conn = multiprocessing.Pipe(True)
            process = multiprocessing.Process(
                target = self.__interleave_trim_read_process,
                args = (trim, number, conn)
            )
            process.start()
            conn[1].close()
            # Store process data
            self.read_processes.append((process, conn[0]))
        # Extract data and count reads
        metrics = {'total':0, 'short':0, 'trim1':0, 'trim2':0}
        with writeFastq(outFastq, None, self.shell) as fastqOut:
            while True:
                try:
                    data = self.__read_processes_recv()
                except EOFError:
                    break
                metrics['total'] += 1
                # Check names and read number
                if data[0][0] != data [1][0]:
                    self.__read_process_stop()
                    raise ValueError('Read name mismatch')
                # Check read number
                if (data[0][1] != '1'
                    or data[1][1] != '2'):
                    self.__read_process_stop()
                    print(data)
                    raise ValueError('Unexpected read numbers')
                # Count and skip short reads
                if (-1 < data[0][2] < minLength
                    or -1 < data[1][2] < minLength):
                    metrics['short'] += 1
                    continue
                # Write output fastq
                if data[0][2] != -1:
                    metrics['trim1'] += 1
                if data[1][2] != -1:
                    metrics['trim2'] += 1
                fastqOut.write(data[0][3])
                fastqOut.write(data[1][3])
        # Stop read processes and return count
        self.__read_process_stop()
        return(metrics)
    
class writeFastq(object):
    ''' An object that uses multiprocessing processes to parralelize the
    writing of FASTQ files. It assumes that the fastq files to write will
    be generated with the FastqGeneralIterator from the ??? package.
    Therefore read names will not have a leading '@' symbol which is
    therefore added by the write function.
    '''
    
    def __init__(self, fastq1, fastq2 = None, shell = True):
        ''' Function to initialise object. Function takes three arguments:
        
        1)  fastq1 - Full path to FASTQ file.
        2)  fastq2 - Full path to paired FASTQ file (optional).
        3)  shell - Boolean indicating whether to use shell gzip command
            to write gzipped output.
        
        '''
        # Store fastq files
        if not fastq2 is None:
            self.fastq_list = [fastq1, fastq2]
            self.pair = True
        else:
            self.fastq_list = [fastq1]
            self.pair = False
        ## Check output directories
        for fastq in self.fastq_list:
            outDir = os.path.dirname(fastq)
            if not os.path.isdir(outDir):
                raise IOError('Could not find output directory')
        # Store shell argument and process list
        self.shell = shell
        self.process_list = []
    
    def __write_process(self, fastq, conn):
        ''' Function to generate a process to write FASTQ files. FASTQ reads
        received from the multiprocessing pipe will be written to file. Receipt
        of None will cause the termination of the process. If self.shell is True
        then gzipped output files will be written using the gzip command in the
        shell. Funcion takes two arguments:
        
        1)  fastq - Full path to the FASTQ file to be created
        2)  pend - End of multiprocessing pipe from which reads will be
            extracted.
        '''
        # Process connections
        conn[1].close
        # Write gzip file using shell
        if self.shell and fastq.endswith('.gz'):
            # Create process
            command = 'gzip -c > %s' %(fastq)
            sp = subprocess.Popen(command, shell=True,
                stdin = subprocess.PIPE)
            fh = sp.stdin
        # Write file using python
        elif fastq.endswith('.gz'):
            sp = None
            fh = gzip.open(fastq, 'w')
        else:
            sp = None
            fh = open(fastq, 'w')
        # Extract reads from pipe and write to file
        while True:
            # Extract read
            read = conn[0].recv()
            # Break loop if read is none
            if read is None:
                break
            # Write read to file
            readString = '@{}\n{}\n+\n{}\n'.format(read[0], read[1], read[2])
            fh.write(readString)
        # Close files, pipes and subprocesses
        fh.close()
        if sp:
            sp.communicate()
        conn[0].close()
    
    def start(self):
        ''' Function creates processes to write the FASTQ files listed in
        self.fastq_list using self.__write_process. Function creates two
        element tuples consisting of the process and pipe end with which
        to communicate with the process. The tuples are stored in
        self.process_list.
        '''
        # Close active processes
        self.close()
        # Loop through fastq files
        for fastq in self.fastq_list:
            # Create pipe and process
            conn = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = self.__write_process,
                args = (fastq, conn)
            )
            process.start()
            conn[0].close()
            # Store pipe end and processes
            self.process_list.append((process, conn[1]))
    
    def close(self):
        ''' Function terminates the processes and closes the pipe-ends
        listed in self.process_list and generated by self.start. Thelist 
        self.process_list is then emptied.
        '''
        # Extract process and pipes
        for process, conn in self.process_list:
            # Add poisin pill, join process and close pipe
            try:
                conn.send(None)
            except IOError:
                pass
            process.join()
            conn.close()
        # Clear pipe and process list
        self.process_list = []
    
    def __pipe_send(self, read, conn):
        ''' Function to send a FASTQ read down a multiprocessing pipe. If
        the read does not correspond to the expected format all active
        processes are terminated and an IOError is raised. Function takes
        two arguments:
        
        1)  read - FASTQ read. This should consist of a tuple/list of
            three elements: read name, sequence, quality string.
        2)  pend - End of multiprocessing pipe down which read should be
            sent.
        '''
        # Add read1 to pipe
        if isinstance(read, (tuple,list)) and len(read) == 3:
            conn.send(read)
        else:
            self.close()
            raise IOError('Read must be a list/tuple of three elements')
    
    def write(self, reads):
        ''' Function to write paired reads to paired FASTQ files. Function
        takes one argument:
        
        1)  reads - For single FASTQ files a single read should be provided
            while for paired FASTQ files a tuple/list of 2 reads should be
            provided. Each read should consist of a tuple/list of three
            strings: read name, sequence and quality
        '''
        # Process pairs
        if self.pair:
            # Check input consists of two elements
            if len(reads) != 2:
                self.close()
                raise IOError('Write requires two elements for paired FASTQ')
            # Send elements to write process
            self.__pipe_send(reads[0], self.process_list[0][1])
            self.__pipe_send(reads[1], self.process_list[1][1])
        # Process single reads
        else:
            # Send read to pipe
            self.__pipe_send(reads, self.process_list[0][1])
    
    def __enter__(self):
        ''' Start processes upon entry into with scope '''
        self.start()
        return(self)
    
    def __exit__(self, type, value, traceback):
        ''' Terminate processes upon exit of with scope '''
        self.close()
