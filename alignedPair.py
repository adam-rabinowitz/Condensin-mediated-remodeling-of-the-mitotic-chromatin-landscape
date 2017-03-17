# Import required modules
from Bio import SeqIO
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq
import bisect
import collections
import gzip
import multiprocessing
import pysam
import subprocess

class writeFileProcess(object):
    
    ''' Creates a class to handle writing to file in a seperate process '''
    
    # Store parental pipe ends
    pipeSendList = []
     
    def __init__(
            self, fileName, shell = True
        ):
        ''' Initiate object.
        
        Args:
            fileName - Full path to output file.
            shell - Boolean, indicating whether gzip fie should be written
                using shell 'gzip' command and the python subprocess module.
        
        '''
        # Store supplied parameters
        self.fileName = fileName
        self.shell = shell
        # Create pipes
        self.pipes = multiprocessing.Pipe('False')
        # Create and start process
        self.process = multiprocessing.Process(target = self.writeFromPipe)
        self.process.start()
        # Process pipes
        self.pipes[0].close()
        writeFileProcess.pipeSendList.append(self.pipes[1])
    
    def writeFromPipe(self):
        ''' Function to write to gzip files within a python
        multiprocessing process.
        '''
        # Close unused pipes
        self.pipes[1].close()
        for p in writeFileProcess.pipeSendList:
            p.close()
        # Open output file
        if self.fileName.endswith('.gz'):
            outFile = gzip.open(self.fileName, 'wb')
        else:
            outFile = open(self.fileName, 'w')
        # Create output gzip file using subprocess
        if self.shell and self.fileName.endswith('.gz'):
            # Create gzip subprocess
            sp = subprocess.Popen('gzip', stdout = outFile,
                stdin = subprocess.PIPE, close_fds=True)
            # Write to output subprocess
            while True:
                try:
                    line = self.pipes[0].recv()
                except EOFError:
                    break
                sp.stdin.write(line)
            # Terminate subprocess
            sp.communicate()
        # Or create output file using pure python
        else:
            # Write to output file
            while True:
                try:
                    line = self.pipes[0].recv()
                except EOFError:
                    break
                outFile.write(line)
        # Close files and pipes
        outFile.close()
        self.pipes[0].close()
    
    def __enter__(self):
        return(self)
    
    def add(self, data):
        ''' pass data to write process
        
        Args:
            data (str)- data to write.
        
        '''
        if not isinstance(data, str):
            self.close()
            raise ValueError('Object to write must be string')
        self.pipes[1].send(data)
    
    def close(self):
        ''' Close outfile and write process. '''
        self.pipes[1].close()
        self.process.join()
    
    def __del__(self):
        self.close()
    
    def __exit__(self, type, value, traceback):
        self.close()

class alignedPair(object):
    
    ''' An object to handle the generation of aligned pairs from a name
    sorted bam file. '''
    
    def __init__(self, bam):
        self.bam = bam
    
    def __concordant(self, reads, maxSize):
        ''' Function to find concordant pairs. Input is a list/tuple
        that sequentially contains the chromosome, start, end and
        strand. Start is the most 5' base on the genome to which the
        read aligns. End is the most 3' base on the genome to which
        the read aligns.
        '''
        # Only examine reads on same chromsomes and different strands
        if reads[0] == reads[4] and reads[3] != reads[7]:
            # Process read1 forward strand
            if reads[3] == '+':
                # Calculate and check read distance:
                distance = reads[6] - reads[1]
                if distance < maxSize:
                    # Check that no read extends beyond its pair
                    if (reads[5] >= reads[1] and
                        reads[6] >= reads[2]):
                        return(True)
            # Process read1 on reverse strand
            else:
                # Calculate and check read distance:
                distance = reads[2] - reads[5]
                if distance < maxSize:
                    # Check that no read extends beyond its pair
                    if (reads[2] >= reads[6] and
                        reads[1] >= reads[5]):
                        return(True)
        # Retrun return variable
        return(False)
    
    def __processPairs(
            self, pipe, pairOut, rmDup, rmConcord, maxSize
        ):
        ''' Function to output read pairs generated from the extract
        function while processing concordant and duplicate reads.
        
        Args:
            pipe - A multriprocessing.Pipe connection
            pairOut (str)- Path to output file.
            rmDup (bool)- Remove duplicate reads.
            rmConcord (bool)- Remove concirdant reads.
            maxSize (int)- Maximum size for concordant pairs.
        
        Returns:
            pairCount - A dictionary containing processing metrics.
        
        '''
        # Create counter and pair set
        pairCount = collections.defaultdict(int)
        pairSet = set()
        # Open output file process
        outObject = writeFileProcess(fileName = pairOut)
        # Loop through pairs
        while True:
            # Get pair from pipe
            pair = pipe.recv()
            if pair == None:
                break
            # Count and check for duplicates pairs
            pairCount['total'] += 1
            if pair in pairSet:
                dup = True
                pairCount['duplicate'] += 1
            else:
                dup = False
                pairCount['unique'] += 1
                pairSet.add(pair)
            # Count and check for concordant pairs
            concord = self.__concordant(pair, maxSize)
            if concord:
                pairCount['concord'] += 1
                if not dup:
                    pairCount['concorduni'] += 1
            else:
                pairCount['discord'] += 1
                if not dup:
                    pairCount['discorduni'] += 1
            # Process output
            if dup and rmDup:
                continue
            elif concord and rmConcord:
                continue
            else:
                outData = '\t'.join(map(str,pair)) + '\n'
                outObject.add(outData)
        # Close file, return data and close pipe
        outObject.close()
        pipe.send(pairCount)
        pipe.close()
    
    def extractPairs(
            self, pairOut, minMapQ = 10, rmDup = True, rmConcord = True,
            maxSize = 2000
        ):
        ''' Function to output read pairs generated from the extract
        function while processing concordant and duplicate reads.
        Function takes five arguments:
        
        1)  inBam - Path to input BAM file.
        2)  minMapQ - minimum mapping quality for a read to be
            processed,
        
        Function returns two items:
        
        1)  A python dictionary where the key is the read pair and the value
            is the frequency at which the read pair is found.
        2)  A python dictionary listing the alignment metrics.
        
        '''
        # Open bamfile
        bamFile = pysam.AlignmentFile(self.bam, 'rb')
        # Generate dictionaries to store and process data
        alignCount = collections.defaultdict(int)
        strDict = {True: '-', False: '+'}
        chrDict = {}
        for r in bamFile.references:
            chrDict[bamFile.gettid(r)] = r
        # Initialise variables to store read data
        currentName = ""
        readList = []
        # Create process to handle pairs
        pipes = multiprocessing.Pipe(True)
        p = multiprocessing.Process(
            target = self.__processPairs,
            args = (pipes[0], pairOut, rmDup, rmConcord, maxSize)
        )
        p.start()
        pipes[0].close()
        # Loop through BAM file
        while True:
            try:
                read = bamFile.next()
                readName = read.query_name
                alignCount['total'] += 1
            except StopIteration:
                readName = 'EndOfFile'
            # Process completed families
            if readName[:-2] != currentName[:-2]:
                # Count number of reads with identical ID
                readNo = len(readList)
                # Count and process properly mapped read-pairs
                if readNo == 2:
                    # Unpack reads and check for read1 and read2
                    read1, read2 = readList
                    if (read1.query_name.endswith(':1') and 
                        read2.query_name.endswith(':2')):
                        # Count pairs and store data
                        output = (
                            chrDict[read1.reference_id],
                            read1.reference_start + 1,
                            read1.reference_end,
                            strDict[read1.is_reverse],
                            chrDict[read2.reference_id],
                            read2.reference_start + 1,
                            read2.reference_end,
                            strDict[read2.is_reverse],
                        )
                        pipes[1].send(output)
                        alignCount['pairs'] += 2
                    # If not, count as multiple alignments
                    else:
                        alignCount['multiple'] += 2
                # Count single mapped and multi mapped reads
                elif readNo == 1:
                    alignCount['singletons'] += 1
                else:
                    alignCount['multiple'] += readNo
                # Reset read list and current name
                currentName = readName
                readList = []
            # Break loop at end of BAM file
            if readName == 'EndOfFile':
                pipes[1].send(None)
                break
            # Count and skip secondary alignments
            elif (256 & read.flag):
                alignCount['secondary'] += 1
            # Count and skip umapped reads
            elif (4 & read.flag):
                alignCount['unmapped'] += 1
            # Count and skip poorly mapped reads
            elif read.mapping_quality < minMapQ:
                alignCount['poormap'] += 1
            # Process reads of sufficient quality
            else:
                readList.append(read)
        # Close BAM file
        bamFile.close()
        # Extract data from process and terminate
        pairCount = pipes[1].recv()
        pipes[1].close()
        p.join()
        # Output data
        return(alignCount, pairCount)

class fragendPair(object):
    
    ''' Class to generate fragend pair files from aligned pair file. '''
    
    def __init__(self, fasta, resite):
        ''' Initialises fragendPair object.
        
        Args:
            fasta (str)- Path to genome fasta file.
            resite (str)- Recognition site for restriction enzyme.
        
        '''
        self.fasta = fasta
        self.resite = resite
        self.frags = self.__findFragendSites()
    
    def __findFragendSites(self):
        ''' Function creates FragendDict object. The object contains
        the location of all fragends for each strand of each chromosome.
        
        Args:
            fasta (str)- Full path to fasta file
            resite (str)- Recognition sequence of restriction enzyme

        '''
        # Process restriction enzyme size and create output dictionary
        resite = self.resite.upper()
        frags = {'resite': resite}
        # Create sequence object for resite and reverse complent
        standard = Seq(resite)
        revcomp = standard.reverse_complement()
        # Open and parse fasta file
        fastaHandle = open(self.fasta)
        fastaData = SeqIO.parse(fastaHandle, 'fasta')
        # Extract fragend information for each chromosome
        for fasta in fastaData:
            # Extract name and sequence
            fName, fSequence = str(fasta.id), str(fasta.seq).upper()
            # Add re sites to dictionary using 1 based index
            forward = nt_search(fSequence, standard)[1:]
            if forward:
                frags[(fName,'+')] = [x + len(resite) for x in forward]
            else:
                frags[(fName,'+')] = []
            reverse = nt_search(fSequence, revcomp)[1:]
            if reverse:
                frags[(fName,'-')] = [x + 1 for x in reverse]
            else:
                frags[(fName,'-')] = []
        # Close input file and return data
        fastaHandle.close()
        return(frags)
    
    def __downstream(self, reads):
        ''' Function to find the associated fragend for a read. For
        reads aligned to the '+' strand the fragend will be downstream
        of the start of the read. For reads aligned to the '-' strand
        the fragend will be upstream of the start of the read. Function
        takes three arguments:
        
        Args:
            reads - An iterator returning a list/tuple with the
                following four elements:
                1) chrom - Chromosome to which the read is aligned.
                2) start - 
                and strand. Start is the most 5' base in the genome to which
                the read is aligned. End is the most 3' base in the genome to
                which the read is aligned
        
        Returns:
            output - A list of tuples containing the following four elements:
                1) chrom - Chromsome to which the read is aligned.
                2) fragLoc - Location of the centre of the fragend site.
                3) strand - Strand of the fragend.
                4) distance - Distance between the start of the read and the
        
        '''
        # Extract resite and create output variable
        output = []
        # Sequentially process in data
        for read in reads:
            # Process variables
            chrom, start, end, strand = read
            start = int(start)
            end = int(end)
            fragLoc = None
            distance = None
            # Check relative location of start and end
            if end < start:
                raise IOError("End must be 'right' of the start")
            # Extract fragends for chromosome and strand
            try:
                fragends = self.frags[(chrom, strand)]
            except KeyError:
                raise IOError('No data for %s strand on %s chromosome' %(
                    strand, chrom))
            # Find location of fragend for forward strand read
            if strand == '+':
                fragIndex = bisect.bisect_left(
                    fragends,
                    end
                )
                # Return default data for invalid index
                if fragIndex < len(fragends):
                    # Find fragend location
                    fragLoc = fragends[fragIndex] - (len(self.resite) / 2)
                    # Find distance between read and fragend
                    distance = (fragLoc - start) + 1
            # Find location of fragend for reverse strand read
            elif strand == '-':
                # Find potential index of downstream fragend
                fragIndex = bisect.bisect_right(
                    fragends,
                    start
                ) - 1
                # Return default data for invalid index
                if fragIndex >= 0:
                    # Find fragend location
                    fragLoc = fragends[fragIndex] + (len(self.resite) / 2)
                    # Find distance between read and fragend
                    distance = (end - fragLoc) + 1
            # Add output data
            if fragLoc == None:
                output.append(None)
            else:
                output.append((chrom, fragLoc, strand, distance))
        # Close IO and return data
        return(output)
    
    def extractFragends(self, pairIn, maxDistance, fragendOut):
        ''' Function identifies fragend pairs.
        
        Args:
            pairIn - Paired read input file
            maxDistance - Maximum distance between read start and RE site.
            pairOut - Name of gzipped output file.
        
        '''
        # Create fragend dictionary and metrics dictionary
        fragendCounts = collections.defaultdict(int)
        fragendCounts['fragDist'] = []
        fragendCounts['ligDist'] = []
        # Open input and output file
        if pairIn.endswith('.gz'):
            inFile = gzip.open(pairIn, 'r')
        else:
            inFile = open(pairIn, 'r')
        outFile = writeFileProcess(fragendOut)
        for pair in inFile:
            pair = pair.strip().split('\t')
            # Count entries
            fragendCounts['total'] += 1
            # Create output containg fragend data
            output = self.__downstream([pair[0:4], pair[4:8]])
            # Skip reads without identified fragends
            if output[0] == None or output[1] == None:
                fragendCounts['none'] += 1
                continue
            # Add fragend distance data for pairs with fragends
            fragendCounts['fragDist'].extend([
                output[0][3],
                output[1][3] 
            ])
            # Count and skip reads too distant from the fragend
            if output[0][3] > maxDistance or output[1][3] > maxDistance:
                fragendCounts['distant'] += 1
                continue
            # Save to file accepted ligation pairs
            outData = '\t'.join(map(str,output[0][0:3] + output[1][0:3]))
            outFile.add(outData + '\n')
            # Count interchromosomal ligations 
            if output[0][0] != output[1][0]:
                fragendCounts['interchromosomal'] += 1
                # Count intrachromosomal ligations and store distance
            else:
                fragendCounts['intrachromosomal'] += 1
                fragendCounts['ligDist'].append(
                    abs(output[0][1] - output[1][1])
                )
        # Close files and return data
        inFile.close()
        outFile.close()
        return(fragendCounts)
    
    
#def concordant(reads, maxSize):
#    ''' Function to find concordant pairs. Input is a list/tuple
#    that sequentially contains the chromosome, start, end and
#    strand. Start is the most 5' base on the genome to which the
#    read aligns. End is the most 3' base on the genome to which
#    the read aligns.
#    '''
#    # Initialise return variable
#    returnVariable = False
#    # Only examine reads on same chromsomes and different strands
#    if reads[0] == reads[4] and reads[3] != reads[7]:
#        if reads[3] == '+':
#            # Calculate and check read distance:
#            distance = reads[6] - reads[1]
#            if distance < maxSize:
#                # Check that no read extends beyond its pair
#                if (reads[5] >= reads[1] and
#                    reads[6] >= reads[2]):
#                    returnVariable = True
#        elif reads[3] == '-':
#            # Calculate and check read distance:
#            distance = reads[2] - reads[5]
#            if distance < maxSize:
#                # Check that no read extends beyond its pair
#                if (reads[2] >= reads[6] and
#                    reads[1] >= reads[5]):
#                    returnVariable = True
#    # Retrun return variable
#    return(returnVariable)
#
#def processPairs(pipe, pairOut, rmDup, rmConcord, maxSize):
#    ''' Function to output read pairs generated from the extract
#    function while processing concordant and duplicate reads.
#    
#    Args:
#        pipe - A multriprocessing.Pipe connection
#        pairOut (str)- Path to output file.
#        rmDup (bool)- Remove duplicate reads.
#        rmConcord (bool)- Remove concirdant reads.
#        maxSize (int)- Maximum size for concordant pairs.
#    
#    Returns:
#        pairCount - A dictionary containing processing metrics.
#    
#    '''
#    # Create counter and pair set
#    pairCount = collections.defaultdict(int)
#    pairSet = set()
#    # Open output file process
#    outObject = writeFileProcess(fileName = pairOut)
#    # Loop through pairs
#    while True:
#        # Get pair from pipe
#        pair = pipe.recv()
#        if pair == None:
#            break
#        # Count and check for duplicates pairs
#        pairCount['total'] += 1
#        if pair in pairSet:
#            dup = True
#            pairCount['duplicate'] += 1
#        else:
#            dup = False
#            pairCount['unique'] += 1
#            pairSet.add(pair)
#        # Count and check for concordant pairs
#        concord =  concordant(pair, maxSize)
#        if concord:
#            pairCount['concord'] += 1
#            if not dup:
#                pairCount['concorduni'] += 1
#        else:
#            pairCount['discord'] += 1
#            if not dup:
#                pairCount['discorduni'] += 1
#        # Process output
#        if dup and rmDup:
#            continue
#        elif concord and rmConcord:
#            continue
#        else:
#            outData = '\t'.join(map(str,pair)) + '\n'
#            outObject.add(outData)
#    # Close file, return data and close pipe
#    outObject.close()
#    pipe.send(pairCount)
#    pipe.close()
#
#def extractPairs(inBam, pairOut, minMapQ, rmDup, rmConcord, maxSize):
#    ''' Function to output read pairs generated from the extract
#    function while processing concordant and duplicate reads.
#    Function takes five arguments:
#    
#    1)  inBam - Path to input BAM file.
#    2)  minMapQ - minimum mapping quality for a read to be
#        processed,
#    
#    Function returns two items:
#    
#    1)  A python dictionary where the key is the read pair and the value
#        is the frequency at which the read pair is found.
#    2)  A python dictionary listing the alignment metrics.
#    
#    '''
#    # Open bamfile
#    bamFile = pysam.AlignmentFile(inBam, 'rb')
#    # Generate dictionaries to store and process data
#    alignCount = collections.defaultdict(int)
#    strDict = {True: '-', False: '+'}
#    chrDict = {}
#    for r in bamFile.references:
#        chrDict[bamFile.gettid(r)] = r
#    # Initialise variables to store read data
#    currentName = ""
#    readList = []
#    # Create process to handle pairs
#    pipes = multiprocessing.Pipe(True)
#    p = multiprocessing.Process(
#        target = processPairs,
#        args = (pipes[0], pairOut, rmDup, rmConcord, maxSize)
#    )
#    p.start()
#    pipes[0].close()
#    # Loop through BAM file
#    while True:
#        try:
#            read = bamFile.next()
#            readName = read.query_name
#            alignCount['total'] += 1
#        except StopIteration:
#            readName = 'EndOfFile'
#        # Process completed families
#        if readName[:-2] != currentName[:-2]:
#            # Count number of reads with identical ID
#            readNo = len(readList)
#            # Count and process properly mapped read-pairs
#            if readNo == 2:
#                # Unpack reads and check for read1 and read2
#                read1, read2 = readList
#                if (read1.query_name.endswith(':1') and 
#                    read2.query_name.endswith(':2')):
#                    # Count pairs and store data
#                    output = (
#                        chrDict[read1.reference_id],
#                        read1.reference_start + 1,
#                        read1.reference_end,
#                        strDict[read1.is_reverse],
#                        chrDict[read2.reference_id],
#                        read2.reference_start + 1,
#                        read2.reference_end,
#                        strDict[read2.is_reverse],
#                    )
#                    pipes[1].send(output)
#                    alignCount['pairs'] += 2
#                # If not, count as multiple alignments
#                else:
#                    alignCount['multiple'] += 2
#            # Count single mapped and multi mapped reads
#            elif readNo == 1:
#                alignCount['singletons'] += 1
#            else:
#                alignCount['multiple'] += readNo
#            # Reset read list and current name
#            currentName = readName
#            readList = []
#        # Break loop at end of BAM file
#        if readName == 'EndOfFile':
#            pipes[1].send(None)
#            break
#        # Count and skip secondary alignments
#        elif (256 & read.flag):
#            alignCount['secondary'] += 1
#        # Count and skip umapped reads
#        elif (4 & read.flag):
#            alignCount['unmapped'] += 1
#        # Count and skip poorly mapped reads
#        elif read.mapping_quality < minMapQ:
#            alignCount['poormap'] += 1
#        # Process reads of sufficient quality
#        else:
#            readList.append(read)
#    # Close BAM file
#    bamFile.close()
#    # Extract data from process and terminate
#    pairCount = pipes[1].recv()
#    pipes[1].close()
#    p.join()
#    # Output data
#    return(alignCount, pairCount)
