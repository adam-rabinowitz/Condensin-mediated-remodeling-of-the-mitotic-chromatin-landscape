import collections
import multiprocessing
import os
import subprocess
import gzip
import cStringIO
import numpy as np
import re

class genomeBin(object):
    
    ''' A class to generate raw count matrices from fragend ligation files '''
    
    def __init__(self, binData):
        ''' Function to initialise genomeBin object.
        
        Args:
            binData - Either a string or a tuple. The string should be a path
                to a tad-delimited text file containing the following three
                columns:
                    1) Chromosome name
                    2) Bin start position. 1-based numbering. Closed interval.
                    3) Bin end position. 1-based numbering. Closed interval.
                The tuple should contain the following elements.
                    1) Tab delimited file of chromosome sizes.
                    2) Bin size.
                    3) Bool indicating bins sizes should be equal.
        
        '''
        self.binData = binData
        if isinstance(self.binData, str):
            self.binDict, self.binNames, self.binCount = self.__readIntervals()
        elif isinstance(self.binData, tuple):
            (self.binDict, self.binNames, self.binCount = (
                self.__createIntervals())
        else:
            raise IOError('tuple or string must be supplied')
    
    def __createIntervals(self):
        # Unpack tuple
        chrFile, maxWidth, fixedWidth = self.binData
        # Create output data variables
        binDict = collections.OrderedDict()
        binNames = []
        count = 0
        # Extract chromosome name and length
        chrData = []
        with open(chrFile) as inFile:
            for line in inFile:
                chrName, chrLength = line.strip().split('\t')
                chrData.append((chrName, chrLength))
        # Loop through chromsome and generate bin data frame
        for chrName, chrLength in chrData:
            # Generate bin start and end for equal width bins
            if fixedWidth:
                # Calculate bin number
                remainder = float(int(chrLength) % maxWidth)
                # Calculate start and stop of first bin
                start = int(np.floor(remainder/2)) 
                end = start + maxWidth
                # Calculate start and stop of remianing bins
                binEnd = np.arange(end, int(chrLength) + 1, maxWidth,
                    dtype = np.uint32)
                binStart = (binEnd - maxWidth) + 1
            # Generate bin start and end sites for unequal width bins
            else:
                # Calculate number of bins
                binNo = np.ceil(float(chrLength) / maxWidth)
                # Calculate excess number of elements within bins
                excess = (binNo * maxWidth) - int(chrLength)
                # Calculate bin size 
                largeBin = int(maxWidth - np.floor(excess / binNo))
                smallBin = int(maxWidth - np.ceil(excess / binNo))
                # Calculate bin numbers
                smallBinNo = int(excess % binNo )
                largeBinNo = int(binNo - smallBinNo)
                # Generate bin width
                binWidth = np.array([largeBin] * largeBinNo + 
                    [smallBin] * smallBinNo)
                # Generate bins
                binEnd = np.cumsum(binWidth, dtype = np.uint32)
                binStart = (binEnd - binWidth) + 1
            # Generate arrays containing bin index and chromosome name
            binIndex = np.arange(count, count + len(binStart),
                dtype = np.uint32)
            count += len(binStart)
            binChr = np.array([chrName] *  len(binStart))
            # Store output data
            for bin in np.transpose([binChr, binStart, binEnd]):
                binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
            binDict[chrName] = {'start' : binStart, 'end' : binEnd,
                'index' : binIndex, 'count' : len(binStart)}
        # Return dataframe
        return(binDict, binNames, count)
    
    def __readIntervals(self):
        # Read in bed file to dataframe
        binChr = np.loadtxt(self.binData, usecols = (0,), dtype = str)
        binStart = np.loadtxt(self.binData, usecols = (1,), dtype = np.uint32)
        binEnd = np.loadtxt(self.binData, usecols = (2,), dtype = np.uint32)
        binIndex = np.arange(0, len(binStart), 1, dtype = np.uint32)
        if len(binChr) != len(binStart) or len(binStart) != len(binEnd):
            raise IOError('Invalid bed file supplied')
        # Create output data variables
        binDict = collections.OrderedDict()
        binNames = []
        count = len(binStart)
        # Create arrays for every chromosome
        for chrom in np.unique(binChr):
            # Extract arrays for chromosome
            indices = binChr == chrom
            chrChr = binChr[indices]
            chrStart = binStart[indices]
            chrEnd = binEnd[indices]
            chrIndex = binIndex[indices]
            # Check for overlapping bins
            binArray = np.empty(2 * sum(indices), dtype = np.uint32)
            binArray[0::2] = chrStart
            binArray[1::2] = chrEnd
            for i in xrange(len(binArray) - 1):
                if binArray[i + 1] <= binArray[i]:
                    raise IOError('Bed file contains overlapping bins')
            # Store output data
            for bin in np.transpose([chrChr, chrStart, chrEnd]):
                binNames.append('%s:%s-%s' %(bin[0], bin[1], bin[2]))
            binDict[chrom] = {'start' : chrStart, 'end' : chrEnd,
                'index' : chrIndex, 'count' : len(chrStart)}
        # Return data
        return(binDict, binNames, count)
    
    def __findBinIndex(self, chrom, position):
        # Convert position format
        position = np.uint32(position)
        # Extract bin location
        try:
            binLocation =  self.binDict[chrom]['end'].searchsorted(position)
        except KeyError:
            return('nochr')
        # Extract bin data
        try:
            binStart = self.binDict[chrom]['start'][binLocation]
            binIndex = self.binDict[chrom]['index'][binLocation]
        except IndexError:
            return('nobin')
        # Extract bin name/number
        if binStart <= position:
            return(binIndex)
        else:
            return('nobin')
    
    def writeIntervals(self, fileName):
        ''' Writes intervals to file.
        
        Args:
            fileName - Path to output file.
        
        '''
        with open(fileName, 'w') as outFile:
            for name in self.binNames:
                chrom, interval = name.split(':')
                start, end = interval.split('-')
                outFile.write('%s\t%s\t%s\n' %(chrom, start, end))
    
    def __matrixProcess(self, inputQueue, outPipe):
        # Create matrix and log array
        matrix = np.zeros((self.binCount, self.binCount),
            dtype = np.uint32)
        logData = np.zeros(4, dtype = np.uint32)
        # Extract data from pipe 
        for fragPair in iter(inputQueue.get, None):
            # Count total
            logData[0] += 1
            # Strip and process line
            fragPair = fragPair.strip().split('\t')
            # Extract location of fragends
            indices = []
            fragends = [fragPair[0:2], fragPair[3:5]]
            for fragChr, fragLoc in fragends:
                index = self.__findBinIndex(fragChr, fragLoc)
                if isinstance(index, np.uint32):
                    indices.append(index)
                else:
                    indices = index
                    break
            # Check that two bin indeces have been identified
            if isinstance(indices, list):
                # Count accepted ligations and data to array
                logData[3] += 1
                matrix[indices[0]][indices[1]] += 1
                matrix[indices[1]][indices[0]] += 1
            # Count incorrect indices
            elif indices == 'nochr':
                logData[1] += 1
            elif indices == 'nobin':
                logData[2] += 1
            else:
                raise ValueError('unrecognised bin index')
        outPipe.send((matrix, logData))
    
    def generateMatrix(self, fragendFile, threads=1):
        ''' Function to generate count matrices
        
        Args:
            fragendFile (str)- Path to fragendFile.
            threads (int)- Number of threads to use.
        
        Returns:
            finalMatrix - A 2D numpy array of interaction counts.
            finalLog - A dictionary of the 
        
        '''
        # Manage thread number
        if threads > 2:
            threads -= 1
        # Open input file
        if fragendFile.endswith('.gz'):
            sp = subprocess.Popen(["zcat", fragendFile],
                stdout = subprocess.PIPE)
            fh = cStringIO.StringIO(sp.communicate()[0])
        else:
            fh = open(fragendFile)
        # Create queue
        fragQueue = multiprocessing.Queue()
        # Start processes to count interactions
        processData = []
        for _ in range(threads):
            # Create pipes and process
            pipeReceive, pipeSend = multiprocessing.Pipe(False)
            process = multiprocessing.Process(
                target = self.__matrixProcess,
                args = (fragQueue, pipeSend)
            )
            process.start()
            pipeSend.close()
            # Strore process and pipe data
            processData.append((process,pipeReceive))
        # Add input data to queue
        for line in fh:
            fragQueue.put(line)
        # Add termination values to queue and close
        for _ in processData:
            fragQueue.put(None)
        fragQueue.close()
        # Extract data from processes and terminate
        for count, (process, pipe) in enumerate(processData):
            # Extract data and close process and pipes
            processMatrix, processLog = pipe.recv()
            process.join()
            pipe.close()
            # Add matrix count data
            if count:
                finalMatrix += processMatrix
                finalLog += processLog
            else:
                finalMatrix = processMatrix
                finalLog = processLog
        # Close input file
        if not fragendFile.endswith('.gz'):
            fh.close()
        # Return data
        return(finalMatrix, finalLog)

class normaliseCountMatrices(object):
    
    ''' Function to generate normalised matrices for multiple input count
    matrices. Function ensures that any bin with low counts in a single matrix
    is removed from all matrices. '''
    
    def __init__(self, matrixList, regionFile):
        ''' Initialise matrix normalisation object.
        
        Args:
            matrixList (list)- List of paths to count matries.
            regionFile (str)- Path to tabdelimited region file. Two acceptable
                formats:
                    1) A four column file listing chromosome name,
                    start position, end position and name of regions.
                    2) A two column file listing region name and a comma
                    seperated list of inidices.
        
        '''
        # Store supplied data
        self.matrixList = matrixList
        # Extract column number and extract indices
        with open (regionFile) as infile:
            colnumber = len(infile.next().split('\t'))
        if colnumber == 4:
            self.regionIndices = self.__generateRegionIndices(regionFile)
        elif colnumber == 2:
            self.regionIndices = self.__extractRegionIndices(regionFile)
    
    def __extractBinData(self):
        # Check bin names are identical across all files
        for count, infile in enumerate(self.matrixList):
            if not infile.endswith('countMatrix.gz'):
                raise IOError("Input files must end '.countMatrix.gz'")
            with gzip.open(infile, 'r') as openFile:
                header = openFile.next().strip().split('\t')
            if count:
                if not header == binNames:
                    raise IOError('Input files must have identical headers')
            else:
                binNames = header
        # Extract and store bin chromosome, start and end
        splitBin = [re.split('[:-]',x) for x in binNames]
        binChr = np.array([x[0] for x in splitBin])
        binStart = np.array([x[1] for x in splitBin], dtype = np.uint32)
        binEnd = np.array([x[2] for x in splitBin], dtype = np.uint32)
        # Check for overlapping bins
        for chrom in np.unique(binChr):
            # Extract indices for chromsome and check they are identical
            indices = np.where(binChr == chrom)[0]
            indexDiff = np.ediff1d(indices)
            if not np.all(indexDiff == 1):
                raise IOError('Count files contains unsorted bins')
            # Extract chromosomal bin start and stop sites and interpolate
            chrStart = binStart[indices]
            chrEnd = binEnd[indices]
            binArray = np.empty(2 * len(indices), dtype = np.uint32)
            binArray[0::2] = chrStart
            binArray[1::2] = chrEnd
            # Check bins are ordered and non-overlapping
            for i in xrange(len(binArray) - 1):
                if binArray[i + 1] <= binArray[i]:
                    raise IOError('Matrices contain anomalous bins')
        # Return data
        return(binChr, binStart, binEnd)
    
    def __generateRegionIndices(self, regionFile):
        ''' Function generates, and returns, a dictionary of indices for
        regions in self.regionFile.
        
        '''
        # Extract bin data
        binChr, binStart, binEnd = self.__extractBinData()
        # Create variable to store region indices and loop through file
        regionIndices = {}
        with open(regionFile, 'r') as infile:
            for line in infile:
                # Extract region data and find indices of matching bins
                chrom, start, end, region = line.strip().split('\t')
                acceptableBins = (
                    (binChr == chrom)
                    & (binStart >= np.uint32(start))
                    & (binEnd <= np.uint32(end))
                )
                indices = np.where(acceptableBins)[0]
                # Add region indices to dictionary
                if region in regionIndices:
                    regionIndices[region] = np.concatenate(
                        (regionIndices[region], indices))
                else:
                    regionIndices[region] = indices
        # Check region dictionary for absent or erronous regions
        for region in regionIndices:
            # Extract indices, sort, and update dictionary
            indices = regionIndices[region]
            indices.sort()
            regionIndices[region] = indices
            # Check for absent or duplicate indices
            if len(indices) == 0:
                raise IOError('{} has no bins'.format(region))
            if len(set(indices)) != len(indices):
                raise IOError('{} has overlapping segments'.format(region))
        # Return region index data
        return(regionIndices)
    
    def __extractRegionIndices(self, regionFile):
        ''' Function generates, and returns, a dictionary of indices for
        regions in self.regionFile.
        
        '''
        regionIndices = {}
        with open(regionFile, 'r') as infile:
            for line in infile:
                region, indices = line.strip().split('\t')
                regionIndices[region] = map(int, indices.split(','))
        return(regionIndices)
    
    def __colSumsMatrix(self, rmDiag, inQueue, outQueue):
        for matrix in iter(inQueue.get, None):
            # Read in counts and remove diagonal if required
            inMatrix = np.loadtxt(matrix, dtype = np.uint32,
                delimiter = '\t', skiprows = 1)
            if rmDiag:
                np.fill_diagonal(inMatrix, 0)
            # Check matrices are symetrical
            m, n = inMatrix.shape
            if m != n:
                raise ValueError('number of columns and rows must be equal')
            if not np.all(inMatrix == inMatrix.T):
                raise ValueError('matrix must be symetrical')
            # Create output data
            binSums = np.zeros(m, dtype=np.uint64)
            # Loop through regions and extract counts
            for indices in self.regionIndices.values():
                subMatrix = inMatrix[np.ix_(indices, indices)]
                binSums[indices] = subMatrix.sum(axis=1)
            # Return column data
            outQueue.put(binSums)
    
    def extractLowBins(
            self, minCount, rmDiag, threads
        ):
        ''' Function calclates low bins for each region '''
        # Check arguments
        if not isinstance(minCount, int) or minCount < 0:
            raise ValueError('minCount must be a non-negative integer')
        if not isinstance(rmDiag, bool):
            raise ValueError('rmDiag must be boolean')
        if not isinstance(threads, int) or threads < 1:
            raise ValueError('threads must a positive integer')
        # Create queues and processes
        inQueue = multiprocessing.Queue()
        outQueue = multiprocessing.Queue()
        processList = []
        for _ in range(threads):
            process = multiprocessing.Process(
                target = self.__colSumsMatrix,
                args = (rmDiag, inQueue, outQueue)
            )
            process.start()
            processList.append(process)
        # Add data to queues and extract results
        colSumList = []
        for matrix in self.matrixList:
            inQueue.put(matrix)
        for matrix in self.matrixList:
            colSumList.append(outQueue.get())
        # Close process and join processes and queues
        for _ in range(threads):
            inQueue.put(None)
        for process in processList:
            process.join()
        inQueue.close()
        outQueue.close()
        # Calculate low bins
        colSumArray = np.array(colSumList)
        lowBins = colSumArray.min(axis=0) < minCount
        # Generate and store low bins
        lowBinsDict = {}
        for region, indices in self.regionIndices.items():
            lowBinsDict[region] = lowBins[indices]
        return(lowBinsDict)
    
    def __saveSubMatrixProcess(
            self, lowBins, outDir, minCount, rmDiag, inQueue, outQueue
        ):
        # Extract matrices and output directory from index
        for matrix in iter(inQueue.get, None):
            # Extract sample name
            nameSearch = re.search('([^/]*).countMatrix.gz$', matrix)
            if not nameSearch:
                raise ValueError('Unrecognised matrix file name')
            sampleName = nameSearch.group(1)
            # Extract bin names
            with gzip.open(matrix) as inFile:
                inHeader = np.array(inFile.next().strip().split('\t'))
            # Create variable to store files and read in counts
            fileList = []
            counts = np.loadtxt(matrix, dtype = np.uint32, delimiter = '\t',
                skiprows = 1)
            # Open log file and loop through regions
            for region, indices in self.regionIndices.items():
                # Create output file names
                outSuffix = os.path.join(outDir, '{}.{}.{}'.format(
                    sampleName, region, minCount))
                if rmDiag:
                    outSuffix += '.noself'
                matrixFile = outSuffix + '.countMatrix.gz'
                logFile = outSuffix + '.log'
                # Extract sub matrix and extract metrics
                regionCounts = counts[np.ix_(indices, indices)]
                binNumber = len(indices)
                totalCounts = regionCounts.sum()
                # Count self ligation bins and remove if required
                selfLigationCount = np.diag(regionCounts).sum()
                if rmDiag:
                    np.fill_diagonal(regionCounts, 0)
                # Extract bins to be excluded from matrix
                regionLowBins = lowBins[region]
                if len(regionLowBins) != len(indices):
                    raise ValueError('Number of bins uncertain')
                regionCounts[regionLowBins,:] = 0
                regionCounts[:,regionLowBins] = 0
                lowBinCount = np.sum(regionLowBins)
                # Save output matrix 
                outHeader = '\t'.join(inHeader[indices])
                np.savetxt(matrixFile, regionCounts, fmt='%d',
                    delimiter='\t', header=outHeader, comments='')
                fileList.append(matrixFile)
                # Add data to log file
                with open(logFile, 'w') as outFile:
                    outFile.write('Grouped samples:\n')
                    for m in self.matrixList:
                        outFile.write('  {}\n'.format(m))
                    outFile.write('\nMatrix statistics:\n')
                    outFile.write('  total counts: {}\n'.format(totalCounts))
                    outFile.write('  bin count: {}\n'.format(binNumber))
                    outFile.write('  low bins: {}\n'.format(lowBinCount))
                    outFile.write('  self-ligation frequency: {:.2f}'\
                        '\n'.format(float(selfLigationCount) / totalCounts))
                    outFile.write('  remove self-ligation: {}\n'.format(
                        rmDiag))
            # Return file list
            outQueue.put(fileList)
    
    def saveSubMatrices(self, outDir, minCount, rmDiag, threads):
        # Check arguments
        if not isinstance(outDir, str) or not os.path.isdir(outDir):
            raise ValueError('directory {} not found'.format(outDir))
        if not isinstance(minCount, int) or minCount < 0:
            raise ValueError('minCount must be a non-negative integer')
        if not isinstance(rmDiag, bool):
            raise ValueError('rmDiag must be boolean')
        if not isinstance(threads, int) or threads < 1:
            raise ValueError('threads must a positive integer')
        # Extract low bins
        if minCount == 0:
            lowBins = {}
            for region, indices in self.regionIndices.items():
                lowBins[region] = np.zeros(len(indices), dtype=bool)
        else:
            lowBins = self.extractLowBins(
                minCount=minCount, rmDiag=rmDiag, threads=threads)
        # Create queue and processes:
        inQueue = multiprocessing.Queue()
        outQueue = multiprocessing.Queue()
        processList = []
        for _ in range(threads):
            process = multiprocessing.Process(
                target = self.__saveSubMatrixProcess,
                args = (lowBins, outDir, minCount, rmDiag, inQueue, outQueue)
            )
            process.start()
            processList.append(process)
        # Add data to queue
        for matrix in self.matrixList:
            inQueue.put(matrix)
        # Extract files from queue
        subMatrixList = []
        for matrix in self.matrixList:
            subMatrixList.extend(outQueue.get())
        # Cleanup queues and processes
        for _ in range(threads):
            inQueue.put(None)
        for process in processList:
            process.join()
        inQueue.close()
        outQueue.close()
        # Store sub matrices
        return(subMatrixList)
    
    def __calculateBias(
        self, matrix, max_iter, max_dev
    ):
        # Check matrix is symetrical
        m, n = matrix.shape
        if m != n:
            raise ValueError('Matrix must be square')
        if not np.all(matrix == matrix.T):
            raise ValueError('Matrix values must be symetrical')
        # Set initial variable for bias calculation
        bias = np.ones((m, 1), dtype=np.float64)
        nonzero = np.where(matrix.sum(axis=0) > 0)[0]
        # Reiteratively correct matrix
        for it in xrange(max_iter):
            # Extract current column sums
            adjColSums = (matrix / (bias * bias.T)).sum(axis=0)
            # Calculate deviation and halt iterations if acceptable
            dev = np.abs(adjColSums[nonzero] - 1).max()
            if dev <= max_dev:
                break
            # Adjust bias based on columns sums
            adjBias = np.sqrt(adjColSums)
            adjBias[adjBias == 0] = 1
            bias *= adjBias.reshape((m, 1))
        # Check normalisation has worked
        if dev > max_dev:
            raise ValueError('max deviation of {} is too high'.format(dev))
        # Return bias
        return(bias)
    
    def __iceNormalisationProcess(
            self, max_iter, max_dev, inQueue
        ):
        # Extract matrices and output directory from index
        for matrix in iter(inQueue.get, None):
            # Create output files
            outPrefix = matrix[:-15]
            biasFile = outPrefix + '.bias.gz'
            normMatrix = outPrefix + '.normMatrix.gz'
            # Extract header
            with gzip.open(matrix) as inFile:
                header = inFile.next().strip()
            # Open input matrix as integer array and check symetry
            inMatrix = np.loadtxt(matrix, dtype=np.float64, skiprows=1,
                delimiter='\t')
            bias = self.__calculateBias(inMatrix, max_iter, max_dev)
            np.savetxt(biasFile, bias,  delimiter='\t')
            # Generate normalised matrix and save data
            adjMatrix = inMatrix / (bias * bias.T)
            adjMatrix = np.round(adjMatrix, 9)
            np.savetxt(normMatrix, adjMatrix, '%.9f', delimiter='\t',
                header=header, comments='')
    
    def iceNormalisation(
            self, outDir, minCount, rmDiag, threads = 1, max_iter = 10000,
            max_dev = 1e-12
        ):
        ''' Function normalise count matrices.
        
        Args:
            outdir (str)- Full path to output directory.
            minCount (int)- Minimum number of interactions counts in bin for
                it to be included in the normalisation process.
            rmDiag (bool)- Whether to remove self-ligaton counts prior to
                nomralisation.
            threads (int)- Number of threads to use in normalisation process.
            max_iter (int)- Maximum number of iterations for normalisation.
            max_dev (float)- Deviation below which normalisation ends.
        
        '''
        # Check arguments
        if not os.path.isdir(outDir):
            raise ValueError('Output directory could not be found')
        if not isinstance(minCount, int) or minCount < 0:
            raise ValueError('minCount must be a non-negative integer')
        if not isinstance(rmDiag, bool):
            raise ValueError('rmDiag must be boolean')
        if not isinstance(threads, int) or threads < 1:
            raise ValueError('threads must a positive integer')
        if not isinstance(max_iter, int) or max_iter < 100:
            raise ValueError('iterations must a positive integer')
        # Create sub matrices
        subMatrices = self.saveSubMatrices(
            outDir=outDir, minCount=minCount, rmDiag=rmDiag, threads=threads)
        # Create queue and processes:
        inQueue = multiprocessing.Queue()
        processList = []
        for _ in range(threads):
            process = multiprocessing.Process(
                target = self.__iceNormalisationProcess,
                args = (max_iter, max_dev, inQueue)
            )
            process.start()
            processList.append(process)
        # Add data to queue
        for matrix in subMatrices:
            inQueue.put(matrix)
        # Cleanup queues and processes
        for _ in range(threads):
            inQueue.put(None)
        for process in processList:
            process.join()
        inQueue.close()
