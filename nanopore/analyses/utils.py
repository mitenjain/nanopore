import pysam
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, cigarReadFromString, PairwiseAlignment, system, fastaWrite, fastqWrite, cigarRead, logger, nameValue
import os
import sys
from cactus.bar import cactus_expectationMaximisation

class AlignedPair:
    """Represents an aligned pair of positions.
    """
    def __init__(self, refPos, refSeq, readPos, isReversed, readSeq, pPair):
        assert refPos >= 0 and refPos < len(refSeq)
        self.refPos = refPos
        self.refSeq = refSeq
        assert readPos >= 0 and readPos < len(readSeq)
        self.readPos = readPos
        self.isReversed = isReversed
        self.readSeq = readSeq
        self.pPair = pPair #Pointer to the previous aligned pair
        
    def isMatch(self):
        return self.getRefBase().upper() == self.getReadBase().upper() and self.getRefBase().upper() in "ACTG"
    
    def isMismatch(self):
        return self.getRefBase().upper() != self.getReadBase().upper() and self.getRefBase().upper() in "ACTG" and self.getReadBase().upper() in "ACTG"
    
    def getRefBase(self):
        return self.refSeq[self.refPos]
    
    def getReadBase(self):
        if self.isReversed:
            return reverseComplement(self.readSeq[self.readPos]) 
        return self.readSeq[self.readPos]
    
    def getSignedReadPos(self):
        if self.isReversed:
            return -self.readPos
        return self.readPos
    
    def getPrecedingReadInsertionLength(self, globalAlignment=False):
        if self.pPair == None:
            if globalAlignment:
                if self.isReversed:
                    assert len(self.readSeq) - self.readPos - 1 >= 0
                    return len(self.readSeq) - self.readPos - 1
                return self.readPos
            return 0
        return self._indelLength(self.readPos, self.pPair.readPos)
    
    def getPrecedingReadDeletionLength(self, globalAlignment=False):
        if self.pPair == None:
            if globalAlignment:
                return self.refPos
            return 0
        return self._indelLength(self.refPos, self.pPair.refPos)
    
    @staticmethod
    def _indelLength(pos, pPos):
        length = abs(pPos - pos) - 1
        assert length >= 0
        return length

    @staticmethod
    def iterator(alignedRead, refSeq, readSeq):
        """Generates aligned pairs from a pysam.AlignedRead object.
        """
        readOffset = getAbsoluteReadOffset(alignedRead, refSeq, readSeq)
        pPair = None
        assert len(alignedRead.seq) <= len(readSeq)
        for readPos, refPos in alignedRead.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedRead.pos and refPos < alignedRead.aend
                if refPos >= len(refSeq): #This is masking an (apparently minor?) one off error in the BWA sam files?
                    logger.critical("Detected an aligned reference position out of bounds! Reference length: %s, reference coordinate: %s" % (len(refSeq), refPos))
                    continue
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), alignedRead.is_reverse, readSeq, pPair)
                if aP.getReadBase().upper() != alignedRead.query[readPos].upper():
                    logger.critical("Detected a discrepancy between the absolute read sequence and the aligned read sequence. Bases: %s %s, read-position: %s, is reversed: %s, absolute read offset: %s, length absolute read sequence %s, length aligned read sequence %s, length aligned read sequence plus soft clipping %s, read name: %s, cigar string %s" % (aP.getReadBase().upper(), alignedRead.query[readPos].upper(), readPos, alignedRead.is_reverse, readOffset, len(readSeq), len(alignedRead.query), len(alignedRead.seq), alignedRead.qname, alignedRead.cigarstring))
                assert aP.getReadBase().upper() == alignedRead.query[readPos].upper()
                pPair = aP
                yield aP

def getAbsoluteReadOffset(alignedRead, refSeq, readSeq):
    """Gets the absolute starting coordinate of the first non-clipped position in the read.
    """
    if alignedRead.cigar[0][0] == 5: #Translate the read position to the original coordinates by removing hard clipping
        readOffset = alignedRead.cigar[0][1]
    else:
        readOffset = 0
    if alignedRead.is_reverse: #SEQ is reverse complemented
        readOffset = -(len(readSeq) - 1 - readOffset)
    readOffset += alignedRead.qstart #This removes any soft-clipping
    return readOffset

def getExonerateCigarFormatString(alignedRead, sam):
    """Gets a complete exonerate like cigar-string describing the sam line
    """
    for op, length in alignedRead.cigar:
        assert op in (0, 1, 2, 4, 5)
    translation = { 0:"M", 1:"I", 2:"D" }
    cigarString = " ".join([ "%s %i" % (translation[op], length) for op, length in alignedRead.cigar if op in translation ]) 
    completeCigarString = "cigar: %s %i %i + %s %i %i + 1 %s" % (
    alignedRead.qname, 0, alignedRead.qend - alignedRead.qstart, 
    sam.getrname(alignedRead.rname), alignedRead.pos, alignedRead.aend, cigarString)
    pA = cigarReadFromString(completeCigarString) #This checks it's an okay cigar
    assert sum([ op.length for op in pA.operationList if op.type == PairwiseAlignment.PAIRWISE_MATCH ]) == len([ readPos for readPos, refPos in alignedRead.aligned_pairs if readPos != None and refPos != None ])
    return completeCigarString

def samToBamFile(samInputFile, bamOutputFile):
    """Converts a sam file to a bam file (sorted)
    """
    samfile = pysam.Samfile(samInputFile, "r" )
    bamfile = pysam.Samfile(bamOutputFile, "wb", template=samfile)
    for line in samfile:
        bamfile.write(line)

    samfile.close()
    bamfile.close()
    
def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding fasta sequence
    """
    names = map(lambda x : x[0].split()[0], fastaRead(open(fastaFile, 'r')))
    assert len(names) == len(set(names)) #Check all the names are unique
    return dict(map(lambda x : (x[0].split()[0], x[1]), fastaRead(open(fastaFile, 'r')))) #Hash of names to sequences

def getFastqDictionary(fastqFile):
    """Returns a dictionary of the first words of fastq headers to their corresponding fastq sequence
    """
    names = map(lambda x : x[0].split()[0], fastqRead(open(fastqFile, 'r')))
    assert len(names) == len(set(names)) #Check all the names are unique
    return dict(map(lambda x : (x[0].split()[0], x[1]), fastqRead(open(fastqFile, 'r')))) #Hash of names to sequences

def makeFastaSequenceNamesUnique(inputFastaFile, outputFastaFile):
    """Makes a fasta file with unique names
    """
    names = set()
    fileHandle = open(outputFastaFile, 'w')
    for name, seq in fastaRead(open(inputFastaFile, 'r')):
        while name in names:
            logger.critical("Got a duplicate fasta sequence name: %s" % name)
            name += "i"
        names.add(name)
        fastaWrite(fileHandle, name, seq)
    fileHandle.close()
    return outputFastaFile

def makeFastqSequenceNamesUnique(inputFastqFile, outputFastqFile):
    """Makes a fastq file with unique names
    """
    names = set()
    fileHandle = open(outputFastqFile, 'w')
    for name, seq, quals in fastqRead(open(inputFastqFile, 'r')):
        while name in names:
            logger.critical("Got a duplicate fastq sequence name: %s" % name)
            name += "i"
        names.add(name)
        fastqWrite(fileHandle, name, seq, quals)
    fileHandle.close()
    return outputFastqFile

def normaliseQualValues(inputFastqFile, outputFastqFile):
    """Makes a fastq with valid qual values
    """
    fileHandle = open(outputFastqFile, 'w')
    for name, seq, quals in fastqRead(open(inputFastqFile, 'r')):
        if quals == None:
            quals = [33] * len(seq)
        fastqWrite(fileHandle, name, seq, quals)
    fileHandle.close()
    return outputFastqFile

def samIterator(sam):
    """Creates an iterator over the aligned reads in a sam file, filtering out
    any reads that have no reference alignment.
    """
    for aR in sam:
        if aR.rname != -1:
            yield aR

def mergeChainedAlignedReads(chainedAlignedReads, refSequence, readSequence):
    """Makes a global aligment for the given chained reads.
    From doc on building pysam line
    a = pysam.AlignedRead()
    a.qname = "read_28833_29006_6945"
    a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.flag = 99
    a.rname = 0
    a.pos = 32
    a.mapq = 20
    a.cigar = ( (0,10), (2,1), (0,25) )
    a.mrnm = 0
    a.mpos=199
    a.isize=167
    a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
    a.tags = ( ("NM", 1),
               ("RG", "L1") )
    """
    cAR = pysam.AlignedRead()
    aR = chainedAlignedReads[0]
    cAR.qname = aR.qname
    
    #Parameters we don't and therefore set properly
    #cAR.flag = aR.flag
    #cAR.mapq = aR.mapq
    #cAR.mrnm = 0
    #cAR.mpos=0
    #cAR.isize=0
    #cAR.qual = "<" * len(readSequence)
    #cAR.tags = aR.tags 
    cAR.rnext = -1
    cAR.pos = 0
    cAR.is_reverse = aR.is_reverse
    if cAR.is_reverse:
        cAR.seq = reverseComplement(readSequence)
    else:
        cAR.seq = readSequence
    cAR.rname = aR.rname
    cigarList = []
    pPos = 0
    if cAR.is_reverse: #Iterate from the other end of the sequence
        pQPos = -(len(readSequence)-1)
    else:
        pQPos = 0
        
    for aR in chainedAlignedReads:
        assert cAR.is_reverse == aR.is_reverse
        #Add a deletion representing the preceding unaligned reference positions
        assert aR.pos >= pPos
        if aR.pos > pPos:
            cigarList.append((2, aR.pos - pPos))
            pPos = aR.pos 
    
        #Add an insertion representing the preceding unaligned read positions
        qPos = getAbsoluteReadOffset(aR, refSequence, readSequence)
        assert qPos >= pQPos
        if qPos > pQPos:
            cigarList.append((1, qPos - pQPos)) 
            pQPos = qPos
        
        #Add the operations of the cigar, filtering hard and soft clipping
        for op, length in aR.cigar:
            assert op in (0, 1, 2, 4, 5)
            if op in (0, 1, 2):
                cigarList.append((op, length))
            if op in (0, 2): #Is match or deletion
                pPos += length
            if op in (0, 1): #Is match or insertion
                pQPos += length
        
    #Now add any trailing deletions/insertions
    assert pPos <= len(refSequence)
    if pPos < len(refSequence):
        cigarList.append((2, len(refSequence) - pPos))
    
    if cAR.is_reverse:
        assert pQPos <= 1
        if pQPos < 1:
            cigarList.append((1, -pQPos + 1))
    else:
        assert pQPos <= len(readSequence)
        if pQPos < len(readSequence):
            cigarList.append((1, len(readSequence) - pQPos))
    
    #Check coordinates
    #print cAR.is_reverse, sum([ length for op, length in cigarList if op in (0, 2)]),  len(refSequence), sum([ length for op, length in cigarList if op in (0, 1)]), len(readSequence), cAR.qname
    assert sum([ length for op, length in cigarList if op in (0, 2)]) == len(refSequence)
    assert sum([ length for op, length in cigarList if op in (0, 1)]) == len(readSequence)
    
    cAR.cigar = tuple(cigarList)
    
    return cAR

def chainFn(alignedReads, refSeq, readSeq, scoreFn=lambda alignedRead, refSeq, readSeq : len(list(AlignedPair.iterator(alignedRead, refSeq, readSeq)))):
    """Gets the highest scoring chain of alignments on either the forward or reverse strand. Score is (by default) number of aligned positions.
    """
    def getStartAndEndCoordinates(alignedRead):
        """Gets the start and end coordinates in both the reference and query
        """
        alignedPairs = list(AlignedPair.iterator(alignedRead, refSeq, readSeq))
        return alignedPairs[0].refPos, alignedPairs[0].getSignedReadPos(), alignedPairs[-1].refPos, alignedPairs[-1].getSignedReadPos()
    alignedReadToScores = dict([ (aR, scoreFn(aR, refSeq, readSeq)) for aR in alignedReads])
    alignedReadToCoordinates = dict([ (aR, getStartAndEndCoordinates(aR)) for aR in alignedReads])
    alignedReadPointers = {}
    
    #Currently uses sloppy quadratic algorithm to find highest chain
    alignedReads = sorted(alignedReads, key=lambda aR : alignedReadToCoordinates[aR][0]) #Sort by reference coordinate
    for i in xrange(len(alignedReads)):
        aR = alignedReads[i]
        rStart, qStart, rEnd, qEnd = alignedReadToCoordinates[aR]
        score = alignedReadToScores[aR]
        for j in xrange(i): #Look at earlier alignments in list
            aR2 = alignedReads[j]
            rStart2, qStart2, rEnd2, qEnd2 = alignedReadToCoordinates[aR2]
            assert rStart2 <= rStart
            if rStart > rEnd2 and qStart > qEnd2 and aR.is_reverse == aR2.is_reverse and \
            score + alignedReadToScores[aR2] > alignedReadToScores[aR]: #Conditions for a chain
                alignedReadToScores[aR] = score + alignedReadToScores[aR2]
                alignedReadPointers[aR] = aR2
    
    #Now find highest scoring alignment
    aR = sorted(alignedReads, key=lambda aR : alignedReadToScores[aR])[-1]
    
    #Construct chain of alignedReads
    chain = [ aR ]
    while aR in alignedReadPointers:
        aR = alignedReadPointers[aR]
        chain.append(aR)
    chain.reverse()
    
    return chain

def chainSamFile(samFile, outputSamFile, readFastqFile, referenceFastaFile, chainFn=chainFn):
    """Chains together the reads in the SAM file so that each read is covered by a single maximal alignment
    """
    sam = pysam.Samfile(samFile, "r" )
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    readSequences = getFastqDictionary(readFastqFile) #Hash of names to sequences
    readsToAlignedReads = {}
    for aR in samIterator(sam): #Iterate on the sam lines and put into buckets by read
        if aR.qname not in readSequences:
            raise RuntimeError("Aligned read name: %s not in read sequences names: %s" % (aR.qname, readSequences.keys()))
        key = (aR.qname,aR.rname)
        if key not in readsToAlignedReads:
            readsToAlignedReads[key] = []
        readsToAlignedReads[key].append(aR)
    #Now write out the sam file
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    
    #Chain together the reads
    chainedAlignedReads = []
    for readName, refID in readsToAlignedReads.keys():
        alignedReads = readsToAlignedReads[(readName, refID)]
        refSeq = refSequences[sam.getrname(refID)]
        readSeq = readSequences[readName]
        chainedAlignedReads.append(mergeChainedAlignedReads(chainFn(alignedReads, refSeq, readSeq), refSeq, readSeq))
    chainedAlignedReads.sort() #Sort by reference coordinate
    for cAR in chainedAlignedReads:
        outputSam.write(cAR)
    sam.close()
    outputSam.close()
    
def learnModelFromSamFileTargetFn(target, samFile, readFastqFile, referenceFastaFile, outputModel):
    """Does expectation maximisation on sam file to learn the hmm for the sam file.
    """
    #Convert the read file to fasta
    reads = os.path.join(target.getGlobalTempDir(), "temp.fa")
    fH = open(reads, 'w')
    for name, seq, quals in fastqRead(readFastqFile):
        fastaWrite(fH, name, seq)
    fH.close()
    
    #Get cigars file
    cigars = os.path.join(target.getGlobalTempDir(), "temp.cigar")
    fH = open(cigars, 'w')
    sam = pysam.Samfile(samFile, "r" )
    for aR in sam: #Iterate on the sam lines realigning them in parallel
        #Exonerate format Cigar string
        fH.write(getExonerateCigarFormatString(aR, sam) + "\n")
    fH.close()
    
    #Run cactus_expectationMaximisation
    options = cactus_expectationMaximisation.Options()
    options.modelType="fiveStateAsymmetric" #"threeStateAsymmetric"
    options.optionsToRealign="--diagonalExpansion=10 --splitMatrixBiggerThanThis=3000" 
    options.randomStart = True
    options.trials = 1
    options.iterations = 30
    options.numberOfAlignmentsPerJob=1
    #options.updateTheBand = True
    options.useDefaultModelAsStart = True
    #options.setJukesCantorStartingEmissions=0.1
    options.trainEmissions=True
    target.setFollowOnTargetFn(cactus_expectationMaximisation.expectationMaximisationTrials, args=(" ".join([reads, referenceFastaFile ]), cigars, outputModel, options))

def realignSamFileTargetFn(target, samFile, outputSamFile, readFastqFile, 
                           referenceFastaFile, gapGamma, hmmFileToTrain=None, chainFn=chainFn):
    """Chains and then realigns the resulting global alignments, using jobTree to do it in parallel on a cluster.
    Optionally runs expectation maximisation.
    """
    #Chain the sam file
    tempSamFile = os.path.join(target.getGlobalTempDir(), "temp.sam")
    chainSamFile(samFile, tempSamFile, readFastqFile, referenceFastaFile, chainFn)
    
    #If we do expectation maximisation we split here:
    if hmmFileToTrain != None:
        target.addChildTargetFn(learnModelFromSamFileTargetFn, args=(tempSamFile, readFastqFile, referenceFastaFile, hmmFileToTrain))
    
    target.setFollowOnTargetFn(realignSamFile2TargetFn, args=(tempSamFile, outputSamFile, readFastqFile, referenceFastaFile, hmmFileToTrain, gapGamma))

def realignSamFile2TargetFn(target, samFile, outputSamFile, readFastqFile, referenceFastaFile, hmmFile, gapGamma):
    #Load reference sequences
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    
    #Read through the SAM file
    sam = pysam.Samfile(samFile, "r" )
    tempCigarFiles = []
    for aR, index in zip(samIterator(sam), xrange(sys.maxint)): #Iterate on the sam lines realigning them in parallel
        #Exonerate format Cigar string
        cigarString = getExonerateCigarFormatString(aR, sam)
        
        #Temporary cigar file
        tempCigarFiles.append(os.path.join(target.getGlobalTempDir(), "rescoredCigar_%i.cig" % index))
        
        #Add a child target to do the alignment
        target.addChildTargetFn(realignCigarTargetFn, args=(getExonerateCigarFormatString(aR, sam), sam.getrname(aR.rname), refSequences[sam.getrname(aR.rname)], aR.qname, aR.query, tempCigarFiles[-1], hmmFile, gapGamma))
    
    target.setFollowOnTargetFn(realignSamFile3TargetFn, args=(samFile, outputSamFile, tempCigarFiles))
    #Finish up
    sam.close()
    
def realignCigarTargetFn(target, exonerateCigarString, referenceSequenceName, referenceSequence, querySequenceName, querySequence, outputCigarFile, hmmFile, gapGamma):
    #Temporary files
    tempRefFile = os.path.join(target.getGlobalTempDir(), "ref.fa")
    tempReadFile = os.path.join(target.getGlobalTempDir(), "read.fa")
    
    #Write the temporary files.
    fastaWrite(tempRefFile, referenceSequenceName, referenceSequence) 
    fastaWrite(tempReadFile, querySequenceName, querySequence)

    #Call to cactus_realign
    loadHmm = nameValue("loadHmm", hmmFile)
    system("echo %s | cactus_realign %s %s --diagonalExpansion=30 --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s > %s" % (exonerateCigarString, tempRefFile, tempReadFile, loadHmm, gapGamma, outputCigarFile))
    assert len([ pA for pA in cigarRead(open(outputCigarFile)) ]) > 0
    assert len([ pA for pA in cigarRead(open(outputCigarFile)) ]) == 1

def realignSamFile3TargetFn(target, samFile, outputSamFile, tempCigarFiles):
    #Setup input and output sam files
    sam = pysam.Samfile(samFile, "r" )
    
    #Replace the cigar lines with the realigned cigar lines
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    for aR, tempCigarFile in zip(samIterator(sam), tempCigarFiles): #Iterate on the sam lines realigning them in parallel
        #Load the cigar
        pA = [ i for i in cigarRead(open(tempCigarFile)) ][0]
        
        #Convert to sam line
        aR.cigar = tuple([ (op.type, op.length) for op in pA.operationList ])
        
        #Write out
        outputSam.write(aR)
    
    #Finish up
    sam.close()
    outputSam.close()

