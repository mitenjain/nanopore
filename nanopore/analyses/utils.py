import pysam
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, cigarReadFromString, PairwiseAlignment

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
    
    def getPrecedingReadInsertionLength(self):
        if self.pPair == None:
            return 0
        return self._indelLength(self.readPos, self.pPair.readPos)
    
    def getPrecedingReadDeletionLength(self):
        if self.pPair == None:
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
        if alignedRead.cigar[0][0] == 5: #Translate the read position to the original coordinates by removing hard clipping
            readOffset = alignedRead.cigar[0][1]
        else:
            readOffset = 0
        if alignedRead.is_reverse: #SEQ is reverse complemented
            readOffset = -(len(readSeq) - 1 - readOffset)
        readOffset += alignedRead.qstart
        pPair = None
        for readPos, refPos in alignedRead.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedRead.pos and refPos < alignedRead.aend
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), alignedRead.is_reverse, readSeq, pPair)
                assert aP.getReadBase().upper() == alignedRead.query[readPos].upper()
                pPair = aP
                yield aP

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
    return dict([ (name.split()[0], seq) for name, seq in fastaRead(open(fastaFile, 'r'))]) #Hash of names to sequences

def getFastqDictionary(fastqFile):
    """Returns a dictionary of the first words of fastq headers to their corresponding fastq sequence
    """
    return dict([ (name.split()[0], seq) for name, seq, quals in fastqRead(open(fastqFile, 'r'))]) #Hash of names to sequences

