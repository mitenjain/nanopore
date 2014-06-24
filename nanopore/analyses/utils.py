import pysam
from jobTree.src.bioio import reverseComplement, fastaRead

class AlignedPair:
    """Represents an aligned pair of positions.
    """
    def __init__(self, refPos, refSeq, readPos, isReversed, readSeq):
        assert refPos >= 0 and refPos < len(refSeq)
        self.refPos = refPos
        self.refSeq = refSeq
        assert readPos >= 0 and readPos < len(readSeq)
        self.readPos = readPos
        self.isReversed = isReversed
        self.readSeq = readSeq
    
    def getRefBase(self):
        return self.refSeq[self.refPos]
    
    def getReadBase(self):
        if self.isReversed:
            return reverseComplement(self.readSeq[self.readPos]) 
        return self.readSeq[self.readPos]

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
        for readPos, refPos in alignedRead.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedRead.pos and refPos < alignedRead.aend
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), alignedRead.is_reverse, readSeq)
                assert aP.getReadBase().upper() == alignedRead.seq[readPos].upper()
                yield aP
