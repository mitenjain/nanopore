from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import pysam
import numpy
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system 

class IndelCounter():
    def __init__(self, readSeqName, refSeqName):
        self.readInsertionLengths = []
        self.readDeletionLengths = []
        self.blockLengths = []
        self.readSeqName = readSeqName
        self.refSeqName = refSeqName
    
    def addReadAlignment(self, alignedRead, refSeq, readSeq):
        blockLength = 0
        for aP in AlignedPair.iterator(alignedRead, refSeq, readSeq): 
            if aP.getPrecedingReadInsertionLength() > 0:
                self.readInsertionLengths.append(aP.getPrecedingReadInsertionLength())
            if aP.getPrecedingReadDeletionLength() > 0:
                self.readDeletionLengths.append(aP.getPrecedingReadDeletionLength())
            if aP.getPrecedingReadInsertionLength() > 0 or aP.getPrecedingReadDeletionLength() > 0:
                assert blockLength > 0
                self.blockLengths.append(blockLength)
                blockLength = 1
            else:
                blockLength += 1
    
    def getXML(self):
        return ET.Element("indels", { "refSeqName":self.refSeqName, 
                                     "readSeqName":self.readSeqName,
                                     "totalReadInsertions":str(len(self.readInsertionLengths)),
                                     "totalReadDeletions":str(len(self.readDeletionLengths)),
                                     "avgReadInsertionLength":str(numpy.average(self.readInsertionLengths)),
                                     "avgReadDeletionLength":str(numpy.average(self.readDeletionLengths)),
                                     "medianReadInsertionLength":str(numpy.median(self.readInsertionLengths)),
                                     "medianReadDeletionLength":str(numpy.median(self.readDeletionLengths)),
                                     "readInsertionLengths":" ".join([ str(i) for i in self.readInsertionLengths ]),
                                     "readDeletionLengths":" ".join([ str(i) for i in self.readDeletionLengths ]) })
    
class Indels(AbstractAnalysis):
    """Calculates stats on indels.
    """
    def run(self):
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        overallIndelCounter = IndelCounter("overall", "overall")
        for aR in samIterator(sam): #Iterate on the sam lines
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = readSequences[aR.qname]
            overallIndelCounter.addReadAlignment(aR, refSeq, readSeq)
        sam.close()
        #Write out the substitution info
        open(os.path.join(self.outputDir, "indels.xml"), "w").write(prettyXml(overallIndelCounter.getXML()))
        stats = dict(ET.parse(os.path.join(self.outputDir, "indels.xml")).findall(".")[0].items())
        outf = open(os.path.join(self.outputDir, "stats.tsv"), "w")
        for key in stats:
            if key == "readInsertionLengths":
                open(os.path.join(self.getLocalTempDir(), "r_insert.txt"), "w").write(stats[key])
            elif key == "readDeletionLengths":
                open(os.path.join(self.getLocalTempDir(), "r_delete.txt"), "w").write(stats[key])
            else:
                outf.write("{}\t{}\n".format(key, str(stats[key])))
        outf.close()
        system("Rscript nanopore/analyses/indel_plot.R {} {} {} {}".format(os.path.join(self.outputDir, "stats.tsv"), os.path.join(self.getLocalTempDir(), "r_insert.txt"), os.path.join(self.getLocalTempDir(), "r_delete.txt"), os.path.join(self.outputDir, "indel_hist.pdf")))
