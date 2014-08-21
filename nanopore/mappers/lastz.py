from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system, fastaRead
from nanopore.analyses.utils import normaliseQualValues
import pysam
import os

class Lastz(AbstractMapper):
    def run(self, args=""):
        tempFastqFile = os.path.join(self.getLocalTempDir(), "temp.fastq")
        normaliseQualValues(self.readFastqFile, tempFastqFile)
        system("lastz %s %s %s --format=sam > %s" % (self.referenceFastaFile, tempFastqFile, args, self.outputSamFile))
        try:
            pysam.Samfile(self.outputSamFile, "r" ).close()
        except ValueError:
            #Hack to make lastz work, creating SQ lines when no alignments are found
            fH = open(self.outputSamFile, 'a')
            for name, seq in fastaRead(open(self.referenceFastaFile, 'r')):
                fH.write("@SQ\tSN:%s\tLN:%s\n" % (name.split()[0], len(seq)))
            fH.close()
        
class LastzChain(Lastz):
    def run(self):
        Lastz.run(self)
        self.chainSamFile()
        
class LastzRealign(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile()
        
class LastzRealignEm(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(doEm=True)
