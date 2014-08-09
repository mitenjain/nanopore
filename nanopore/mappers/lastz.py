from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system, fastaRead
from nanopore.analyses.utils import normaliseQualValues
import pysam
import os

class Lastz(AbstractMapper):
    def run(self):
        tempFastqFile = os.path.join(self.getLocalTempDir(), "temp.fastq")
        normaliseQualValues(self.readFastqFile, tempFastqFile)
        system("lastz %s %s --format=sam > %s" % (self.referenceFastaFile, tempFastqFile, self.outputSamFile))
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

class LastzRealign_GapGamma0(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(gapGamma=0.0)
        
class LastzRealign_GapGamma0_Em(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(doEm=True, gapGamma=0.0)
        
class LastzRealign_GapGamma2(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(gapGamma=0.2)
        
class LastzRealign_GapGamma2_Em(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(doEm=True, gapGamma=0.2)
        
class LastzRealign_GapGamma5(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(gapGamma=0.5)

class LastzRealign_GapGamma5_Em(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(doEm=True, gapGamma=0.5)

class LastzRealign_GapGamma9(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(gapGamma=0.9)

class LastzRealign_GapGamma9_Em(Lastz):
    def run(self):
        Lastz.run(self)
        self.realignSamFile(doEm=True, gapGamma=0.9)
