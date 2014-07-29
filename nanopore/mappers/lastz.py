from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system, fastaRead
import pysam

class Lastz(AbstractMapper):
    def run(self):
        system("lastz %s %s --format=sam > %s" % (self.referenceFastaFile, self.readFastqFile, self.outputSamFile))
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