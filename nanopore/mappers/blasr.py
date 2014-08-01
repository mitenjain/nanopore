from nanopore.mappers.abstractMapper import AbstractMapper
from nanopore.analyses.utils import getFastqDictionary
from sonLib.bioio import system
import os
import pysam

class Blasr(AbstractMapper):
    def run(self):
        tempSamFile = os.path.join(self.getLocalTempDir(), "temp.sam")
        system("blasr %s %s -sam > %s" % (self.readFastqFile, self.referenceFastaFile, tempSamFile))
        #Blasr seems to corrupt the names of read sequences, so lets correct them.
        sam = pysam.Samfile(tempSamFile, "r" )
        outputSam = pysam.Samfile(self.outputSamFile, "wh", template=sam)
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        for aR in sam: #Iterate on the sam lines and put into buckets by read
            if aR.qname not in readSequences:
                newName = '/'.join(aR.qname.split('/')[:-1])
                if newName not in readSequences:
                    raise RuntimeError("Tried to deduce correct read name: %s, %s" % (newName, readSequences.keys()))
                aR.qname = newName
            outputSam.write(aR)
        outputSam.close()
        

class BlasrChain(Blasr):
    def run(self):
        Blasr.run(self)
        self.chainSamFile()

class BlasrRealign(Blasr):
    def run(self):
        Blasr.run(self)
        self.realignSamFile()
