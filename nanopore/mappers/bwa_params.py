from nanopore.mappers.bwa import Bwa
from sonLib.bioio import system
import os

class BwaParams(Bwa):
    def run(self):
        Bwa.run(self, args="-x pacbio")

class BwaParamsChain(BwaParams):
    def run(self):
        BwaParams.run(self)
        self.chainSamFile()
        
class BwaParamsRealign(BwaParams):
    def run(self):
        BwaParams.run(self)
        self.realignSamFile()

class BwaParamsRealignEm(BwaParams):
    def run(self):
        BwaParams.run(self)
        self.realignSamFile(doEm=True)
        
class BwaParamsRealignTrainedModel(BwaParams):
    def run(self):
        BwaParams.run(self)
        self.realignSamFile(useTrainedModel=True)