from nanopore.mappers.lastz import Lastz
import pysam
import os

class LastzParams(Lastz):
    def run(self):
        Lastz.run(self, args="--hspthresh=1800 --gap=100,100")
        
class LastzParamsChain(LastzParams):
    def run(self):
        LastzParams.run(self)
        self.chainSamFile()
        
class LastzParamsRealign(LastzParams):
    def run(self):
        LastzParams.run(self)
        self.realignSamFile()
        
class LastzParamsRealignEm(LastzParams):
    def run(self):
        LastzParams.run(self)
        self.realignSamFile(doEm=True)

