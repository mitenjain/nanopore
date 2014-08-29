from nanopore.mappers.lastz import Lastz
import pysam
import os
from nanopore.analyses.utils import pathToBaseNanoporeDir


class LastzParams(Lastz):
    def run(self):
        #scoreFile = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "last_em_575_M13_2D_scores.txt")
        #Lastz.run(self, args="--hspthresh=1200 --gappedthresh=1500  --seed=match12 --scores=%s" % scoreFile)
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
        
class LastzParamsRealignTrainedModel(LastzParams):
    def run(self):
        LastzParams.run(self)
        self.realignSamFile(useTrainedModel=True)
