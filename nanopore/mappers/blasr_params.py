from nanopore.mappers.blasr import Blasr
from sonLib.bioio import system
import os

class BlasrParams(Blasr):
    def run(self):
        Blasr.run(self, args="-sdpTupleSize 8 -bestn 1 -m 0")
        #system("blasr %s %s -sdpTupleSize 8 -bestn 1 -clipping hard -nproc 8 -sam -out %s -m 0" % (self.readFastqFile, self.referenceFastaFile, self.outputSamFile))

class BlasrParamsChain(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.chainSamFile()

class BlasrParamsRealign(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.realignSamFile()
        
class BlasrParamsRealignEm(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.realignSamFile(doEm=True)

class BlasrParamsRealignTrainedModel(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.realignSamFile(useTrainedModel=True)

class BlasrParamsRealignTrainedModelAMAP(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.realignSamFile(useTrainedModel=True, gapGamma=0.5, matchGamma=0.0)

class BlasrParamsRealignTrainedModelMatchGamma5(BlasrParams):
    def run(self):
        BlasrParams.run(self)
        self.realignSamFile(useTrainedModel=True, gapGamma=0.0, matchGamma=0.5)