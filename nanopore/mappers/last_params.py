from nanopore.mappers.abstractMapper import AbstractMapper
from nanopore.mappers.last import Last
from sonLib.bioio import system, fastaRead, fastqRead, fastaWrite
import os

class LastParams(Last):
    def run(self):
        Last.run(self, params="-s 2 -T 0 -Q 0 -a 1")
        
class LastParamsChain(LastParams):
    def run(self):
        LastParams.run(self)
        self.chainSamFile()
        
class LastParamsRealign(LastParams):
    def run(self):
        LastParams.run(self)
        self.realignSamFile()

class LastParamsRealignEm(LastParams):
    def run(self):
        LastParams.run(self)
        self.realignSamFile(doEm=True, gapGamma=0.5, matchGamma=0.0)

class LastParamsRealignTrainedModel(LastParams):
    def run(self):
        LastParams.run(self)
        self.realignSamFile(useTrainedModel=True)

class LastParamsRealignTrainedModel20(LastParams):
    def run(self):
        LastParams.run(self)
        self.realignSamFile(useTrainedModel=True, trainedModelFile="BLASR_DD_575_R7_M13_08_03_14_R72D_V1.3.1_hmm_20.txt")
        
class LastParamsRealignTrainedModel40(LastParams):
    def run(self):
        LastParams.run(self)
        self.realignSamFile(useTrainedModel=True, trainedModelFile="BLASR_DD_575_R7_M13_08_03_14_R72D_V1.3.1_hmm_40.txt")
     
