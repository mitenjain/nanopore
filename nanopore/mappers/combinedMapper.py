from nanopore.mappers.abstractMapper import AbstractMapper

from nanopore.mappers.last_params import LastParams
from nanopore.mappers.lastzParams import LastzParams
from nanopore.mappers.bwa_params import BwaParams
from nanopore.mappers.blasr_params import BlasrParams

import pysam
import os
from nanopore.analyses.utils import combineSamFiles

class CombinedMapper(AbstractMapper):
    def run(self):
        #Make child mapping job for each mapper
        mappers=[ LastParams, LastzParams, BwaParams, BlasrParams]
        tempSamFiles = []
        tempHmmFiles = []
        for mapper in mappers:
            tempSamFiles.append(os.path.join(self.getGlobalTempDir(), "mapping_%s.sam" % mapper.__name__))
            tempHmmFiles.append(os.path.join(self.getGlobalTempDir(), "mapping_%s.hmm" % mapper.__name__))
            self.addChildTarget(mapper(self.readFastqFile, self.readType, self.referenceFastaFile, tempSamFiles[-1], tempHmmFiles[-1]))
        #Now make a follow on to cat together the different alignments
        self.setFollowOnFn(combineSamFiles, args=(tempSamFiles[0], tempSamFiles[1:], self.outputSamFile))

class CombinedMapperChain2(AbstractMapper):
    def run(self):
        self.chainSamFile()

class CombinedMapperChain(AbstractMapper):
    def run(self, followOnMapper=CombinedMapperChain2):
        self.addChildTarget(CombinedMapper(self.readFastqFile, self.readType, self.referenceFastaFile, self.outputSamFile, self.emptyHmmFile))
        self.setFollowOnTarget(followOnMapper(self.readFastqFile, self.readType, self.referenceFastaFile, self.outputSamFile, self.emptyHmmFile))

class CombinedMapperRealign2(AbstractMapper):
    def run(self):
        self.realignSamFile()

class CombinedMapperRealign(CombinedMapperChain):
    def run(self):
        CombinedMapperChain.run(self, followOnMapper=CombinedMapperRealign2)

class CombinedMapperRealignEm2(AbstractMapper):
    def run(self):
        self.realignSamFile(doEm=True)

class CombinedMapperRealignEm(CombinedMapperChain):
    def run(self):
        CombinedMapperChain.run(self, followOnMapper=CombinedMapperRealignEm2)

class CombinedMapperRealignTrainedModel2(AbstractMapper):
    def run(self):
        self.realignSamFile(useTrainedModel=True)

class CombinedMapperRealignTrainedModel(CombinedMapperChain):
    def run(self):
        CombinedMapperChain.run(self, followOnMapper=CombinedMapperRealignTrainedModel2)

