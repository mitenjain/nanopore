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