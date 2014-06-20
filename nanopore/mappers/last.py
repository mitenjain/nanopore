from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class Last(AbstractMapper):
    def run(self):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because we don't want to have any crufty files created in the local temp dir.
        indexFile = os.path.join(self.getLocalTempDir(), "my-index") #Index file
        mafFile = os.path.join(self.getLocalTempDir(), "out.maf") #MAF file
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile)) #Copy across the ref file
        system("lastdb %s %s" % (indexFile, localReferenceFastaFile)) #Build the index
        system("lastal %s %s > %s" % (indexFile, self.readFastaFile, mafFile)) #Build the alignment
        system("maf-convert.py sam %s > %s" % (mafFile, self.outputSamFile)) #Now convert sam file