from nanopore.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system, fastaRead, fastqRead, fastaWrite
import os

class Last(AbstractMapper):
    def run(self):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because we don't want to have any crufty files created in the local temp dir.
        indexFile = os.path.join(self.getLocalTempDir(), "my-index") #Index file
        mafFile = os.path.join(self.getLocalTempDir(), "out.maf") #MAF file
        #Hack to make last work, creating SQ line
        fH = open(self.outputSamFile, 'w')
        for name, seq in fastaRead(open(self.referenceFastaFile, 'r')):
            fH.write("@SQ\tSN:%s\tLN:%s\n" % (name.split()[0], len(seq)))
        fH.close()
        
        #Make fasta file, as last fastq seems broken
        localReadFile = os.path.join(self.getLocalTempDir(), "reads.fa") #Index file
        fH = open(localReadFile, 'w')
        for name, seq, quals in fastqRead(self.readFastqFile):
            fastaWrite(fH, name, seq)
        fH.close()
        
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile)) #Copy across the ref file
        system("lastdb %s %s" % (indexFile, localReferenceFastaFile)) #Build the index
        system("lastal %s %s > %s" % (indexFile, localReadFile, mafFile)) #Build the alignment
        system("maf-convert.py sam %s >> %s" % (mafFile, self.outputSamFile)) #Now convert sam file

class LastChain(Last):
    def run(self):
        Last.run(self)
        self.chainSamFile()
        
class LastRealign(Last):
    def run(self):
        Last.run(self)
        self.realignSamFile()