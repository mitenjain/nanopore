import os
from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import xml.etree.cElementTree as ET
from jobTree.src.bioio import *

class Hmm(AbstractAnalysis):
    """Calculates stats on indels.
    """
    def run(self):
        #Call base method to do some logging
        AbstractAnalysis.run(self) 
        
        #Hmm file
        hmmFile = os.path.join(os.path.split(self.samFile)[0], "hmm.txt.xml")
        if os.path.exists(hmmFile):
            #Load the hmm
            hmmsNode = ET.parse(hmmFile).getroot()

            #Plot graphviz version of nanopore hmm, showing transitions and variances.
            fH = open(os.path.join(self.outputDir, "hmm.dot"), 'w')
            setupGraphFile(fH)
            #Make states
            addNodeToGraph("n0n", fH, "match")
            addNodeToGraph("n1n", fH, "short insert")
            addNodeToGraph("n2n", fH, "short delete")
            addNodeToGraph("n3n", fH, "long insert")
            addNodeToGraph("n4n", fH, "long delete")

            #Make edges with labelled transition probs.
            for transition in hmmsNode.findall("transition"):
                addEdgeToGraph("n%sn" % transition.attrib["from"], 
                               "n%sn" % transition.attrib["to"], 
                               fH, dir="arrow", style='""',
                               label="%.3f,%.3f" % (float(transition.attrib["avg"]), float(transition.attrib["std"])))

            #Finish up
            finishGraphFile(fH)
            fH.close()

            #Plot match emission data
            emissions = dict([ ((emission.attrib["x"], emission.attrib["y"]), emission.attrib["avg"]) \
                  for emission in hmmsNode.findall("emission") if emission.attrib["state"] == '0' ])
            
            matchEmissionsFile = os.path.join(self.outputDir, "subst.tsv")
            outf = open(matchEmissionsFile, "w")
            bases = "ACGT"
            outf.write("\t".join(bases) + "\n")
            for base in bases:
                outf.write("\t".join([ base] + map(lambda x : emissions[(base, x)], bases)) + "\n")
            outf.close()
            system("Rscript nanopore/analyses/emissions_plot.R %s %s" % (matchEmissionsFile, os.path.join(self.outputDir, "substitution_plot.pdf")))

            #Plot indel info
            #Get the sequences to contrast the neutral model.
            refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
            readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
            #Need to do plot of insert and deletion gap emissions

            #Plot convergence data
            for hmmNode in hmmsNode.findall("hmm"):
                runningLikelihoods = map(float, hmmNode.attrib["runningLikelihoods"].split())
                #Need to do plot of iterations vs. running likelihoods.
            
        self.finish() #Indicates the batch is done