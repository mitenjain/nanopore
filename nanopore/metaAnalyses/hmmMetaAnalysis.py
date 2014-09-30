import os
from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import xml.etree.cElementTree as ET
from jobTree.src.bioio import *
import itertools
import numpy as np

class HmmMetaAnalysis(AbstractMetaAnalysis):
    def run(self):
        #Call base method to do some logging
        AbstractMetaAnalysis.run(self) 
        bases = "ACGT"
        readTypes = set([ readType for readFastqFile, readType in self.readFastqFiles ])
        for readType in readTypes:
            #This is the info we collate
            transitions = {}
            insertEmissions = { "A":[], 'C':[], 'G':[], 'T':[] }
            deleteEmissions = { "A":[], 'C':[], 'G':[], 'T':[] }
            substitutionEmissions = dict(zip(itertools.product(bases, bases), map(lambda i : [], range(len(bases)**2))))
            
            for referenceFastaFile in self.referenceFastaFiles:
                for readFastqFile, readFileReadType in self.readFastqFiles:
                    if readFileReadType == readType:
                        for mapper in self.mappers:
                            analyses, resultsDir = self.experimentHash[((readFastqFile, readType), referenceFastaFile, mapper)]
                            if os.path.exists(os.path.join(resultsDir, "hmm.txt.xml")):
                                hmmsNode = ET.parse(os.path.join(resultsDir, "hmm.txt.xml")).getroot()
                                
                                #Aggregate transition expectations
                                for transition in hmmsNode.findall("transition"):
                                    if float(transition.attrib["avg"]) > 0.0:
                                        key= (transition.attrib["from"], transition.attrib["to"])
                                        if key not in transitions:
                                            transitions[key] = []
                                        transitions[key].append((float(transition.attrib["avg"]), float(transition.attrib["std"])))
                                
                                #Aggregate substitution expectations
                                #Plot match emission data
                                insertEmissions2 = dict(zip(bases, [0.0]*len(bases)))
                                deleteEmissions2 = dict(zip(bases, [0.0]*len(bases)))
                                for emission in hmmsNode.findall("emission"): 
                                    if emission.attrib["state"] == '0':
                                        substitutionEmissions[(emission.attrib["x"], emission.attrib["y"])].append((float(emission.attrib["avg"]), float(emission.attrib["std"])))
                                    elif emission.attrib["state"] == '1':
                                        deleteEmissions2[emission.attrib["x"]] += float(emission.attrib["avg"])
                                    elif emission.attrib["state"] == '2':
                                        insertEmissions2[emission.attrib["y"]] += float(emission.attrib["avg"])
                                for base in bases:
                                    insertEmissions[base].append(insertEmissions2[base])
                                    deleteEmissions[base].append(deleteEmissions2[base])
            
            #Write out a dot file representing the avg HMM transitions with std errors
            
            #Plot graphviz version of nanopore hmm, showing transitions and variances.
            fH = open(os.path.join(self.outputDir, "hmm_%s.dot" % readType), 'w')
            setupGraphFile(fH)
            #Make states
            addNodeToGraph("n0n", fH, "match")
            addNodeToGraph("n1n", fH, "short delete")
            addNodeToGraph("n2n", fH, "short insert")
            addNodeToGraph("n3n", fH, "long insert")
            addNodeToGraph("n4n", fH, "long delete")

            #Make edges with labelled transition probs.
            for fromState, toState in transitions:
                t = transitions[(fromState, toState)]
                addEdgeToGraph("n%sn" % fromState, 
                               "n%sn" % toState, fH, dir="arrow", style='""',
                                label="%.3f,%.3f" % (np.average(map(lambda x : x[0], t)), np.std(map(lambda x : x[0], t))))

            #Finish up
            finishGraphFile(fH)
            fH.close()       
            
            #Write out a matrix representing the substitutions with std errors, normalised by reference base
            matchEmissionsFile = os.path.join(self.outputDir, "matchEmissionsNormalisedByReference_%s.tsv" % readType)
            outf = open(matchEmissionsFile, "w")
            bases = "ACGT"
            outf.write("\t".join(bases) + "\n")
            for base in bases:
                d = sum(map(lambda x : np.average(substitutionEmissions[(base, x)][0]), bases))
                outf.write("\t".join([ base ] + map(lambda x : str(np.average(substitutionEmissions[(base, x)][0])/d), bases)) + "\n")
            outf.close()
            system("Rscript nanopore/analyses/substitution_plot.R %s %s %s" % (matchEmissionsFile, os.path.join(self.outputDir, "substitutionPlotNormalisedByReference_%s.pdf" % readType), "Avg. of ML substitution rates given the reference base"))
            
            #Write out the substitution errors. 
            #Write out a matrix representing the substitutions with std errors, normalised by reference base
            matchEmissionsFile = os.path.join(self.outputDir, "matchEmissionsUnnormalised_%s.tsv" % readType)
            outf = open(matchEmissionsFile, "w")
            bases = "ACGT"
            outf.write("\t".join(bases) + "\n")
            for base in bases:
                outf.write("\t".join([ base ] + map(lambda x : str(np.average(substitutionEmissions[(base, x)][0])), bases)) + "\n")
            outf.close()
            system("Rscript nanopore/analyses/substitution_plot.R %s %s %s" % (matchEmissionsFile, os.path.join(self.outputDir, "substitutionPlotUnnormalised_%s.pdf" % readType), "Avg. ML substitution estimates"))
            
            #Write out a matrix representing the substitutions with std errors, normalised by reference base
            matchEmissionsFile = os.path.join(self.outputDir, "matchEmissionsUnnormalisedStdErrors_%s.tsv" % readType)
            outf = open(matchEmissionsFile, "w")
            bases = "ACGT"
            outf.write("\t".join(bases) + "\n")
            for base in bases:
                outf.write("\t".join([ base ] + map(lambda x : str(np.average(substitutionEmissions[(base, x)][1])), bases)) + "\n")
            outf.close()
            system("Rscript nanopore/analyses/substitution_plot.R %s %s %s" % (matchEmissionsFile, os.path.join(self.outputDir, "substitutionPlotUnnormalisedStdErrors_%s.pdf" % readType), "Avg. ML substitution estimates"))
            
        