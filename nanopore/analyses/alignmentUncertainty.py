from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator
import os
import pysam
import numpy
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system, fastaWrite, cigarRead, PairwiseAlignment, cigarReadFromString

class AlignmentUncertainty(AbstractAnalysis):
    """Calculates stats on indels.
    """
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        
        #The data we collect
        avgPosteriorMatchProbabilityInCigar = []
        alignedPairsInCigar = []
        posteriorMatchProbabilities = []

        for aR in samIterator(sam): #Iterate on the sam lines
            #Exonerate format Cigar string
            cigarString = getExonerateCigarFormatString(aR, sam)
            
            #Temporary files
            tempCigarFile = os.path.join(self.getLocalTempDir(), "rescoredCigar.cig")
            tempRefFile = os.path.join(self.getLocalTempDir(), "ref.fa")
            tempReadFile = os.path.join(self.getLocalTempDir(), "read.fa")
            tempPosteriorProbsFile = os.path.join(self.getLocalTempDir(), "probs.tsv")
            
            #Write the temporary files.
            fastaWrite(tempRefFile, sam.getrname(aR.rname), refSequences[sam.getrname(aR.rname)]) 
            fastaWrite(tempReadFile, aR.qname, aR.query)
            
            #Call to cactus_realign
            system("echo %s | cactus_realign %s %s --rescoreByPosteriorProbIgnoringGaps --rescoreOriginalAlignment --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputPosteriorProbs=%s> %s" % (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, tempCigarFile))
            
            #Load the cigar and get the posterior prob
            assert len([ pA for pA in cigarRead(open(tempCigarFile)) ]) == 1
            pA = [ i for i in cigarRead(open(tempCigarFile)) ][0]
            avgPosteriorMatchProbabilityInCigar.append(pA.score)
            
            #Calculate the number of aligned pairs in the cigar
            alignedPairsInCigar.append(sum([ op.length for op in pA.operationList if op.type == PairwiseAlignment.PAIRWISE_MATCH ]))
            assert alignedPairsInCigar[-1] == len([ readPos for readPos, refPos in aR.aligned_pairs if readPos != None and refPos != None ])
            
            #Get the posterior probs
            #posteriorMatchProbabilities += [ float(line.split()[2]) for line in open(tempPosteriorProbsFile) ]
            
        sam.close()
        #Write out the substitution info
        node = ET.Element("alignmentUncertainty", { 
                "averagePosteriorMatchProbabilityPerRead":str(self.formatRatio(sum(avgPosteriorMatchProbabilityInCigar), len(avgPosteriorMatchProbabilityInCigar))),
                "averagePosteriorMatchProbability":str(self.formatRatio(float(sum([ avgMatchProb*alignedPairs for avgMatchProb, alignedPairs in zip(avgPosteriorMatchProbabilityInCigar, alignedPairsInCigar) ])),sum(alignedPairsInCigar))),
                "averagePosteriorMatchProbabilitesPerRead":",".join([ str(i) for i in avgPosteriorMatchProbabilityInCigar ]), 
                "alignedPairsInCigar":",".join([ str(i) for i in alignedPairsInCigar ]) })
        open(os.path.join(self.outputDir, "alignmentUncertainty.xml"), "w").write(prettyXml(node))
        if len(avgPosteriorMatchProbabilityInCigar) > 0:
            outf = open(os.path.join(self.getLocalTempDir(), "tmp_uncertainty"), "w")
            outf.write("\t".join([ str(i) for i in avgPosteriorMatchProbabilityInCigar ])); outf.write("\n")
            outf.close()
            system("Rscript nanopore/analyses/match_hist.R {} {}".format(os.path.join(self.getLocalTempDir(), "tmp_uncertainty"), os.path.join(self.outputDir, "posterior_prob_hist.pdf")))
        #Indicate everything is all done
        self.finish()

