from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator, pathToBaseNanoporeDir
import os
import pysam
import numpy
import math
import random
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system, fastaWrite, cigarRead, PairwiseAlignment, cigarReadFromString
from itertools import product
from cactus.bar.cactus_expectationMaximisation import Hmm

bases = "ACGT"

def getProb(subMatrix, start, end):
    return subMatrix[(start, end)]

def calcBasePosteriorProbs(baseObservations, refBase, 
                  evolutionarySubstitionMatrix, errorSubstutionMatrix):
    logBaseProbs = map(lambda missingBase : math.log(getProb(evolutionarySubstitionMatrix, refBase, missingBase)) + 
            reduce(lambda x, y : x + y, map(lambda observedBase : math.log(getProb(errorSubstutionMatrix, missingBase, observedBase))*baseObservations[observedBase], bases)), bases)
    totalLogProb = reduce(lambda x, y : x + math.log(1 + math.exp(y-x)), logBaseProbs)
    return dict(zip(bases, map(lambda logProb : math.exp(logProb - totalLogProb), logBaseProbs)))

def loadHmmErrorSubstitutionMatrix(hmmFile):
    hmm = Hmm.loadHmm(hmmFile)
    m = hmm.emissions[:len(bases)**2]
    m = map(lambda i : m[i] / sum(m[4*(i/4):4*(1 + i/4)]), range(len(m))) #Normalise m
    return dict(zip(product(bases, bases), m))

def getNullSubstitutionMatrix():
    return dict(zip(product(bases, bases), [1.0]*len(bases)**2))

def getJukesCantorTypeSubstitutionMatrix():
    return dict(zip(product(bases, bases), map(lambda x : 0.75 if x[0] == x[1] else (0.25/3), product(bases, bases))))

class MarginAlignSnpCaller(AbstractAnalysis):
    """Calculates stats on snp calling.
    """
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        
        node = ET.Element("marginAlignComparison")
        for hmmType in ("cactus", "trained_0",  "trained_20", "trained_40"): 
            for coverage in (1000000, 500, 120, 60, 30, 10): 
                for replicate in xrange(3 if coverage < 1000000 else 1): #Do replicates, unless coverage is all
                    sam = pysam.Samfile(self.samFile, "r" )
                    
                    #Trained hmm file to use.q
                    hmmFile0 = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "BLASR_DD_575_R7_M13_08_03_14_R72D_V1.3.1_hmm_0.txt")
                    hmmFile20 = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "BLASR_DD_575_R7_M13_08_03_14_R72D_V1.3.1_hmm_20.txt")
                    hmmFile40 = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "BLASR_DD_575_R7_M13_08_03_14_R72D_V1.3.1_hmm_40.txt")
              
                    #Get substitution matrices
                    nullSubstitionMatrix = getNullSubstitutionMatrix()
                    flatSubstitutionMatrix = getJukesCantorTypeSubstitutionMatrix()
                    hmmErrorSubstitutionMatrix = loadHmmErrorSubstitutionMatrix(hmmFile0)
                
                    #Load the held out snps
                    snpSet = {}
                    referenceAlignmentFile = self.referenceFastaFile + "_Index.txt"
                    if os.path.exists(referenceAlignmentFile):
                        seqsAndMutatedSeqs = getFastaDictionary(referenceAlignmentFile)
                        count = 0
                        for name in seqsAndMutatedSeqs:
                            if name in refSequences:
                                count += 1
                                trueSeq = seqsAndMutatedSeqs[name]
                                mutatedSeq = seqsAndMutatedSeqs[name + "_mutated"]
                                assert mutatedSeq == refSequences[name]
                                for i in xrange(len(trueSeq)):
                                    if trueSeq[i] != mutatedSeq[i]:
                                        snpSet[(name, i)] = trueSeq[i] 
                            else:
                                assert name.split("_")[-1] == "mutated"
                        assert count == len(refSequences.keys())
                    
                    #The data we collect
                    expectationsOfBasesAtEachPosition = {}
                    frequenciesOfAlignedBasesAtEachPosition = {}
                    
                    totalSampledReads = 0
                    totalAlignedPairs = 0
                    totalReadLength = 0
                    totalReferenceLength = sum(map(len, refSequences.values()))
                    
                    #Get a randomised ordering for the reads
                    reads = [ aR for aR in samIterator(sam) ]
                    random.shuffle(reads)
                    
                    for aR in reads: #Iterate on the sam lines
                        if totalReadLength/totalReferenceLength >= coverage: #Stop when coverage exceeds the quota
                            break
                        totalReadLength += len(readSequences[aR.qname])
                        totalSampledReads += 1
                        
                        #Temporary files
                        tempCigarFile = os.path.join(self.getLocalTempDir(), "rescoredCigar.cig")
                        tempRefFile = os.path.join(self.getLocalTempDir(), "ref.fa")
                        tempReadFile = os.path.join(self.getLocalTempDir(), "read.fa")
                        tempPosteriorProbsFile = os.path.join(self.getLocalTempDir(), "probs.tsv")
                        
                        #Ref name
                        refSeqName = sam.getrname(aR.rname)
                        
                        #Sequences
                        refSeq = refSequences[sam.getrname(aR.rname)]
                        
                        #Walk through the aligned pairs to collate the bases of aligned positions
                        for aP in AlignedPair.iterator(aR, refSeq, readSequences[aR.qname]): 
                            totalAlignedPairs += 1 #Record an aligned pair
                            key = (refSeqName, aP.refPos)
                            if key not in frequenciesOfAlignedBasesAtEachPosition:
                                frequenciesOfAlignedBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases))) 
                            readBase = aP.getReadBase() #readSeq[aP.readPos].upper() #Use the absolute read, ins
                            if readBase in bases:
                                frequenciesOfAlignedBasesAtEachPosition[key][readBase] += 1
                        
                        #Write the temporary files.
                        readSeq = aR.query #This excludes bases that were soft-clipped and is always of positive strand coordinates
                        fastaWrite(tempRefFile, refSeqName, refSeq) 
                        fastaWrite(tempReadFile, aR.qname, readSeq)
                        
                        #Exonerate format Cigar string, which is in readSeq coordinates (positive strand).
                        assert aR.pos == 0
                        assert aR.qstart == 0
                        assert aR.qend == len(readSeq)
                        assert aR.aend == len(refSeq)

                        cigarString = getExonerateCigarFormatString(aR, sam)
                        
                        #Call to cactus_realign
                        if hmmType == "trained_0":
                            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s > %s" % \
                                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, hmmFile0, tempCigarFile))
                        elif hmmType == "trained_20":
                            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s > %s" % \
                                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, hmmFile20, tempCigarFile))
                        elif hmmType == "trained_40":
                            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s > %s" % \
                                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, hmmFile40, tempCigarFile))
                        else:
                            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s > %s" % \
                                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, tempCigarFile))
                        
                        #Now collate the reference position expectations
                        for refPosition, readPosition, posteriorProb in map(lambda x : map(float, x.split()), open(tempPosteriorProbsFile, 'r')):
                            key = (refSeqName, int(refPosition))
                            if key not in expectationsOfBasesAtEachPosition:
                                expectationsOfBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases))) 
                            readBase = readSeq[int(readPosition)].upper()
                            if readBase in bases:
                                expectationsOfBasesAtEachPosition[key][readBase] += posteriorProb
                        
                        #Collate aligned positions from cigars
            
                    sam.close()
                    
                    totalHeldOut = len(snpSet)
                    totalNotHeldOut = totalReferenceLength - totalHeldOut
                    
                    class SnpCalls:
                        def __init__(self):
                            self.snpCalls = dict(zip(product("ACTGN", "ACTGN", "ACTGN"), [0.0]*(5**3)))
                            self.falsePositives = []
                            self.truePositives = []
                            self.falseNegatives = []
                        
                        @staticmethod
                        def bucket(calls):
                            calls = calls[:]
                            calls.sort()
                            buckets = [0.0]*101
                            for prob in calls: #Discretize
                                buckets[int(round(prob*100))] += 1
                            for i in xrange(len(buckets)-2, -1, -1): #Make cumulative
                                buckets[i] += buckets[i+1]
                            return buckets
                        
                        def getPrecisionByProbability(self):
                            tPs = self.bucket(map(lambda x : x[0], self.truePositives)) 
                            fPs = self.bucket(map(lambda x : x[0], self.truePositives))
                            return map(lambda i : float(tPs[i]) / (tPs[i] + fPs[i]) if tPs[i] + fPs[i] != 0 else 0, xrange(len(tPs)))
                        
                        def getFalsePositiveRatesByProbability(self):
                            return map(lambda i : i / (totalHeldOut + totalNotHeldOut), self.bucket(map(lambda x : x[0], self.falsePositives)))
                        
                        def getTruePositiveRatesByProbability(self):
                            return map(lambda i : i/totalHeldOut if totalHeldOut != 0 else 0, self.bucket(map(lambda x : x[0], self.truePositives)))
                        
                        def getTruePositiveLocations(self):
                            return map(lambda x : x[1], self.truePositives)
                        
                        def getFalsePositiveLocations(self):
                            return map(lambda x : x[1], self.falsePositives)
                        
                        def getFalseNegativeLocations(self):
                            return map(lambda x : x[0], self.falseNegatives)
            
                    #The different call sets
                    marginAlignMaxExpectedSnpCalls = SnpCalls()
                    marginAlignMaxExpectedSnpCallsBias95 = SnpCalls()
                    marginAlignMaxExpectedSnpCallsBias8 = SnpCalls()
                    marginAlignMaxExpectedSnpCallsBias6 = SnpCalls()
                    marginAlignMaxLikelihoodSnpCalls = SnpCalls()
                    marginAlignMaxLikelihoodSnpCallsBias95 = SnpCalls()
                    marginAlignMaxLikelihoodSnpCallsBias8 = SnpCalls()
                    marginAlignMaxLikelihoodSnpCallsBias6 = SnpCalls()
                    maxFrequencySnpCalls = SnpCalls()
                    maximumLikelihoodSnpCalls = SnpCalls()
                    
                    #Now calculate the calls
                    for refSeqName in refSequences:
                        refSeq = refSequences[refSeqName]
                        for refPosition in xrange(len(refSeq)):
                            mutatedRefBase = refSeq[refPosition].upper()
                            trueRefBase = (mutatedRefBase if not (refSeqName, refPosition) in snpSet else snpSet[(refSeqName, refPosition)]).upper()
                            key = (refSeqName, refPosition)
                            
                            #Get base calls
                            for errorSubstitutionMatrix, evolutionarySubstitutionMatrix, baseExpectations, biasFactor, snpCalls in \
                            ((flatSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 1.0, marginAlignMaxExpectedSnpCalls),
                             (flatSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.95, marginAlignMaxExpectedSnpCallsBias95),
                             (flatSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.8, marginAlignMaxExpectedSnpCallsBias8),
                             (flatSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.6, marginAlignMaxExpectedSnpCallsBias6),
                             (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 1.0, marginAlignMaxLikelihoodSnpCalls),
                             (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.95, marginAlignMaxLikelihoodSnpCallsBias95),
                             (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.8, marginAlignMaxLikelihoodSnpCallsBias8),
                             (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, 0.6, marginAlignMaxLikelihoodSnpCallsBias6),
                             (flatSubstitutionMatrix, nullSubstitionMatrix, frequenciesOfAlignedBasesAtEachPosition, 1.0, maxFrequencySnpCalls),
                             (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, frequenciesOfAlignedBasesAtEachPosition, 1.0, maximumLikelihoodSnpCalls)):
                                chosenBase = 'N'
                                if key in baseExpectations:
                                    #Get posterior likelihoods
                                    expectations = baseExpectations[key]
                                    totalExpectation = sum(expectations.values())
                                    if totalExpectation > 0.0: #expectationCallingThreshold:
                                        posteriorProbs = calcBasePosteriorProbs(dict(zip(bases, map(lambda x : float(expectations[x])/totalExpectation, bases))), trueRefBase, 
                                                               evolutionarySubstitutionMatrix, errorSubstitutionMatrix)
                                        posteriorProbs[mutatedRefBase] *= biasFactor
                                        maxPosteriorProb = max(posteriorProbs.values())
                                        chosenBase = random.choice([ base for base in posteriorProbs if posteriorProbs[base] == maxPosteriorProb ]).upper() #Very naive way to call the base
                                        if trueRefBase != mutatedRefBase:
                                            if trueRefBase == chosenBase:
                                                snpCalls.truePositives.append((maxPosteriorProb, refPosition)) #True positive
                                        elif chosenBase != mutatedRefBase: #This is a false positive as does not match either
                                            snpCalls.falsePositives.append((maxPosteriorProb, refPosition)) #False positive
                                    if trueRefBase != mutatedRefBase and trueRefBase != chosenBase:
                                         snpCalls.falseNegatives.append((refPosition, trueRefBase, mutatedRefBase, [ posteriorProbs[base] for base in "ACGT" ])) #False negative
                                        
                                #Add to margin-align max expected snp calls - "N" is a no-call.
                                snpCalls.snpCalls[(trueRefBase, mutatedRefBase, chosenBase)] += 1
                            
                    
                    for snpCalls, tagName in ((marginAlignMaxExpectedSnpCalls, "marginAlignMaxExpectedSnpCalls"), 
                                              (marginAlignMaxExpectedSnpCallsBias95, "marginAlignMaxExpectedSnpCallsBias95"), 
                                              (marginAlignMaxExpectedSnpCallsBias8, "marginAlignMaxExpectedSnpCallsBias8"), 
                                              (marginAlignMaxExpectedSnpCallsBias6, "marginAlignMaxExpectedSnpCallsBias6"), 
                                              (marginAlignMaxLikelihoodSnpCallsBias95, "marginAlignMaxLikelihoodSnpCallsBias95"), 
                                              (marginAlignMaxLikelihoodSnpCallsBias8, "marginAlignMaxLikelihoodSnpCallsBias8"), 
                                              (marginAlignMaxLikelihoodSnpCallsBias6, "marginAlignMaxLikelihoodSnpCallsBias6"), 
                                              (marginAlignMaxLikelihoodSnpCalls, "marginAlignMaxLikelihoodSnpCalls"),
                                              (maxFrequencySnpCalls, "maxFrequencySnpCalls"),
                                              (maximumLikelihoodSnpCalls, "maximumLikelihoodSnpCalls")):
                        fraction = lambda selectionFn : sum([ snpCalls.snpCalls[(true, mut, call)] for true, mut, call in snpCalls.snpCalls if selectionFn(true, mut, call) ])
                        
                        #Total hold outs
                        totalHeldOut = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true != mut)
                        #Total non-held out 
                        totalNonHeldOut = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true == mut)
                        #How many of held out cases do we predict correctly?
                        totalHeldOutCallsTrue = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true != mut and true == call)
                        #How many of held out cases do we predict wrongly and choose the reference?
                        totalHeldOutCallsFalseAndReference = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true != mut and mut == call)
                        #How many of held out cases do we predict wrongly and not choose the reference 
                        totalHeldOutCallsFalseAndNonReference = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and call != 'N' and true != mut and mut != call and true != call)
                        #How many of held out cases do we not call?
                        totalHeldOutNotCalled = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true != mut and call == 'N')
                        #How many of non-held out cases do we predict correctly?
                        totalNonHeldOutCallsTrue = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true == mut and true == call)
                        #How many of non-held out cases do we predict wrongly?
                        totalNonHeldOutCallsFalse = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and call != 'N' and true == mut and true != call)
                        #How many/proportion of non-held out cases do we not call?
                        totalNonHeldOutNotCalled = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true == mut and call == 'N')
                        
                        #Write out the substitution info
                        node2 = ET.SubElement(node, tagName + "_" + hmmType, {  
                                "coverage":str(coverage),
                                "actualCoverage":str(float(totalAlignedPairs)/totalReferenceLength),
                                "totalAlignedPairs":str(totalAlignedPairs),
                                "totalReferenceLength":str(totalReferenceLength),
                                "replicate":str(replicate),
                                "totalReads":str(len(reads)),
                                "avgSampledReadLength":str(float(totalReadLength)/totalSampledReads),
                                "totalSampledReads":str(totalSampledReads),
                                "totalHeldOut":str(totalHeldOut),
                                "totalNonHeldOut":str(totalNonHeldOut),
                                "totalHeldOutCallsTrue":str(totalHeldOutCallsTrue),
                                "totalFalsePositives":str(len(snpCalls.falsePositives)),
                                "totalHeldOutCallsFalseAndReference":str(totalHeldOutCallsFalseAndReference),
                                "totalHeldOutCallsFalseAndNonReference":str(totalHeldOutCallsFalseAndNonReference),
                                "totalHeldOutNotCalled":str(totalHeldOutNotCalled),
                                "totalNonHeldOutCallsTrue":str(totalNonHeldOutCallsTrue),
                                "totalNonHeldOutCallsFalse":str(totalNonHeldOutCallsFalse),
                                "totalNonHeldOutNotCalled":str(totalNonHeldOutNotCalled),
                                "falsePositiveRatesByProbability":" ".join(map(str, snpCalls.getFalsePositiveRatesByProbability())),
                                "truePositiveRatesByProbability":" ".join(map(str, snpCalls.getTruePositiveRatesByProbability())),
                                "precisionByProbability":" ".join(map(str, snpCalls.getPrecisionByProbability())),
                                "falsePositiveLocations":" ".join(map(str, snpCalls.getFalsePositiveLocations())),
                                "falseNegativeLocations":" ".join(map(str, snpCalls.getFalseNegativeLocations())),
                                "truePositiveLocations":" ".join(map(str, snpCalls.getTruePositiveLocations())) })
                        for snpCall in snpCalls.snpCalls:
                            ET.SubElement(node2, "%s_%s_%s" % snpCall, { "total":str(snpCalls.snpCalls[snpCall])})
                        for refPosition, trueRefBase, mutatedRefBase, posteriorProbs in snpCalls.falseNegatives:
                            ET.SubElement(node2, "falseNegative_%s_%s" % (trueRefBase, mutatedRefBase), { "posteriorProbs":" ".join(map(str, posteriorProbs))})
                        for falseNegativeBase in bases:
                            for mutatedBase in bases:
                                posteriorProbsArray = [ posteriorProbs for refPosition, trueRefBase, mutatedRefBase, posteriorProbs in snpCalls.falseNegatives if (trueRefBase.upper() == falseNegativeBase.upper() and mutatedBase.upper() == mutatedRefBase.upper() ) ]
                                if len(posteriorProbsArray) > 0:
                                    summedProbs = reduce(lambda x, y : map(lambda i : x[i] + y[i], xrange(len(x))), posteriorProbsArray)
                                    summedProbs = map(lambda x : float(x)/sum(summedProbs), summedProbs)
                                    ET.SubElement(node2, "combinedFalseNegative_%s_%s" % (falseNegativeBase, mutatedBase), { "posteriorProbs":" ".join(map(str, summedProbs))})
                        
        open(os.path.join(self.outputDir, "marginaliseConsensus.xml"), "w").write(prettyXml(node))
        
        
        #Indicate everything is all done
        self.finish()

