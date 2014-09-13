from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator, pathToBaseNanoporeDir
import os
import pysam
import numpy
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
    baseProbs = map(lambda missingBase : getProb(evolutionarySubstitionMatrix, refBase, missingBase) * 
            sum(map(lambda observedBase : baseObservations[observedBase] * getProb(errorSubstutionMatrix, missingBase, observedBase), bases)), bases)
    return dict(zip(bases, map(lambda prob : prob/sum(baseProbs), baseProbs)))

def loadHmmErrorSubstitutionMatrix(hmmFile):
    hmm = Hmm.loadHmm(hmmFile)
    m = hmm.emissions[:len(bases)**2]
    m = map(lambda i : m[i] / sum(m[4*(i/4):4*(1 + i/4)]), range(len(m))) #Normalise m
    return dict(zip(product(bases, bases), m)) 

def getNullSubstitutionMatrix():
    return dict(zip(product(bases, bases), [1.0]*len(bases)**2))

def getIdentitySubstitutionMatrix():
    return dict(zip(product(bases, bases), map(lambda x : 1.0 if x[0] == x[1] else 0.0, product(bases, bases))))

def invertSubstitutionMatrix(substitutionMatrix):
    return dict(map(lambda x : ((x[1], x[0]), substitutionMatrix[x]), substitutionMatrix.keys()))

class MarginAlignSnpCaller(AbstractAnalysis):
    """Calculates stats on snp calling.
    """
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        
        #Trained hmm file to use.q
        hmmFile = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "last_em_575_M13_2D_hmm.txt")
 
        #Get substitution matrices
        nullSubstitionMatrix = getNullSubstitutionMatrix()
        identitySubstitutionMatrix = getIdentitySubstitutionMatrix()
        hmmErrorSubstitutionMatrix = loadHmmErrorSubstitutionMatrix(hmmFile)
    
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
        
        for aR in samIterator(sam): #Iterate on the sam lines
            #Exonerate format Cigar string
            cigarString = getExonerateCigarFormatString(aR, sam)
            
            #Temporary files
            tempCigarFile = os.path.join(self.getLocalTempDir(), "rescoredCigar.cig")
            tempRefFile = os.path.join(self.getLocalTempDir(), "ref.fa")
            tempReadFile = os.path.join(self.getLocalTempDir(), "read.fa")
            tempPosteriorProbsFile = os.path.join(self.getLocalTempDir(), "probs.tsv")
            
            #Ref name
            refSeqName = sam.getrname(aR.rname)
            
            #Sequences
            refSeq = refSequences[sam.getrname(aR.rname)]
            readSeq = aR.query #This excludes bases that were soft-clipped
            
            #Walk through the aligned pairs to collate the bases of aligned positions
            for aP in AlignedPair.iterator(aR, refSeq, readSequences[aR.qname]): 
                key = (refSeqName, aP.refPos)
                if key not in frequenciesOfAlignedBasesAtEachPosition:
                    frequenciesOfAlignedBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases))) 
                readBase = readSeq[aP.readPos].upper() #Use the absolute read, ins
                if readBase in bases:
                    frequenciesOfAlignedBasesAtEachPosition[key][readBase] += 1
            
            #Write the temporary files.
            fastaWrite(tempRefFile, refSeqName, refSeq) 
            fastaWrite(tempReadFile, aR.qname, readSeq)
            
            #Call to cactus_realign
            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s > %s" % \
                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, hmmFile, tempCigarFile))
            
            #Now collate the reference position expectations
            for refPosition, readPosition, posteriorProb in map(lambda x : map(float, x.split()), open(tempPosteriorProbsFile, 'r')):
                key = (refSeqName, int(refPosition))
                if key not in expectationsOfBasesAtEachPosition:
                    expectationsOfBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases))) 
                readBase = readSeq[int(readPosition)].upper()
                if readBase in bases:
                    expectationsOfBasesAtEachPosition[key][readBase] += posteriorProb

        sam.close()
        
        class SnpCalls:
            def __init__(self):
                self.snpCalls = dict(zip(product("ACTGN", "ACTGN", "ACTGN"), [0.0]*(5**3)))
                self.falsePositives = []
                self.truePositives = []
            
            @staticmethod
            def bucket(calls):
                calls = calls[:]
                calls.sort()
                buckets = [0.0]*101
                for prob in calls: #Discretize
                    buckets[int(round(prob*100))] += 1
                for i in xrange(len(buckets)-2, -1, -1): #Make cumulative
                    buckets[i] += buckets[i+1]
                return map(lambda x : 0 if buckets[0] == 0 else x/buckets[0], buckets) #Return normalised buckets
            
            def getCumulativeFalsePositives(self):
                return self.bucket(self.falsePositives)
            
            def getCumulativeTruePositives(self):
                return self.bucket(self.falsePositives)

        #The different call sets
        marginAlignMaxExpectedSnpCalls = SnpCalls()
        marginAlignMaxLikelihoodSnpCalls = SnpCalls()
        marginAlignMaxLikelihoodSnpCallsInvertCheck = SnpCalls()
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
                for errorSubstitutionMatrix, evolutionarySubstitutionMatrix, baseExpectations, snpCalls in \
                ((identitySubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, marginAlignMaxExpectedSnpCalls),
                 (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, expectationsOfBasesAtEachPosition, marginAlignMaxLikelihoodSnpCalls),
                 (invertSubstitutionMatrix(hmmErrorSubstitutionMatrix), nullSubstitionMatrix, expectationsOfBasesAtEachPosition, marginAlignMaxLikelihoodSnpCallsInvertCheck),
                 (identitySubstitutionMatrix, nullSubstitionMatrix, frequenciesOfAlignedBasesAtEachPosition, maxFrequencySnpCalls),
                 (hmmErrorSubstitutionMatrix, nullSubstitionMatrix, frequenciesOfAlignedBasesAtEachPosition, maximumLikelihoodSnpCalls)):
                    chosenBase = 'N'
                    if key in baseExpectations:
                        #Get posterior likelihoods
                        expectations = baseExpectations[key]
                        totalExpectation = sum(expectations.values())
                        if totalExpectation > 0.0: #expectationCallingThreshold:
                            posteriorProbs = calcBasePosteriorProbs(dict(zip(bases, map(lambda x : float(expectations[x])/totalExpectation, bases))), trueRefBase, 
                                                   evolutionarySubstitutionMatrix, errorSubstitutionMatrix)
                            maxPosteriorProb = max(posteriorProbs.values())
                            chosenBase = random.choice([ base for base in posteriorProbs if posteriorProbs[base] == maxPosteriorProb ]).upper() #Very naive way to call the base
                            
                            if trueRefBase == chosenBase:
                                snpCalls.truePositives.append(maxPosteriorProb)
                            elif chosenBase != mutatedRefBase: #This is a false positive as does not match either
                                snpCalls.falsePositives.append(maxPosteriorProb)
                            
                    #Add to margin-align max expected snp calls - "N" is a no-call.
                    snpCalls.snpCalls[(trueRefBase, mutatedRefBase, chosenBase)] += 1
                
        
        node = ET.Element("marginAlignComparison")
        
        for snpCalls, tagName in ((marginAlignMaxExpectedSnpCalls, "marginAlignMaxExpectedSnpCalls"), 
                                  (marginAlignMaxLikelihoodSnpCalls, "marginAlignMaxLikelihoodSnpCalls"),
                                  (marginAlignMaxLikelihoodSnpCallsInvertCheck, "marginAlignMaxLikelihoodSnpCallsInvertCheck"),
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
            node2 = ET.SubElement(node, tagName, {  
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
                    "cumulativeFalsePositives":" ".join(map(str, snpCalls.getCumulativeFalsePositives())),
                    "cumulativeTruePositives":" ".join(map(str, snpCalls.getCumulativeTruePositives())) })
            for snpCall in snpCalls.snpCalls:
                ET.SubElement(node2, "%s_%s_%s" % snpCall, { "total":str(snpCalls.snpCalls[snpCall])})
            
        open(os.path.join(self.outputDir, "marginaliseConsensus.xml"), "w").write(prettyXml(node))
        
        ####Put in ROC curves here
        
        
        #Indicate everything is all done
        self.finish()

