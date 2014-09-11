from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, getExonerateCigarFormatString, samIterator, pathToBaseNanoporeDir
import os
import pysam
import numpy
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system, fastaWrite, cigarRead, PairwiseAlignment, cigarReadFromString

class SnpCalling(AbstractAnalysis):
    """Calculates stats on snp calling.
    """
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        
        #Load the held out snps
        snpSet = {}
        referenceAlignmentFile = self.referenceFastaFile + "_Index.txt"
        if os.path.exists(referenceAlignmentFile):
            #Get the true and mutated reference sequences
            lines = [ i for i in open(referenceAlignmentFile, 'r') ]
            assert len(lines) == 2
            trueRefSeq = lines[0]
            mutatedRefSeq = lines[1]
            assert len(trueRefSeq) == len(mutatedRefSeq)
            #Check there is only one reference sequence (This is just a check for now / we may change the way we do this later)
            assert len(refSequences) == 1
            assert refSequences.values()[0] == mutatedRefSeq #This is true while there are no indel variants
            for i in xrange(len(trueRefSeq)):
                if trueRefSeq[i] != mutatedRefSeq[i]:
                   snpSet[(refSequences.keys()[0], i)] = trueRefSeq[i] 
        
        #The data we collect
        posteriorProbsOfBasesAtEachPosition = {}
        
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
            readSeq = aR.query
            
            #Write the temporary files.
            fastaWrite(tempRefFile, refSeqName, refSeq) 
            fastaWrite(tempReadFile, aR.qname, readSeq)
            
            #Trained hmm file to use.
            hmmFile = os.path.join(pathToBaseNanoporeDir(), "nanopore", "mappers", "last_em_575_M13_2D_hmm.txt")
            
            #Call to cactus_realign
            system("echo %s | cactus_realign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s > %s" % \
                   (cigarString, tempRefFile, tempReadFile, tempPosteriorProbsFile, hmmFile, tempCigarFile))
            
            #Now collate the reference positions
            for refPosition, readPosition, posteriorProb in map(lambda x : map(float, x.split()), open(tempPosteriorProbsFile, 'r')):
                key = (refSeqName, int(refPosition))
                if key not in posteriorProbsOfBasesAtEachPosition:
                    posteriorProbsOfBasesAtEachPosition[key] = { 'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0 }
                readBase = readSeq[int(readPosition)]
                if readBase in "ACGTacgt":
                    posteriorProbsOfBasesAtEachPosition[key][readBase] += posteriorProb

        sam.close()
        
        snpCalls = dict(zip(product("ACTGN", "ACTGN", "ACTGN"), [0.0]*(5**3)))
        
        #Now calculate the consensus sequence
        for refSeqName in refSequences:
            refSeq = refSequences[refSeqName]
            for refPosition in xrange(len(refSeq)):
                mutatedRefBase = refSeq[refPosition].upper()
                trueRefBase = (mutatedRefBase if not (refSeqName, refPosition) in snpSet else refSeq).upper()
                key = (refSeqName, refPosition)
                
                #Get posterior call
                maxExpectationBase = 'N'
                if key in posteriorProbsOfBasesAtEachPosition:
                    expectations = posteriorProbsOfBasesAtEachPosition[(refSeqName, refPosition)]
                    totalExpectation = sum(expectations.values())
                    maxExpectation = max(expectations.values())
                    if totalExpectation > expectationCallingThreshold:
                        maxExpectationBase = random.choice([ base for base in expectations if expectations[base] == maxExpectation ]).upper() #Very naive way to call the base
                
                key = (trueRefBase, mutatedRefBase, maxExpectationBase)
                snpCalls[(trueRefBase, mutatedRefBase, maxExpectationBase)] += 1
        
        fraction = lambda selectionFn : sum([ snpCalls(true, mut, call) for true, mut, call in snpCalls if selectionFn(true, mut, call) ])
        
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
        totalHeldOutNotCalled = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true != mut and true == 'N')
        #How many of non-held out cases do we predict correctly?
        totalNonHeldOutCallsTrue = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true == mut and true == call)
        #How many/proportion of non-held out cases do we not call?
        totalNonHeldOutNotCalled = fraction(lambda true, mut, call : true != 'N' and mut != 'N' and true == mut and call == 'N')
        
        #Write out the substitution info
        node = ET.Element("marginaliseConsensus", {  
                "totalHeldOut":str(totalHeldOut),
                "totalNonHeldOut":str(totalNonHeldOut),
                "totalHeldOutCallsTrue":str(totalHeldOutCallsTrue),
                "totalHeldOutCallsFalseAndReference":str(totalHeldOutCallsFalseAndReference),
                "totalHeldOutCallsFalseAndNonReference":str(totalHeldOutCallsFalseAndNonReference),
                "totalHeldOutNotCalled":str(totalHeldOutNotCalled),
                "totalNonHeldOutCallsTrue":str(totalNonHeldOutCallsTrue),
                "totalNonHeldOutNotCalled":str(totalNonHeldOutNotCalled) })
        for snpCall in snpCalls:
            ET.SubElement(node, "%s_%s_%s" % snpCall, { "total":snpCalls[snpCall]})
            
        open(os.path.join(self.outputDir, "marginaliseConsensus.xml"), "w").write(prettyXml(node))
        
        #Indicate everything is all done
        self.finish()

