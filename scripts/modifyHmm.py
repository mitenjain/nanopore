###Modify HMMs output from EM training to normalise for background nucleotide frequencies.
from cactus.bar.cactus_expectationMaximisation import Hmm, SYMBOL_NUMBER
import sys
from nanopore.analyses.utils import setHmmIndelEmissionsToBeFlat, normaliseHmmByReferenceGCContent, modifyHmmEmissionsByExpectedVariationRate, toMatrix
import numpy as np

def main():
    print "ARGS", sys.argv
    #Load HMM
    hmm = Hmm.loadHmm(sys.argv[1])

    #setHmmIndelEmissionsToBeFlat(hmm)

    #Normalise background emission frequencies, if requested to GC% given
    gcContent = float(sys.argv[2])
    print "Got GC content", gcContent
    normaliseHmmByReferenceGCContent(hmm, gcContent)
    
    #Modify match emissions by proposed substitution rate
    substitutionRate = float(sys.argv[3])
    print "Got substitution rate", substitutionRate
    modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate)
    
    for state in range(0, hmm.stateNumber):
        n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)])
        print "For state, ref frequencies", map(sum, n)
        print "For state, read frequencies", map(sum, np.transpose(n))
    
    #Write out HMM
    hmm.write(sys.argv[4])

if __name__ == '__main__':
    main()

