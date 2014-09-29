###Modify HMMs output from EM training to normalise for background nucleotide frequencies.
from cactus.bar.cactus_expectationMaximisation import Hmm, SYMBOL_NUMBER
import numpy as np
import sys

def main():
    print "ARGS", sys.argv
    #Load HMM
    hmm = Hmm.loadHmm(sys.argv[1])

    #Convert to matrix
    toMatrix = lambda e : map(lambda i : e[SYMBOL_NUMBER*i:SYMBOL_NUMBER*(i+1)], xrange(SYMBOL_NUMBER))
    fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)
    
    #Set indel emissions to all be normalised
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2

    #Normalise background emission frequencies, if requested to frequencies in given file
    gcContent = float(sys.argv[2])
    print "Got GC content", gcContent
    for state in range(hmm.stateNumber):
        n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)])
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), range(SYMBOL_NUMBER))) #Normalise

    #Modify match emissions by proposed substitution rate
    substitutionRate = float(sys.argv[3])
    print "Got substitution rate", substitutionRate
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), xrange(SYMBOL_NUMBER**2)))
    print "Got substitution matrix", n
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

    #Modify 

    #Write out HMM
    hmm.write(sys.argv[4])

if __name__ == '__main__':
    main()

