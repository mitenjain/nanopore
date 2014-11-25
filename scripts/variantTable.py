import sys
from tex import *
import xml.etree.ElementTree as ET
import os.path

fileHandle = open(sys.argv[1], 'w')

writeDocumentPreliminaries(fileHandle)

recall, precision, fscore = {}, {}, {}
fileHandle2 = open(sys.argv[2], 'r')

tableNumber=1

while 1:
    l = fileHandle2.readline()
    if l == '':
        break
    tokens = l.split()
    if len(tokens) == 0 or tokens[0] != "readType":
        continue
    
    writePreliminaries(6, fileHandle)
    mutFreq = [ 1, 5, 10, 20 ]
    for i in mutFreq:
        l = fileHandle2.readline().split()
        readType = l[0]
        mapper = l[1]
        algorithm = l[2]
        recall[i] = l[4:16][1::3]
        precision[i] = l[16:28][1::3]
        fscore[i] = l[28:40][1::3]
    
    def fn(x):
        return "%.2f" % (100*float(x))
    
    #writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)
    writeLine(6, 1, (("SNV detection using %s reads" % readType, 0, 5, 0, 0),), fileHandle)
    writeLine(6, 2, (("Metric", 0, 0, 0, 1), 
                     ("Mut. Freq.", 1, 1, 0, 1), 
                  ("Coverage", 2, 5, 0, 0), 
                  ("30", 2, 2, 1, 1),
                  ("60", 3, 3, 1, 1),
                  ("120", 4, 4, 1, 1),
                  ("ALL", 5, 5, 1, 1)), fileHandle)
    
    def makeLine(word, hash):
        line = [ (word, 0, 0, 0, len(mutFreq)-1)]
        for i in xrange(len(mutFreq)):
            line.append((str(mutFreq[i]), 1, 1, i, i))
            for j in xrange(len(hash[mutFreq[i]])):
                line.append((fn(hash[mutFreq[i]][j]), j+2, j+2, i, i))
        writeLine(6, len(mutFreq), line, fileHandle)
    
    makeLine("Recall", recall)
    makeLine("Precision", precision)
    makeLine("F-score", fscore)
    
    #Format the algorithm/mapper name
    if mapper == "BlasrParamsChain":
        mapperName = "tuned Blasr (run using the `-x pacbio' flags)"
    elif mapper == "LastParamsChain":
        mapperName = "tuned Last (run using the `-s 2 -T 0 -Q 0 -a 1' flags)"
    elif mapper == "BlasrParamsRealignTrainedModel40":
        mapperName = "trained Blasr (Blasr run using the `-x pacbio' flags, realignment done using the described trained EM model)"
    
    mapperDescription = [ ]
    if "marginAlign" in algorithm:
        mapperDescription.append("Variant calling was performed using posterior match probabilities to integrate over every possible read alignment to the mutated reference sequence, using the initial guide alignment to band the calculations.")
    else:
        mapperDescription.append("Variant calling was performed conditioned on the fixed input alignment.")
    if "ikelihood" in algorithm:
        mapperDescription.append("Variant calling used a trained substitution matrix to calculate the maximum likelihood base (see method description).")
    else:
        mapperDescription.append("Variant calling corresponds to choosing the maximum-frequency/expectation of a non-reference base.")
    if "cactus" in algorithm:
        mapperDescription.append("Cactus realignment done using its stock pair-HMM.")
    if "trained_0" in algorithm:
        mapperDescription.append("Posterior match probabilities calculated using the EM trained HMM model, without accounting for substitution differences between the given reference and true underlying reference.")
    if "trained_20" in algorithm:
        mapperDescription.append("Posterior match probabilities calculated using the EM trained HMM model, accounting for substitution differences between the mutated reference and true underlying reference, assuming 20\% divergence.")
    if "trained_40" in algorithm:
        mapperDescription.append("Posterior match probabilities calculated using the EM trained HMM model, accounting for substitution differences between the mutated reference and true underlying reference, assuming 40\% divergence.")
    mapperDescription.append("Variant calling results shown for a posterior base calling probability threshold that gives the optimal F-score.")
    
    mapperDescription.append("Mutation frequency is the approximate proportion of sites mutated in the reference to which reads where aligned, and for which variants were called.")
    mapperDescription.append("Coverage is the total length of reads sampled divided by the length of the reference. ALL corresponds to using all the reads for a given experiment.")
    mapperDescription.append("Results shown are across three replicate experiments, and, at each coverage value, three different samplings of the reads. Raw results are available in the supplementary spread-sheet.")
    
    writeEnd(fileHandle, "variantCallingTable%i" % tableNumber, "Variant calling on M13 using %s reads starting with the %s mapping algorithm. %s" % (readType, mapperName, " ".join(mapperDescription)))
    
    algorithmShrunk = "".join(algorithm.split("_"))
    
    caption = "Precision/recall curves showing variant calling performance for four different mutation frequencies: 1, 5, 10 and 20 percent. Variant calling performed using %s reads starting with the %s mapping algorithm. %s" % (readType, mapperName, " ".join(mapperDescription))
    
    ##Write image
    writeFigure(fileHandle, os.path.join(sys.argv[3], "plots/%s%s/%sROCcurves.pdf" % (readType, mapper, algorithmShrunk)), caption, "variantCallingFigure%i" % tableNumber, width=15)
    
    tableNumber += 1
    
    
    
writeDocumentEnd(fileHandle)
fileHandle2.close()
fileHandle.close()    
