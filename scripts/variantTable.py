import sys
from tex import *
import xml.etree.ElementTree as ET

fileHandle = open(sys.argv[1], 'w')

writeDocumentPreliminaries(fileHandle)
writePreliminaries(7, fileHandle)

recall, precision, fscore = {}, {}, {}
fileHandle2 = open(sys.argv[2], 'r')
mutFreq = [ 1, 5, 10, 20 ]
for i in mutFreq:
    l = fileHandle2.readline().split()
    readType = l[0]
    mapper = l[1]
    algorithm = l[2]
    recall[i] = l[4:19][1::3]
    precision[i] = l[19:34][1::3]
    fscore[i] = l[34:49][1::3]
fileHandle2.close()

def fn(x):
    return "%.2f" % (100*float(x))

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)
writeLine(7, 1, (("SNV detection using %s reads" % readType, 0, 6, 0, 0),), fileHandle)
writeLine(7, 2, (("Metric", 0, 0, 0, 1), 
                 ("Mut. Freq.", 1, 1, 0, 1), 
              ("Coverage", 2, 6, 0, 0), 
              ("10", 2, 2, 1, 1),
              ("30", 3, 3, 1, 1),
              ("60", 4, 4, 1, 1),
              ("120", 5, 5, 1, 1),
              ("ALL", 6, 6, 1, 1)), fileHandle)

def makeLine(word, hash):
    line = [ (word, 0, 0, 0, len(mutFreq)-1)]
    for i in xrange(len(mutFreq)):
        line.append((str(mutFreq[i]), 1, 1, i, i))
        for j in xrange(len(hash[mutFreq[i]])):
            line.append((fn(hash[mutFreq[i]][j]), j+2, j+2, i, i))
    writeLine(7, len(mutFreq), line, fileHandle)

makeLine("Recall", recall)
makeLine("Precision", precision)
makeLine("F-score", fscore)

writeEnd(fileHandle, "variantCallingTable", "Variant calling using %s reads with the %s mapping algorithm and the %s variant calling algorithm." % (readType, mapper, algorithm))
writeDocumentEnd(fileHandle)