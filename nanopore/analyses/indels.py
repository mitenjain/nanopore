from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import AlignedPair, getFastaDictionary, getFastqDictionary, samIterator
import os
import pysam
import numpy
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system 

class IndelCounter():
    def __init__(self, refSeqName, refSeq, readSeqName, readSeq, alignedRead):
        self.readInsertionLengths = []
        self.readDeletionLengths = []
        self.blockLengths = []
        self.readSeqName = readSeqName
        self.readSeq = readSeq
        self.refSeqName = refSeqName
        self.refSeq = refSeq

        #Now add the aligned read
        blockLength = 0
        for aP in AlignedPair.iterator(alignedRead, self.refSeq, self.readSeq): 
            if aP.getPrecedingReadInsertionLength() > 0:
                self.readInsertionLengths.append(aP.getPrecedingReadInsertionLength())
            if aP.getPrecedingReadDeletionLength() > 0:
                self.readDeletionLengths.append(aP.getPrecedingReadDeletionLength())
            if aP.getPrecedingReadInsertionLength() > 0 or aP.getPrecedingReadDeletionLength() > 0:
                assert blockLength > 0
                self.blockLengths.append(blockLength)
                blockLength = 1
            else:
                blockLength += 1

    def getXML(self):
        return ET.Element("indels", { "refSeqName":self.refSeqName, 
                                      "refSeqLength":str(len(self.refSeq)),
                                     "readSeqName":self.readSeqName,
                                     "readSeqLength":str(len(self.readSeq)),
                                     "numberReadInsertions":str(len(self.readInsertionLengths)),
                                     "numberReadDeletions":str(len(self.readDeletionLengths)),
                                     "avgReadInsertionLength":str(numpy.average(self.readInsertionLengths)),
                                     "avgReadDeletionLength":str(numpy.average(self.readDeletionLengths)),
                                     "medianReadInsertionLength":str(numpy.median(self.readInsertionLengths)),
                                     "medianReadDeletionLength":str(numpy.median(self.readDeletionLengths)),
                                     "readInsertionLengths":" ".join([ str(i) for i in self.readInsertionLengths ]),
                                     "readDeletionLengths":" ".join([ str(i) for i in self.readDeletionLengths ]) })

def getAggregateIndelStats(indelCounters):
    """Calculates aggregate stats across a set of read alignments.
    """
    
    readInsertionLengths = reduce(lambda x, y : x + y, map(lambda ic : ic.readInsertionLengths, indelCounters))
    readDeletionLengths = reduce(lambda x, y : x + y, map(lambda ic : ic.readDeletionLengths, indelCounters))
    
    attribs = { "numberOfReadAlignments":str(len(indelCounters)), 
                "readInsertionLengths":" ".join(map(str, readInsertionLengths)),
                "readDeletionLengths":" ".join(map(str, readDeletionLengths)) }
    
    readSequenceLengths = map(lambda ic : len(ic.readSeq), indelCounters)
    
    numberReadInsertions = map(lambda ic : len(ic.readInsertionLengths), indelCounters)
    numberReadDeletions = map(lambda ic : len(ic.readDeletionLengths), indelCounters)
    
    medianReadInsertionLengths = map(lambda ic : numpy.median(ic.readInsertionLengths), indelCounters)
    medianReadDeletionLengths = map(lambda ic : numpy.median(ic.readDeletionLengths), indelCounters)
    
    def stats(distribution):
        distribution = distribution[:]
        distribution.sort()
        return distribution[0], numpy.average(distribution), numpy.median(distribution), distribution[-1], " ".join(map(str, distribution))
    
    for name, distribution in [("ReadSequenceLengths", readSequenceLengths),
                               ("NumberReadInsertions", numberReadInsertions),
                               ("NumberReadDeletions", numberReadDeletions),
                               ("MedianReadInsertionLengths", medianReadInsertionLengths),
                               ("MedianReadDeletionLengths", medianReadDeletionLengths) ]:
        for attribName, value in zip([ "min" + name, "avg" + name, "median" + name, "max" + name, "distribution" + name ], list(stats(distribution))):
            attribs[name] = str(value)
            
    parentNode = ET.Element("indels", attribs)
    for indelCounter in indelCounters:
        parentNode.append(indelCounter.getXML())
    return parentNode

class Indels(AbstractAnalysis):
    """Calculates stats on indels.
    """
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        refSequences = getFastaDictionary(self.referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(self.readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(self.samFile, "r" )
        indelCounters = map(lambda aR : IndelCounter(sam.getrname(aR.rname), refSequences[sam.getrname(aR.rname)], aR.qname, readSequences[aR.qname], aR), samIterator(sam)) #Iterate on the sam lines
        sam.close()
        #Write out the substitution info
        if len(indelCounters) > 0:
            indelXML = getAggregateIndelStats(indelCounters)
            open(os.path.join(self.outputDir, "indels.xml"), "w").write(prettyXml(indelXML))
            #tmp = open(os.path.join(self.outputDir, "tmp.tsv"), "w")
            tmp = open(os.path.join(self.getLocalTempDir(), "tmp.tsv"), "w")
            #build list of data as vectors
            data_list = []
            var = ["readInsertionLengths", "readDeletionLengths", "ReadSequenceLengths", "NumberReadInsertions", "NumberReadDeletions", "MedianReadInsertionLengths", "MedianReadDeletionLengths"]
            for x in var:
                data_list.append([x] + indelXML.attrib[x].split())
            #transpose this list so R doesn't take hours to load it using magic
            data_list = map(None, *data_list)
            for line in data_list:
                tmp.write("\t".join(map(str,line))); tmp.write("\n")
            tmp.close()
            #system("Rscript nanopore/analyses/indelPlots.R {} {}".format(os.path.join(self.outputDir, "tmp.tsv"), os.path.join(self.outputDir, "indel_plots.pdf")))
            system("Rscript nanopore/analyses/indelPlots.R {} {}".format(os.path.join(self.getLocalTempDir(), "tmp.tsv"), os.path.join(self.outputDir, "indel_plots.pdf")))
        
        self.finish() #Indicates the batch is done


            #Plots:
            ##Read insertion lengths plot: x-axis insertion legnth, y-axis: frequency (see map(int, indelXML.attrib["readInsertionLengths].split()))
            ##Read deletion lengths plot: x-axis deletion length, y-axis: frequency (see map(int,indelXML.attrib["readDeletionLengths].split()))
            ##Number of read insertions vs. read length (map(int, indelXML.attrib["distributionReadSequenceLengths"].split()) vs. map(int, indelXML.attrib["distributionNumberReadInsertions"].split())
            ##Number of read deletion vs. read length
            ##Distribution of median insertion lengths (map(int, indelXML.attrib["distributionMedianReadInsertionLengths"].split()))
            ##Distribution of deletion  lengths (map(int, indelXML.attrib["distributionMedianReadInsertionLengths"].split()))

