#!/usr/bin/env python

import sys, os, pysam
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 
from jobTree.src.bioio import fastqRead, fastaRead, setLoggingFromOptions, logger, popenCatch, system
from itertools import izip, product
from collections import Counter, defaultdict

"""Used to BLAST results from the combined analysis result in order to find the set of reads that truly unmappable.
Reports the BLAST results as well as the raw hits, and also reports a fasta of the reads that did not map anywhere,
and generates a summary barplot."""

readTypes = ["2D", "template", "complement"]
combinedAnalyses = ["LastzParamsRealignEm","LastParamsRealignEm","BwaParamsRealignEm","BlasrParamsRealignEm"]#,"BwaChain","BlasrChain","LastChain","LastzChain"]

def parse_blast(blast_handle):
    """generator to yield blast results for each read, iterating over blast with outfmt="7 qseqid sseqid sscinames stitle"
    and max_target set to 1"""
    result = None
    for line in blast_handle:
        if "0 hits found" in line:
            yield (query, result)
        elif line.startswith("#") and "Query: " in line:
            query = line.split("Query: ")[-1].rstrip()
        elif result is None and not line.startswith("#"):
            result = line.strip().split("\t")[-3::]
            yield (query, result)
        elif result is not None and line.startswith("#"):
            result = None

def find_analyses(target, unmappedByReadType, outputDir):
    outfiles = dict()
    for readType in unmappedByReadType:
        outfiles[readType] = list()
        records = list()
        for (name, sequence), i in izip(unmappedByReadType[readType].iteritems(), xrange(len(unmappedByReadType[readType]))):
                records.append(">{}\n{}\n".format(name, sequence))
                if i % 10 == 1200 or i == len(unmappedByReadType[readType]) - 1:
                    tmpalign = os.path.join(target.getGlobalTempDir(), str(i) + ".txt")
                    outfiles[readType].append(tmpalign)
                    target.addChildTarget(Target.makeTargetFn(run_blast, args=(records, tmpalign)))
                    records = list()
    target.setFollowOnTargetFn(merge, args=(outfiles, outputDir))

def run_blast(target, records, tmpalign):
    query = "".join(records)
    p = popenCatch('blastn -outfmt "7 qseqid sseqid sscinames stitle" -db nt', stdinString=query)
    outf = open(tmpalign, "w"); outf.write(p); outf.close()


def merge(target, outfiles, outputDir):
    for readType in outfiles:
        outf = os.path.join(outputDir, readType + "_blast_out.txt")
        with open(outf, "w") as outfile:
            for f in outfiles[readType]:
                with open(f) as i:
                    outfile.write(i.read())
        outfile.close()

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    outputDir = "blast_combined/output/"

    if not os.path.exists(outputDir):
        logger.info("Output dir {} does not exist. Creating.")
        os.mkdir(outputDir)
    if len(os.listdir(outputDir)) > 0:
        logger.info("Output dir not empty.")

    #find all read fastq files, load into a dict by read type
    readFastqFiles = dict()
    for readType in readTypes:
        readFastqFiles[readType] = [os.path.join("../output/processedReadFastqFiles/", readType, x) for x in os.listdir(os.path.join("../output/processedReadFastqFiles/", readType)) if x.endswith(".fq") or x.endswith(".fastq")]
    
    #find all reference fasta files
    referenceFastaFiles = [x for x in os.listdir("../referenceFastaFiles") if x.endswith(".fasta") or x.endswith(".fa")]

    #find all sam files that were analyzed using combinedAnalyses
    samFiles = {}
    for readType in readTypes:
        samFiles[readType] = [(readFastqFile, os.path.join("../output", "analysis_" + readType, "experiment_" + os.path.basename(readFastqFile) + "_" + referenceFastaFile + "_" + analysis, "mapping.sam")) for readFastqFile, referenceFastaFile, analysis in product(readFastqFiles[readType], referenceFastaFiles, combinedAnalyses)]

    mappedByReadType = defaultdict(set)
    for readType in readTypes:
        for readFastqFileFullPath, samFile in samFiles[readType]:
            readFastqFile = os.path.basename(readFastqFileFullPath)
            mappedNames = {(x.qname, readFastqFile) for x in pysam.Samfile(samFile) if not x.is_unmapped}
            mappedByReadType[readType] = mappedByReadType[readType].union(mappedNames)

    unmappedByReadType = defaultdict(dict)
    for readType in readTypes:
        for readFastqFileFullPath, samFile in samFiles[readType]:
            readFastqFile = os.path.basename(readFastqFileFullPath)
            for name, seq, qual in fastqRead(readFastqFileFullPath):
                name = name.split(" ")[0]
                if (name, readFastqFile) not in mappedByReadType[readType]:
                    unmappedByReadType[readType][(name, readFastqFile)] = seq
        

    i = Stack(Target.makeTargetFn(find_analyses, args=(unmappedByReadType, outputDir))).startJobTree(options) 

    if i != 0:
        raise RuntimeError("Got {} failed jobs".format(i))

    for readType in readTypes:
        #build a counter of blast hits and set of read names that did not map
        blast_hits, no_hits = Counter(), set()
        for query, result in parse_blast(open(os.path.join(outputDir, readType + "_blast_out.txt"))):
            if result is None:
                no_hits.add(query)
            else:
                blast_hits[tuple(result)] += 1 #count number of times each hit was seen
        #write the unmapped hits to a fasta file
        outf = open(os.path.join(outputDir, readType + "_no_hits.fasta"), "w")
        for (name, readFastqFile), seq in unmappedByReadType[readType].iteritems():
            if name in no_hits:
                outf.write(">{}\n{}\n".format(name, seq))
        outf.close()
        #write the blast report
        blast_out = open(os.path.join(outputDir, readType + "_blast_report.txt"), "w")
        blast_out.write("gi|##|gb|##|\tSpecies\tseqID\tCount\n") #header to output
        for result, count in sorted(blast_hits.items(), key = lambda x: -int(x[-1])):
            blast_out.write("{}\t{}\n".format("\t".join(result), count))
        blast_out.close()
        #calculate percents and make a barplot
        blast_count =  sum(blast_hits.values())
        unmapped_count = len(unmappedByReadType[readType]) - sum(blast_hits.values())
        mapped_count = len(mappedByReadType[readType])
        
        #blast_percent = 1.0 * sum(blast_hits.values()) / (len(mappedByReadType[readType]) + len(unmappedByReadType[readType]))
        #unmapped_percent = (1.0 * len(unmappedByReadType[readType]) - sum(blast_hits.values())) / (len(mappedByReadType[readType]) + len(unmappedByReadType[readType]))
        #mapped_percent = 1.0 * len(mappedByReadType[readType]) / (len(mappedByReadType[readType]) + len(unmappedByReadType[readType]))
        outf = open(os.path.join(outputDir, readType + "percents.txt"),"w")
        outf.write("\n".join(map(str,[blast_count, unmapped_count, mapped_count])))
        outf.close()
        #system("Rscript blast_combined/barplot_blast.R {} {} {} {} {}".format(blast_percent, unmapped_percent, mapped_percent, readType, os.path.join(outputDir, readType + "_blast_barplot.pdf")))
        system("Rscript blast_combined/barplot_blast.R {} {} {} {} {}".format(blast_count, unmapped_count, mapped_count, readType, os.path.join(outputDir, readType + "_blast_barplot.pdf")))

if __name__ == "__main__":
    from scripts.blast_combined.blast_combined import *
    main()
