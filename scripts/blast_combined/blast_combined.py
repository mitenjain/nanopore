#!/usr/bin/env python

import sys, os, pysam
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 
from jobTree.src.bioio import fastqRead, fastaRead, setLoggingFromOptions, logger, popenCatch, system
from itertools import izip, product
from collections import Counter

"""Used to BLAST results from the combined analysis result in order to find the set of reads that truly unmappable.
Reports the BLAST results as well as the raw hits, and also reports a fasta of the reads that did not map anywhere,
and generates a summary barplot."""

readTypes = ["2D", "template", "complement"]
combinedAnalyses = ["CombinedMapper", "CombinedMapperChain", "CombinedMapperRealign", "CombinedMapperRealignEm", "CombinedMapperRealignTrainedModel"]

def parse_blast(self, blast_handle):
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
            result = None

def find_analyses(target, unmappedByReadType, outputDir):
    outfiles = dict()
    for readType in unmappedByReadType:
        outfiles[readType] = list()
        records = list()
        for (sequence, (name, readFastqFile)), i in izip(unmappedByReadType[readType].iteritems(), xrange(len(unmappedByReadType[readType]))):
                records.append(">{} {}\n{}\n".format(name, os.path.basename(readFastqFile), sequence))
                if i % 30 == 0 or i == len(unmappedByReadType[readType]) - 1:
                    tmpalign = os.path.join(target.getGlobalTempDir(), str(i) + ".fasta")
                    outfiles[readType].append(tmpalign)
                    target.addChildTarget(Target.makeTargetFn(run_blast, args=(records, tmpalign)))
                    records = list()
    target.setFollowOnTargetFn(merge, args=(records, tmpalign))

def run_blast(target, records, tmpalign):
    query = "".join(records)
    p = popenCatch('blastn -outfmt "7 qseqid sseqid sscinames stitle" -db nt -max_target_seqs 1', stdinString=query)
    outf = open(tmpalign, "w"); outf.write(p); outf.close()


def merge(target, outfiles, outputDir):
    for readType in outfiles:
        outf = os.path.join(outputDir, readType + "_raw_blast.txt")
        with open(outf, "w") as outfile:
            for f in outfiles[readType]:
                with open(f) as i:
                    outf.write(i.read())
        outf.close()

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
        readFastqFiles[readType] = [os.path.join("../readFastqFiles", readType, x) for x in os.listdir(os.path.join("../readFastqFiles", readType)) if x.endswith(".fq") or x.endswith(".fastq")]
    
    #find all reference fasta files
    referenceFastaFiles = [x for x in os.listdir("../referenceFastaFiles") if x.endswith(".fasta") or x.endswith(".fa")]

    #find all sam files that were analyzed using combinedAnalyses
    samFiles = dict()
    for readType in readTypes:
        samFiles[readType] = [(readFastqFile, referenceFastaFile, os.path.join("../output", "analysis_" + readType, "experiment_" + os.path.basename(readFastqFile) + "_" + referenceFastaFile + "_" + analysis, "mapping.sam")) for readFastqFile, referenceFastaFile, analysis in product(readFastqFiles[readType], referenceFastaFiles, combinedAnalyses)]

    mappedByReadType = dict(); unmappedByReadType = dict()
    for readType in readTypes:
        unmappedReads = dict(); mappedReads = dict()
        for readFastqFile, referenceFastaFile, samFile in samFiles[readType]:
            mappedNames = {x.qname for x in pysam.Samfile(samFile) if not x.is_unmapped}
            for name, seq, qual in fastqRead(readFastqFile):
                if name not in mappedNames:
                    unmappedReads[seq] = (name, readFastqFile)
                else:
                    mappedReads[seq] = (name, readFastqFile)
        unmappedByReadType[readType] = unmappedReads
        mappedByReadType[readType] = mappedReads

    i = Stack(Target.makeTargetFn(find_analyses, args=(unmappedByReadType, outputDir))).startJobTree(options) 

    if i != 0:
        raise RuntimeError("Got {} failed jobs".format(i))

    blast_hits, no_hits = Counter(), set()
    for query, result in parse_blast(open(os.path.join(outputDir, readType + "_blast_out.txt"))):
        if result is None:
            no_hits.add(tuple(query.split(" "))) #need to save both read name and read fastq file 
        else:
            blast_hits[tuple(result)] += 1 #count number of times each hit was seen

    for readType in readTypes:
        outf = open(os.path.join(outputDir, readType + "_no_hits.fasta"), "w")
        for seq, (name, readFastqFile) in unmappedReads.iteritems():
            if (name, readFastqFile) in no_hits:
                outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
        outf.close()

    for readType in readTypes:
        blast_out = open(os.path.join(outputDir, readType + "_blast_report.txt"), "w")
        blast_out.write("gi|##|gb|##|\tSpecies\tseqID\tCount\n") #header to output
        for result, count in sorted(blast_hits.items(), key = lambda x: -int(x[-1])):
            blast_out.write("{}\t{}\n".format("\t".join(result), count))
        blast_out.close()

    for readType in readTypes:
        blast_percent = (1.0 * sum(blast_hits.values()) / len(mappedByReadType[readType]) + len(unmappedByReadType[readType]))
        unmapped_percent = (1.0 * len(unmappedByReadType[readType]) - sum(blast_hits.values())) / (len(mappedByReadType[readType]) + len(unmappedByReadType[readType]))
        system("Rscript blast_combined/barplot_blast.R {} {} {}".format(blast_percent, unmapped_percent, os.path.join(outputDir, "blast_barplot.pdf")))

if __name__ == "__main__":
    from scripts.blast_combined.blast_combined import *
    main()
