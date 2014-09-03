#!/usr/bin/env python

import sys, os, pysam
from optparse import OptionParser
from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack 
from jobTree.src.bioio import fastqRead, fastaRead, setLoggingFromOptions, logger
from subprocess import Popen, PIPE
from itertools import izip

"""
Analysis script to determine if template/complement reads can be 'rescued' by providing the small region of the
reference where the 2D read they combined into aligned.

Should be ran via run_muscle.sh, which has three arguments:
--template_sam, --twoD_sam, --complement_sam

The sam files should be from the same aligner and read/reference set.
This program will dig through the nanopore pipeline's directory structure looking for reference and read files.
This assumes that you are on the version of the pipeline whereby reads are sorted by type into respective folders
(I.E. nanopore/readFastqFiles/template/a_template_file.fastq, etc)
"""


def find_analyses(target, recordsToAnalyze, templateFastqFiles, complementFastqFiles, references, outputDir):
    """takes a set of records to analyze and finds the corresponding sequences and creates alignment targets"""
    files = {"template":[], "complement":[]}

    logger.info("Finding template analyses")
    for fastqFile in templateFastqFiles:
        for name, seq, qual in fastqRead(fastqFile):
            if name in recordsToAnalyze:
                outfile = os.path.join(target.getGlobalTempDir(), "template_" + name)
                files["template"].append(outfile)
                ref_name, ref_start, ref_stop = recordsToAnalyze[name]
                ref_seq = references[ref_name][ref_start : ref_stop]
                analysis = [name, seq, ref_name, ref_seq, outfile]
                target.addChildTarget(Target.makeTargetFn(analyze, args=analysis))

    logger.info("Finding complement analyses")
    for fastqFile in complementFastqFiles:
        for name, seq, qual in fastqRead(fastqFile):
            if name in recordsToAnalyze:
                outfile = os.path.join(target.getGlobalTempDir(), "complement_" + name)
                files["complement"].append(outfile)
                ref_name, ref_start, ref_stop = recordsToAnalyze[name]
                ref_seq = references[ref_name][ref_start : ref_stop]
                analysis = [name, seq, ref_name, ref_seq, outfile]
                target.addChildTarget(Target.makeTargetFn(analyze, args=analysis))

    target.setFollowOnTargetFn(merge, args=(files, outputDir))

def analyze(target, name, seq, ref_name, ref_seq, outfile):
    """main analysis target; runs muscle on each pair of sequences"""
    outf = open(outfile, "w")
    p = Popen(["muscle", "-quiet"], stdout=PIPE, stdin=PIPE).communicate(">{}\n{}\n>{}\n{}\n".format(name, seq, ref_name, ref_seq))[0]
    outf.write(p); outf.close()

def merge(target, files, outputDir):
    """merges all muscle output into one fasta and runs metrics() on each"""
    for typeof in files:
        outmetrics = open(os.path.join(outputDir, typeof + "_metrics.tsv"), "w")
        outmetrics.write("Read\tReference\tMatches\tMismatches\tReadDeletionLength\tReadInsertionLength\tIdentity\tReferenceCoverage\n")
        for f in files[typeof]:
            handle = fastaRead(f)
            name, seq = handle.next()
            ref_name, ref_seq = handle.next()
            name = name.lstrip(">"); ref_name = ref_name.lstrip(">")
            outmetrics.write("\t".join([name, ref_name] + metrics(seq, ref_seq))); outmetrics.write("\n")
        outmetrics.close()
    
def metrics(seq, ref_seq):
    """takes in two aligned fasta sequences and calculates identity/coverage metrics"""
    matches = 0.0; mismatches = 0.0; readDeletionLength = 0.0; readInsertionLength = 0.0
    for s, r in izip(seq, ref_seq):
        if s == "-":
            readDeletionLength += 1
        elif r == "-":
            readInsertionLength += 1
        elif s == r == "-":
            continue #just in case?
        elif s == r:
            matches += 1
        elif s != r:
            mismatches += 1
    identity = matches / (matches + mismatches)
    referenceCoverage = (matches + mismatches) / (matches + mismatches + readDeletionLength)
    return map(str, [matches, mismatches, readDeletionLength, readInsertionLength, identity, referenceCoverage])

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    outputDir = "muscle_compare_2d/output/"

    if not os.path.exists(outputDir):
        logger.info("Output dir {} does not exist. Creating.")
        os.mkdir(outputDir)
    if len(os.listdir(outputDir)) > 0:
        logger.info("Output dir not empty.")

    if len(args) != 3:
        raise RuntimeError("Error: expected three arguments got %s arguments: %s" % (len(args), " ".join(args)))

    templateRecords = {x.qname for x in pysam.Samfile(args[0]) if not x.is_unmapped}
    complementRecords = {x.qname for x in pysam.Samfile(args[1]) if not x.is_unmapped}
    
    twodSamFile = pysam.Samfile(args[2])
    twodRecords = {x.qname : x for x in twodSamFile if not x.is_unmapped}

    recordsToAnalyze = dict()
    for name, record in twodRecords.iteritems():
        if name not in templateRecords and name not in complementRecords:
            ref_name = twodSamFile.getrname(record.tid)
            ref_start, ref_stop = int(record.aend - record.alen), int(record.aend)
            recordsToAnalyze[name] = [ref_name, ref_start, ref_stop]
    if os.path.exists("../readFastqFiles/template/") and os.path.exists("../readFastqFiles/complement"):
        templateFastqFiles = [os.path.join("../readFastqFiles/template/", x) for x in os.listdir("../readFastqFiles/template/") if x.endswith(".fastq") or x.endswith(".fq")]
        complementFastqFiles = [os.path.join("../readFastqFiles/complement/", x) for x in os.listdir("../readFastqFiles/complement/") if x.endswith(".fastq") or x.endswith(".fq")]
    else:
        raise RuntimeError("Error: readFastqFiles does not contain template and/or complement folders")

    referenceFastaFiles = [os.path.join("../referenceFastaFiles", x) for x in os.listdir("../referenceFastaFiles") if x.endswith(".fa") or x.endswith(".fasta")]
    
    if len(referenceFastaFiles) > 0:
        references = { y[0].split(" ")[0] : y[1] for x in referenceFastaFiles for y in fastaRead(x) }
    else:
        raise RuntimeError("Error: no reference fasta files")

    if len(recordsToAnalyze) == 0:
        raise RuntimeError("Error: none of the mappable twoD reads in this set did not map as template/complement.")

    logger.info("Starting to find analyses to run...")
    args = (recordsToAnalyze, templateFastqFiles, complementFastqFiles, references, outputDir)
    i = Stack(Target.makeTargetFn(find_analyses, args=args)).startJobTree(options) 

    if i != 0:
        raise RuntimeError("Got {} failed jobs".format(i))


if __name__ == "__main__":
    from scripts.muscle_compare_2d.muscle_compare_2d import *
    main()
