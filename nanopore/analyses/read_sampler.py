import pysam
import os, sys, glob, random

# using generator and yield function to read 4 lines at a time
def getStanza (infile):
    while True:
        fasta_id = infile.readline().strip()
        fasta_seq = infile.readline().strip()
        qual_id = infile.readline().strip()
        qual_scr = infile.readline().strip()
        if fasta_id != '':
            yield [fasta_id, fasta_seq, qual_id, qual_scr]
        else:
#             print >> sys.stderr, "Warning: End of Sequence"
            break

def SampleReads(workingDir):
    for readType in glob.glob(os.path.join(workingDir + "readFastqFiles", "*")):
        for readFastqFile in glob.glob(os.path.join(readType, "*")):
            if (not "percent" in readFastqFile and not "Consensus" in readFastqFile) and (".fq" in readFastqFile or ".fastq" in readFastqFile):
                for i in range(90, 0, -10):
                    newReadFastqFile = readFastqFile.split(".fastq")[0] + "_" + str(i) + "_percent.fastq"
                    if not os.path.exists(newReadFastqFile):
                        index = 0
                        readFastq = open(readFastqFile, "r")
                        output = len(readFastq.readlines())
                        readFastq.close()
                        # number of sequences
                        total_seqs = output / 4
                        num_seq = total_seqs * i / 100
                        indexes = sorted(random.sample(range(total_seqs), num_seq))
                        readFastq = open(readFastqFile, "r")
                        newreadFastq = open(newReadFastqFile, "w")
                        for stanza in getStanza(readFastq):
                            if index in indexes:
                                newreadFastq.write("\n".join(stanza))
                                newreadFastq.write("\n")
                            index += 1
                        readFastq.close()
                        newreadFastq.close()
