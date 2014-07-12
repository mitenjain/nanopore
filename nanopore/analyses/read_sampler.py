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
    for readFastqFile in glob.glob(os.path.join(workingDir + "readFastqFiles", "*")):
        if not "%" in readFastqFile and (".fq" in readFastqFile or ".fastq" in readFastqFile):
            for i in range(75, 1, -25):
                newReadFastqFile = readFastqFile.split(".fastq")[0] + "_" + str(i) + "%.fastq"
                if not os.path.exists(newReadFastqFile):
                    index = 0
                    readFastq = open(readFastqFile, "rb")
                    output = len(readFastq.readlines())
                    readFastq.close()
                    # number of sequences
                    total_seqs = output / 4
                    num_seq = total_seqs * i / 100
                    indexes = sorted(random.sample(range(total_seqs), num_seq))
                    readFastq = open(readFastqFile, "rb")
                    newreadFastq = open(newReadFastqFile, "wb")
                    for stanza in getStanza(readFastq):
                        if index in indexes:
                            newreadFastq.write("\n".join(stanza))
                            newreadFastq.write("\n")
                        index += 1
                    readFastq.close()
                    newreadFastq.close()

