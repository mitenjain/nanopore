from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from sonLib.bioio import system
from nanopore.analyses.utils import samToBamFile
import os
import pysam

def formatConsensusFastq(inputConsensusFastq, outputConsensusFastq):
    infile = open(inputConsensusFastq, "r")
    outfile = open(outputConsensusFastq, "w")
    fasta_id = None
    fasta_seq = {}
    qual_id = {}
    qual_scr = {}
    flag = None
    for line in infile:
        if not line == "":
            if line.startswith("@") and not flag == "qual":
                fasta_id = line.strip()
                flag = "seq"
                fasta_seq[fasta_id] = []
                continue
            elif line.startswith("+"):
                qual_id[fasta_id] = line.strip()
                flag = "qual"
                qual_scr[fasta_id] = []
                continue
            if flag == "seq":
                fasta_seq[fasta_id].append(line.upper().strip())
            elif flag == "qual":
                qual_scr[fasta_id].append(line.strip())
        else:
            #print >> sys.stderr, "Warning: End of Sequence"
            break

    for fasta_id in fasta_seq.keys():
        if len(fasta_seq[fasta_id]) > 0:
            outfile.write(fasta_id)
            outfile.write("\n")
            outfile.write("".join(fasta_seq[fasta_id]))
            outfile.write("\n")
            outfile.write(qual_id[fasta_id])
            outfile.write("\n")
            outfile.write("".join(qual_scr[fasta_id]))
            outfile.write("\n")
    
    infile.close()
    outfile.close()

class Consensus(AbstractAnalysis):
    def run(self):
        AbstractAnalysis.run(self) #Call base method to do some logging
        localBamFile = os.path.join(self.getLocalTempDir(), "mapping.bam")
        localSortedBamFile = os.path.join(self.getLocalTempDir(), "mapping.sorted")

        samToBamFile(self.samFile, localBamFile)
        pysam.sort(localBamFile, localSortedBamFile)
        pysam.index(localSortedBamFile + ".bam")
        pysam.faidx(self.referenceFastaFile)
        
        file_header = self.readFastqFile.split(".fastq")[0].split("/")[-1] +  "_" + self.referenceFastaFile.split(".fa")[0].split("/")[-1]
        consensus_vcf = os.path.join(self.outputDir, file_header + "_Consensus.vcf")
        consensus_fastq = os.path.join(self.outputDir, file_header + "_Consensus.fastq")

        system("samtools mpileup -Q 0 -uf %s %s | bcftools view -cg - > %s" \
                % (self.referenceFastaFile, localSortedBamFile + ".bam", consensus_vcf))
        system("vcfutils.pl vcf2fq %s > %s" % (consensus_vcf, consensus_fastq))
        system("rm -rf %s" % (self.referenceFastaFile + ".fai"))
        
        formatted_consensus_fastq = os.path.join(self.getLocalTempDir(), "Consensus.fastq")
        
        formatConsensusFastq(consensus_fastq, formatted_consensus_fastq)
        system("mv %s %s" % (formatted_consensus_fastq, consensus_fastq))
        cns_to_reads_folder = "./readFastqFiles/" + file_header + "_Consensus.fastq"
        if os.stat(consensus_fastq).st_size > 4:
            system("cp %s %s" % (consensus_fastq, cns_to_reads_folder))
        
        self.finish()
