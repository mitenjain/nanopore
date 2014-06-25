from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from jobTree.src.bioio import fastaRead, system
import pysam, os

class KmerAnalysis(AbstractAnalysis):
    """Runs Karen's kmer analysis pipeline on a aligned sam and reference
    """
    def convert_sam_record(self, record, target):
        """Takes a single record from a sam and a str representing a sequence
        and converts this toper a aligned fasta (MAF) style format
        """
        read, ref = list(), list()
        pos, ref_pos = 0, record.pos
        for oper, num in record.cigar:
            if oper == 0 or oper == 7:
                #treating M and = equally
                read.append(record.seq[pos : pos + num])
                ref.append(record.seq[pos : pos + num])
                pos += num; ref_pos += num
            elif oper == 1:
                #insertion
                read.append(record.seq[pos : pos + num])
                ref.append("-" * num)
                pos += num
            elif oper == 2 or oper == 3:
                #treating D and N equally
                read.append("-" * num)
                ref.append(target[ref_pos : ref_pos + num])
                ref_pos += num
            elif oper == 8:
                #mismatch
                read.append(record.seq[pos : pos + num])
                ref.append(target[ref_pos : ref_pos + num])
                pos += num; ref_pos += num
            elif oper == 4:
                #soft clip means ignore that part of read
                pos += num
            elif oper == 5:
                #hard-clipped sequences can be ignored here
                continue
            elif oper == 6:
                #ignore padding as well
                continue
        return ("".join(read).upper(), "".join(ref).upper())
    def run(self, kmer_size=5):
        """Run karen's pipeline.
        """
        self.ref = dict(fastaRead(open(self.referenceFastaFile, "r")))
        #outf = open(os.path.join(self.getLocalTempDir(), "karen_tmp_fastaish"), "w")
        outf = open(os.path.join(self.getGlobalTempDir(),"karen_tmp_fastaish"), "w")
        sam = pysam.Samfile(self.samFile, "r" )
        for record in sam:
            rseq = self.ref[sam.getrname(record.tid)]
            seq, ref = self.convert_sam_record(record, rseq)
            outf.write("{}\t{}\n".format(seq, ref))
        outf.close()
        system("submodules/kmer/kmer.pl {} {} {}".format(self.readFastaFile, str(kmer_size), os.path.join(self.outputDir, "read_" + str(kmer_size) + "kmer")))
        system("submodules/kmer/kmer.pl {} {} {}".format(self.referenceFastaFile, str(kmer_size), os.path.join(self.outputDir, "ref_" + str(kmer_size) + "kmer")))
        system("submodules/kmer/cmpKmer.pl {} {} {}".format(os.path.join(self.outputDir, "ref_" + str(kmer_size) + "kmer"), os.path.join(self.outputDir, "read_" + str(kmer_size) + "kmer"), os.path.join(self.outputDir, str(kmer_size) + "kmer_Cmp.out")))
        system("submodules/kmer/kmer_del.pl {} {} {}".format(os.path.join(self.getGlobalTempDir(), "karen_tmp_fastaish"), os.path.join(self.outputDir, str(kmer_size) + "kmer_Cmp.out"), os.path.join(self.outputDir, "kmer_del.out")))
        system("submodules/kmer/kmer_ins.pl {} {} {}".format(os.path.join(self.getGlobalTempDir(), "karen_tmp_fastaish"), os.path.join(self.outputDir, str(kmer_size) + "kmer_Cmp.out"), os.path.join(self.outputDir, "kmer_ins.out")))
        
