from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
import pysam
from collections import Counter


class UnmappedBlastKmer(AbstractMetaAnalysis):
	"""Finds reads that do not map to any reference with any mapper
	blasts them, then does kmer analysis on those that do not have hits.
	Also reports stats on hits"""
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
				result = line.split("\t")[-3::]
				result[-1].rstrip()
				yield (query, result)
				result = None

	def run(self, kmer_size=5):
		for readType in self.info.readTypes:
			unmapped = {(x.name, x.seq, x.readFastqFile) for x in self.reads if x.readType == readType and x.is_mapped is False}
			if len(unmapped) > 1:
				outf = open(os.path.join(self.getLocalTempDir(), "unmapped.fasta"), "w")
				for name, seq, readFastqFile in unmapped:
					outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
				outf.close()

				system('blastn -outfmt "7 qseqid sseqid sscinames stitle" -db nt -max_target_seqs 1 -query {} -out {}'.format(os.path.join(self.getLocalTempDir(), "unmapped.fasta"), os.path.join(self.getLocalTempDir(), "blast_out.txt")))
				
				blast_hits, no_hits = Counter(), set()
				for query, result in self.parse_blast(open(os.path.join(self.getLocalTempDir(), "blast_out.txt"))):
					if result is None:
						no_hits.add(tuple(query.split(" "))) #need to save both read name and read fastq file 
					else:
						blast_hits[tuple(result)] += 1 #count number of times each hit was seen

				#generate a report on the blast hits that were returned
				blast_out = open(os.path.join(self.outputDir, readType + "_blast_report.txt"), "w")
				blast_out.write("gi|##|gb|##|\tSpecies\tseqID\tCount\n") #header to output
				for result, count in blast_hits.iteritems():
					blast_out.write("{}\t{}\n".format("\t".join(blast_hits[record])), count)
				blast_out.close()

				all_mapped = {(x.name, x.seq) for x in self.reads if x.is_mapped or (x.name, x.readFastqFile) not in no_hits}

				outf = open(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), "w")
				for name, seq in all_mapped:
					outf.write(">{}\n{}\n".format(name, seq))
				outf.close()
				outf = open(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), "w")
				for name, seq, readFastqFile in unmapped:
					outf.write(">{}\n{}\n".format(name, seq))
				outf.close()

				system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            	system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            	system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_" + str(kmer_size) + "kmer_Cmp.out")))
                system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(os.path.join(self.outputDir, readType + "_" + str(kmer_size)) + "kmer_Cmp.out"), os.path.join(self.outputDir, readType + "_top_kmers.tsv"), os.path.join(self.outputDir, readType + "_bot_kmers.tsv")))
          		
                
                all_blast = {(x.name, x.seq) for x in self.reads if (x.name, x.readFastqFile) not in no_hits and x.is_mapped is False}
                originally_mapped = {(x.name, x.seq) for x in self.reads if x.is_mapped is True}

                outf = open(os.path.join(self.getLocalTempDir(), "tmp"),"w")
                outf.write("{} {} {}\n".format(len(all_blast), len(originally_mapped), len(unmapped)))
                outf.close()
                system("Rscript nanopore.metaAnalyses/barplot_blast.R {} {}".format(os.path.join(self.getLocalTempDir(), "tmp"), os.path.join(self.outputDir), "blast_counts.pdf"), readType)




