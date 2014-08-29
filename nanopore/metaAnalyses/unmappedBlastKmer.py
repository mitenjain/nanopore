from nanopore.metaAnalyses.abstractUnmappedAnalysis import AbstractUnmappedMetaAnalysis
import os, sys
import xml.etree.cElementTree as ET
from jobTree.src.bioio import system, fastqRead
from nanopore.analyses.utils import samIterator
import pysam
from collections import Counter


class UnmappedBlastKmer(AbstractUnmappedMetaAnalysis):
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
				result = line.strip().split("\t")[-3::]
				yield (query, result)
				result = None

	def run(self, kmer_size=5):
		for readType in self.readTypes:
			unmapped = {(x.name, x.seq, x.readFastqFile) for x in self.reads if x.readType == readType and x.is_mapped is False}
			if len(unmapped) > 1:
				outf = open(os.path.join(self.outputDir, readType + "_unmapped.fasta"), "w")
				for name, seq, readFastqFile in unmapped:
					outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
				outf.close()

				system('blastn -outfmt "7 qseqid sseqid sscinames stitle" -db nt -max_target_seqs 1 -query {} -out {}'.format(os.path.join(self.outputDir, "unmapped.fasta"), os.path.join(self.outputDir, readType + "_blast_out.txt")))
				
				blast_hits, no_hits = Counter(), set()
				for query, result in self.parse_blast(open(os.path.join(self.outputDir, readType + "_blast_out.txt"))):
					if result is None:
						no_hits.add(tuple(query.split(" "))) #need to save both read name and read fastq file 
					else:
						blast_hits[tuple(result)] += 1 #count number of times each hit was seen

				#generate a report on the blast hits that were returned
				blast_out = open(os.path.join(self.outputDir, readType + "_blast_report.txt"), "w")
				blast_out.write("gi|##|gb|##|\tSpecies\tseqID\tCount\n") #header to output
				for result, count in sorted(blast_hits.items(), key = lambda x: -int(x[-1])):
					blast_out.write("{}\t{}\n".format("\t".join(result), count))
				blast_out.close()

				#generate a set of reads that mapped to anything
				all_mapped = {(x.name, x.seq, x.readFastqFile) for x in self.reads if x.is_mapped or (x.name, x.readFastqFile) not in no_hits and x.readType == readType}

				#generate a set of reads that mapped to nothing
				no_mappings = {(x.name, x.seq, x.readFastqFile) for x in self.reads if (x.name, x.readFastqFile) in no_hits and x.readType == readType}


				outf = open(os.path.join(self.outputDir, readType + "_mapped_reads.fasta"), "w")
				for name, seq, readFastqFile in all_mapped:
					outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
				outf.close()

				outf = open(os.path.join(self.outputDir, readType +  "_no_hits.fasta"), "w")
				for name, seq, readFastqFile in no_mappings:
					outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
				outf.close()

				system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_no_hits.fasta"), os.path.join(self.outputDir, "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), str(kmer_size)))
				system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_mapped_reads.fasta"), os.path.join(self.outputDir, "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), str(kmer_size)))
				system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.outputDir, readType + "_mapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_" + str(kmer_size) + "kmer_Cmp.out")))
				system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(os.path.join(self.outputDir, readType + "_" + str(kmer_size)) + "kmer_Cmp.out"), os.path.join(self.outputDir, readType + "_top_kmers.tsv"), os.path.join(self.outputDir, readType + "_bot_kmers.tsv")))

				outf = open(os.path.join(self.outputDir, readType + "_tmp"),"w")
				blast_percent = (1.0 * sum(blast_hits.values())) / len(self.reads)
				unmapped_percent = (1.0 * len(no_mappings)) / len(self.reads)
				mapped_percent = 1 - unmapped_percent
				outf.write("{} {} {}\n".format(blast_percent, mapped_percent, unmapped_percent))
				outf.close()
				system("Rscript nanopore/metaAnalyses/barplot_blast.R {} {} {}".format(os.path.join(self.outputDir, readType + "tmp"), os.path.join(self.outputDir, readType + "blast_counts.pdf"), readType))


