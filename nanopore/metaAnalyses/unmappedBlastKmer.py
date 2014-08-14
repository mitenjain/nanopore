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
				result = line.split("\t")[-3::]
				result[-1].rstrip()
				yield (query, result)
				result = None


	def run(self, kmer_size=5):	
		#build set of (read, fastq) of all mapped reads
		mappedReads, unmappedReads = dict(), dict()
		for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
			if readType not in mappedReads:
				mappedReads[readType] = set()
			for record in samIterator(pysam.Samfile(os.path.join(resultsDir, "mapping.sam"))):
				if not record.is_unmapped: #some aligners save unmapped reads in their samfiles (bwa)
					mappedReads[readType].add((record.qname, readFastqFile))
					#can't save seqs here because sam file does not have full sequence in all aligners (no hard clipping)
		
		#loop again and save only reads that did not map with any mapper
		for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
			if readType not in unmappedReads:
				unmappedReads[readType] = set()
			for name, seq, qual in fastqRead(readFastqFile):
				if (name, readFastqFile) not in mappedReads:
					unmappedReads[readType].add((name, readFastqFile, seq))

		#blast these unmapped reads
		for readType in unmappedReads:
			if len(unmappedReads[readType]) >= 1:
				outf = open(os.path.join(self.getLocalTempDir(), "unmapped.fasta"), "w")
				for name, readFastqFile, seq in unmappedReads[readType]:
					outf.write(">{} {}\n{}\n".format(name, readFastqFile, seq))
				outf.close()
				system('blastn -outfmt "7 qseqid sseqid sscinames stitle" -db nt -max_target_seqs 1 -query {} -out {}'.format(os.path.join(self.getLocalTempDir(), "unmapped.fasta"), os.path.join(self.getLocalTempDir(), "blast_out.txt")))
				blast_hits, no_hits = Counter(), set()
				for query, result in parse_blast(open(os.path.join(self.getLocalTempDir(), "blast_out.txt"))):
					if result is None:
						no_hits.add(tuple(query.split(" "))) #need to save both read name and read fastq file 
					else:
						blast_hits[result] += 1 #count number of times each hit was seen

				#generate a report on the blast hits that were returned
				blast_out = open(os.path.join(self.outputDir, readType + "_blast_report.txt"), "w")
				blast_out.write("gi|##|gb|##|\tSpecies\tseqID\tCount\n") #header to output
				for record in blast_hits:
					blast_out.write("{}\t{}\n".format(record, "\t".join(blast_hits[record])))
				blast_out.close()

				#build a list of no-hit reads to be used in the kmer pipeline
				no_hit_reads = [(name, seq) for readFastqFile, thisReadType, referenceFastaFile, mapper, analyses, resultsDir \
						in self.experiments for name, seq, qual in fastqRead(readFastqFile) if (name, readFastqFile) in no_hits \
						and thisReadType == readType]

				#write to file for Karen's perl scripts
				outf = open(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), "w")
				for name, seq in no_hit_reads:
					outf.write(">{}\n{}\n".format(name, seq))
				outf.close()

				#build an equivalent list of mapped reads (both blast and to reference)
				hit_reads = [(name, seq) for readFastqFile, thisReadType, referenceFastaFile, mapper, analyses, resultsDir \
						in self.experiments for name, seq, qual in fastqRead(readFastqFile) if (name, readFastqFile) not in no_hits \
						or (name, readFastqFile) in mapped_reads[readType] and thisReadType == readType]

				#also write these to a file
				outf = open(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), "w")
				for name, seq in hit_reads:
					outf.write(">{}\n{}\n".format(name, seq))
				outf.close()

				#run kmer scripts
				system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "unmapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            	system("nanopore/analyses/kmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "mapped_reads.fasta"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), str(kmer_size)))
            	system("nanopore/analyses/cmpKmer.pl {} {} {}".format(os.path.join(self.getLocalTempDir(), "readType_" + readType + "_mapped_" + str(kmer_size) + "mer"), os.path.join(self.getLocalTempDir(), "readType_" + readType + "_unmapped_" + str(kmer_size) + "mer"), os.path.join(self.outputDir, readType + "_" + str(kmer_size) + "kmer_Cmp.out")))
                system("Rscript nanopore/analyses/kmer_most_under_over.R {} {} {}".format(os.path.join(os.path.join(self.outputDir, readType + "_" + str(kmer_size)) + "kmer_Cmp.out"), os.path.join(self.outputDir, readType + "_top_kmers.tsv"), os.path.join(self.outputDir, readType + "_bot_kmers.tsv")))
          



