from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
from sonLib.bioio import system
import os, sys, glob
from nanopore.analyses.utils import samToBamFile
import pysam

class Fastaseq():
	"""
	fasta reader
	"""
  def __init__(self):
		self.id = None
		self.seq = ''
		self.length = ''
	  
	@staticmethod 
	def readline(linein):
		seqobj = Fastaseq()
		for line in linein:
			if len(line)==0: 
				print >> sys.stderr, 'empty line'
				continue
			if line.startswith('>'):
				if seqobj.id is None:
					seqobj.id = line.rstrip()
					continue
				else:
					yield seqobj
					seqobj = Fastaseq()
					seqobj.id = line.rstrip()
			else:
				seqobj.seq += line.rstrip('\n\r').upper()
		yield seqobj

class CoverageDepth(AbstractMetaAnalysis):
	"""
	Uses samtoold depth to obtain and plot coverage depth per base across reference
	"""
	def __init__(self, outputDir, experiments):
		AbstractMetaAnalysis.__init__(self, outputDir, experiments)
		parentFolder = self.outputDir + "/"

		experiments = []
		for readFastqFile, readType, referenceFastaFile, mapper, analyses, resultsDir in self.experiments:
			experiment = resultsDir.split("/")[-1]
			experiments.append(experiment)

			# Check and create bam, sorted bam, and indexed bam files
			bamFile = os.path.join(resultsDir, "mapping.bam")
			sortedbamFile = os.path.join(resultsDir, "mapping.sorted")
			depthFile = os.path.join(self.outputDir, experiment + "_Depth.txt")
			if not os.path.isfile(sortedbamFile):
				try:
					samToBamFile(os.path.join(resultsDir, "mapping.sam"), bamFile)
					pysam.sort(bamFile, sortedbamFile)
					pysam.index(sortedbamFile + ".bam")
				except:
					continue
			else:
				continue
			
			# if sorted bam file present then compute coverage and plot
			if os.path.isfile(sortedbamFile + ".bam"):
				try:
					system("samtools depth %s > %s" % (sortedbamFile + ".bam", depthFile))
 					covStats = parentFolder + experiment + "_Stats.out"
 					fastaFile = open(referenceFastaFile, "r")
 					depth_file = open(depthFile, "r")
 					ref_seq = None
 					for seq in Fastaseq.readline(fastaFile):
 						ref_seq = seq.seq
 					fastaFile.close()
 
 					cov_stats = open(covStats, "w")
 					cov_stats.write("Position\tCoverage\tKmer\n")
 					for line in depth_file:
 						line = line.strip().split("\t")
 						if int(line[2]) <= 30:
 							if int(line[1]) >= 5:
 								cov_stats.write(str(line[1]) + "\t" + str(line[2]) + "\t" + str(ref_seq[(int(line[1])) - 5: int(line[1])]) + "\n")
 							else:
 								cov_stats.write(str(line[1]) + "\t" + str(line[2]) + "\t" + str(ref_seq[0: int(line[1])]) + "\n")
 					depth_file.close()
 					cov_stats.close()
 					system("Rscript nanopore/metaAnalyses/coverageDepth_plot.R {} {}".format(depthFile, os.path.join(self.outputDir, experiment + ".pdf")))
				except:
					continue
