from nanopore.metaAnalyses.abstractMetaAnalysis import AbstractMetaAnalysis
from sonLib.bioio import system
import os, sys, glob, subprocess, numpy
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
 
					# macro numbers for coverage
 					all_cov = map(int, subprocess.check_output('''awk '{split($0, a, "\t"); print a[3]}' %s''' % (depthFile), shell=True).strip().split("\n"))
					mean_cov = numpy.mean(all_cov)
					sd_cov = numpy.std(all_cov)
					#Threshold is 2X SD, for 95%
					sd_threshold = 2 * sd_cov

 					cov_stats = open(covStats, "w")
 					cov_stats.write("Position\tCoverage (mu=" + str(mean_cov) + "X, sd=" + str(sd_cov) + "X)\tKmer\n")
					
 					prev_pos_cov = 0 # Counter for previous position coverage
 					for line in depth_file:
 						line = line.strip().split("\t")
 						if int(line[2]) - prev_pos_cov >= sd_threshold:
 							if int(line[1]) >= 5:
 								cov_stats.write(str(line[1]) + "\t" + str(line[2]) + "\t" + str(ref_seq[(int(line[1])) - 5: int(line[1])]) + "\n")
 							else:
 								cov_stats.write(str(line[1]) + "\t" + str(line[2]) + "\t" + str(ref_seq[0: int(line[1])]) + "\n")
 						prev_pos_cov = int(line[2]) # Reassign
 	
 					depth_file.close()
 					cov_stats.close()
 					system("Rscript nanopore/metaAnalyses/coverageDepth_plot.R {} {}".format(depthFile, os.path.join(self.outputDir, experiment + ".pdf")))
				except:
					continue
