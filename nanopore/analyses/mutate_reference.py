import random
import sys, os, glob, itertools

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

def MutateReference(workingDir):
	for referenceFastaFile in glob.glob(os.path.join(workingDir + "referenceFastaFiles", "*")):
		if not "percent" in referenceFastaFile and (".fa" in referenceFastaFile or ".fasta" in referenceFastaFile):
			indel_rates = [0.10]
			mutation_rates = [0.01, 0.05, 0.10, 0.20]
			for indel_rate in indel_rates:
				j = indel_rate * 100
				for mutation_rate in mutation_rates:
					i = mutation_rate * 100
					newreferenceFastaFile = referenceFastaFile.split(".fa")[0] + "_" + str(i) + "_percent_SNPs_" + str(j) + "_percent_InDels.fasta"
					mutationIndexFile = referenceFastaFile.split(".fa")[0] + "_" + str(i) + "_percent_SNPs_" + str(j) + "_percent_InDels_Index.txt"
					if not os.path.exists(newreferenceFastaFile):
						referenceFasta = open(referenceFastaFile, "r")
						newreferenceFasta = open(newreferenceFastaFile, "w")
						mutationIndex = open(mutationIndexFile, "w")
						
						referenceFastaFile_ID = referenceFastaFile.split(".fa")[0].split("/")[-1]
						for seq in Fastaseq.readline(referenceFasta):
							temp_seq = []
							main_seq = list(seq.seq)
							ref_seq = []
							flag = False
							index = None
							new_seq_id = seq.id + " " + str(i) + " percent SNPs " + str(j) + " percent InDels"
							for position in xrange(len(main_seq)):
								char = main_seq[position]
								mut = ""
								ins = ""
								snp = ""
								indel = ""
								inserts = ""
								deletes = ""
								
								ref_seq.append(char)
								if position <= index and flag == True:
									continue
								else:
									index = None
									flag = False

								snp = random.random()
								if snp < mutation_rate and flag == False: # Add SNP
									# choose a random nucleotide that's different
									mut = random.choice([x for x in "ACGT" if x != char])
									temp_seq.append(mut)
								else:
									temp_seq.append(char)

								indel = random.random()
								if indel < indel_rate and flag == False: # Add InDel
									indel_type = random.choice([m for m in "ID"])
									if indel_type == "I":
										inserts = int(random.expovariate(1)) + 1
										ins = "".join([random.choice([x for x in "ACGT"]) for j in range(inserts)])
										temp_seq.append(ins)
										ref_seq.append("-" * inserts)
									elif indel_type == "D":
										deletes = int(random.expovariate(1)) + 1
										temp_seq.append("-" * deletes)
										index = position + deletes
										flag = True

						mutationIndex.write(seq.id)
						mutationIndex.write("\n")
						mutationIndex.write("".join(list(itertools.chain(*ref_seq))))
						mutationIndex.write("\n")
						mutationIndex.write(new_seq_id)
						mutationIndex.write("\n")
						mutationIndex.write("".join(list(itertools.chain(*temp_seq))))
						mutationIndex.write("\n")

						new_ref_seq = "".join(list(itertools.chain(*temp_seq))).replace("-", "")
						newreferenceFasta.write(new_seq_id)
						newreferenceFasta.write("\n")
						newreferenceFasta.write(new_ref_seq)
						newreferenceFasta.write("\n")
						referenceFasta.close()
						newreferenceFasta.close()
						mutationIndex.close()
