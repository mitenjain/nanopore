import os, sys, glob, random, itertools

'''
Inserts SNPs and 20% of SNPs as InDels
Outputs Two files for every file - reference.fasta
1. mutated reference - reference_X_percent_SNPs_Y_percent_InDels.fasta
2. mutation index for true and mutated reference - reference_X_percent_SNPs_Y_percent_InDels_Index.txt
'''

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
			mutation_rates = [0.01, 0.05, 0.10, 0.20]
			for mutation_rate in mutation_rates:
				indel_rate = 0.2 * mutation_rate # indel rate = 20% of Substitution rate
				i = mutation_rate * 100
				j = indel_rate * 100
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
						positions = []
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
								positions.append(position + 1)
								continue
							else:
								index = None
								flag = False

							snp = random.random()
							if snp < mutation_rate and flag == False: # Add SNP
								# choose a random nucleotide that's different
								mut = random.choice([x for x in "ACGT" if x != char])
								temp_seq.append(mut)
								positions.append(position + 1)
							else:
								temp_seq.append(char)
								positions.append(position + 1)

							indel = random.random()
							if indel < indel_rate and flag == False: # Add InDel
								indel_type = random.choice([m for m in "ID"])
								if indel_type == "I":
									inserts = int(random.expovariate(1)) + 1
									ins = "".join([random.choice([x for x in "ACGT"]) for j in range(inserts)])
									temp_seq.append(ins)
									positions.append(",".join("-" * inserts))
									ref_seq.append("-" * inserts)
								elif indel_type == "D":
									deletes = int(random.expovariate(1)) + 1
									temp_seq.append("-" * deletes)
									index = position + deletes
									flag = True

					mutationIndex.write(",".join(map(str, positions)))
					mutationIndex.write("\n")
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
