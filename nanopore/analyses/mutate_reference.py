from random import random, choice
import os, sys, glob

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
				i = 100 - (mutation_rate * 100)
				newreferenceFastaFile = referenceFastaFile.split(".fa")[0] + "_" +  str(i) + "_percent.fasta"
				mutationIndexFile = referenceFastaFile.split(".fa")[0] + "_" + str(i) + "_percent_Mutation_Index.txt"
				if not os.path.exists(newreferenceFastaFile):
					referenceFasta = open(referenceFastaFile, "r")
					newreferenceFasta = open(newreferenceFastaFile, "w")

					mutationIndex = open(mutationIndexFile, "w")
					referenceFastaFile_ID = referenceFastaFile.split(".fa")[0].split("/")[-1]
					mutationIndex.write("##ID\t" + referenceFastaFile_ID)
					mutationIndex.write("\n")
					mutationIndex.write("#POS\tREF\tMUT")
					mutationIndex.write("\n")

					for seq in Fastaseq.readline(referenceFasta):
						temp_seq = list(seq.seq)
						for index, char in enumerate(temp_seq):
							error = random()
							if error < mutation_rate:
								# choose a random nucleotide that's different
								temp_seq[index] = choice([x for x in "ACGT" if x != char])
								mutationIndex.write(str(index + 1) + "\t" + char + "\t" + temp_seq[index])
								mutationIndex.write("\n")

						newreferenceFasta.write(seq.id)
						newreferenceFasta.write("\n")
						newreferenceFasta.write("".join(temp_seq))
						newreferenceFasta.write("\n")
					referenceFasta.close()
					newreferenceFasta.close()
					mutationIndex.close()
