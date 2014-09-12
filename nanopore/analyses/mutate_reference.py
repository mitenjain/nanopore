import os, sys, random
from sonLib.bioio import fastaRead, fastaWrite

'''
Inserts SNPs and 20% of SNPs as InDels
Outputs Two files for every file - reference.fasta
1. mutated reference - reference_X_percent_SNPs_Y_percent_InDels.fasta
2. mutation index for true and mutated reference - reference_X_percent_SNPs_Y_percent_InDels_Index.txt
'''

def mutateSequence(sequence, snpRate): #Does not preserve softmasking
    return "".join(map(lambda x : x.upper() if random.random() >= snpRate else random.choice(list(set(("A", 'C', 'G', 'T')) - set(x.upper()))), sequence))

def mutateReferenceSequences(referenceFastaFiles):
    updatedReferenceFastaFiles = referenceFastaFiles[:]
    for referenceFastaFile in referenceFastaFiles:
        if not "percent" in referenceFastaFile:
            mutation_rates = [0.01, 0.05, 0.10, 0.20]
            for mutation_rate in mutation_rates:
                indel_rate = 0.0 * mutation_rate # indel rate = 20% of Substitution rate
                i = mutation_rate * 100
                j = indel_rate * 100
                newReferenceFastaFile = referenceFastaFile.split(".fa")[0] + "_" + str(i) + "_percent_SNPs_" + str(j) + "_percent_InDels.fasta"
                mutationIndexFile = referenceFastaFile.split(".fa")[0] + "_" + str(i) + "_percent_SNPs_" + str(j) + "_percent_InDels.fasta_Index.txt"
                updatedReferenceFastaFiles.append(newReferenceFastaFile)
                if not os.path.exists(newReferenceFastaFile):
                    fH = open(newReferenceFastaFile, 'w')
                    fH2 = open(mutationIndexFile, 'w')
                    for header, seq in fastaRead(referenceFastaFile):
                        header = header.split()[0]
                        mutatedSeq = mutateSequence(seq, mutation_rate)
                        fastaWrite(fH, header, mutatedSeq)
                        fastaWrite(fH2, header, seq)
                        fastaWrite(fH2, header + "_mutated", mutatedSeq)
                    fH.close()
                    fH2.close()
    return updatedReferenceFastaFiles
