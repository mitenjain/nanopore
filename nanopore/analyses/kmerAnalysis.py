from abstractAnalysis import AbstractAnalysis
import re, subprocess, sys

class KmerAnalysis(AbstractAnalysis):
    """Runs Karen's kmer analysis pipeline on a aligned sam and reference
    """
    def __init__(self):
        self.cigar_regex = re.compile("\d+[MIDNSHPX=]")
    def parse_sam(sam_handle):
        """Generator to yield select entries from a sam file
        """
        for line in sam_handle:
            if not line.startswith("@"):
                line = line.split()
                record = dict()
                record['qname'] = line[0]
                record['rname'] = line[2]
                record['pos'] = line[3] 
                record['cigar'] = re.findall(self.cigar_regex, line[5])
                record['seq'] = line[9]
                yield record
    def read_fasta(fasta_handle):
        """Generator to yield tuple representing one entry in a fasta file
        tuple format: (<id>,<comments>,<sequence>)
        Source: BioPython source code"""
        name, comments, seq = None, None, []
        for line in file_handle:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    yield (name, comments, ''.join(seq))
                line = line[1:].split()
                name, comments, seq = line[0], line[1:], []
            else:
                line = ''.join([x for x in line if not x.isdigit() and not x.isspace()])
                seq.append(line)
        if name:
            yield (name, comments, ''.join(seq))
    def load_reference():
        """Loads reference fasta into memory"""
        self.ref = dict()
        for name, com, seq in read_fasta(open(self.referenceFastaFile)):
            self.ref[name] = seq
    def convert_sam_record(record, target):
        """Takes a single record from parse_sam() and a str representing a sequence
        and converts this to a aligned fasta (MAF) style format
        """
        read, ref = list(), list()
        pos, ref_pos = 0, record['pos']
        for num, char in record['cigar']:
            if char == "M":
                read.append(record['seq'][pos : pos + num])
                ref.append(record['seq'][pos : pos + num])
                pos += num; ref_pos += num
            elif char == "I":
                read.append(record['seq'][pos : pos + num])
                ref.append("-" * num)
                pos += num
            elif char == "D" or char == "N":
                #for this format, deletion is same as skipped
                read.append("-" * num)
                ref.append(target[ref_pos : ref_pos + num])
                ref_pos += num
            elif char == "S":
                #soft clip means ignore that part of read
                pos += num
            elif char == "H":
                #hard-clipped sequences can be ignored here
                continue
            elif char == "P":
                #ignore padding as well
                continue
        return ("".join(read), "".join(ref))
    def sam_to_fa():
        """Converts sam to pairs of aligned fasta records
        """
        for record in parse_sam(open(self.samFile)):
            read, ref = convert_sam_record(record, self.ref[record['rname']])
            yield ">{}\n{}\n>{}\n{}\n".format(record['qname'], read, \
                    record['rname'], ref)
    def sam_to_karen():
        """Same as sam_to_fa() but converts to tab separated sequences - 
        this is how karen's pipeline works
        """
        seqs = list()
        for record in parse_sam(open(self.samFile)):
            read, ref = convert_sam_record(record, self.ref[record['rname']])
            yield "{}\t{}\n".format(read, ref)
    def run(self, kmer_size=5):
        """Run karen's pipeline.
        """
        load_reference()
        outf = open(self.getLocalTempDir()+"karen", "w")
        for line in sam_to_karen():
            outf.write(line)
        outf.close()
        c = subprocess.Popen(["perl", "kmer.pl", self.readFastaFile, str(kmer_size)]).communicate()
        c = subprocess.Popen(["perl", "kmer.pl", self.referenceFastaFile, str(kmer_size)]).communicate()
        c = subprocess.Popen(["perl", "cmpKmer.pl", self.readFastaFile + "." + str(kmer_size) + "mer.fa", self.referenceFastaFile + "." + str(kmer_size) + "mer.fa"]).communicate()
        c = subprocess.Popen(["perl", "kmer_del.pl", "tmp", self.outputDir + "/kmer_del.out"]).communicate()
        c = subprocess.Popen(["perl", "kmer_ins.pl", "tmp", self.outputDir + "/kmer_ins.out"]).communicate()








