#!/usr/bin/env python

import sys, os, pysam, argparse, array
from subprocess import Popen, PIPE

def parse_args(args):
	parser = argparse.ArgumentParser()
	parser.add_argument("--template_sam", type=str, help="samfile from template alignments")
	parser.add_argument("--twoD_sam", type=str, help="samfile from 2d alignments")
	parser.add_argument("--complement_sam", type=str, help="samfile from complement alignments")
	parser.add_argument("--reference", type=argparse.FileType("r"), help="matching reference fasta")
	parser.add_argument("--template", type=argparse.FileType("r"), help="matching template fastq")
	parser.add_argument("--complement", type=argparse.FileType("r"), help="matching complement fastq")
	parser.add_argument("--output", type=argparse.FileType("w"), help="output file")
	parser.add_argument("--muscle", type=argparse.FileType("w"), help="uncounted muscle output")
	return parser.parse_args()

def fastaRead(fileHandle):
    """iteratively a sequence for each '>' it encounters, ignores '#' lines
    """
    line = fileHandle.readline()
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = array.array('c')
            while line != '' and line[0] != '>':
                if line[0] != '#':
                    seq.extend([ i for i in line[:-1] if not i.isspace() ]) #The white-space check is to remove any annoying trailing characters.
                line = fileHandle.readline()
            for i in seq:
                #For safety and sanity I only allows roman alphabet characters in fasta sequences.
                if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'):
                    raise RuntimeError("Invalid FASTA character, ASCII code = \'%d\', found in input sequence %s" % (ord(i), name))
            yield name, seq.tostring()
        else:
            line = fileHandle.readline()

def fastqReadNoQual(fileHandle):
    """Reads a fastq file iteratively - modified to not return quality scores
    """
    line = fileHandle.readline()
    while line != '':
        if line[0] == '@':
            name = line[1:-1]
            seq = fileHandle.readline()[:-1]
            plus = fileHandle.readline()
            if plus[0] != '+':
                raise RuntimeError("Got unexpected line: %s" % plus)
            qualValues = [ ord(i) for i in fileHandle.readline()[:-1] ]
            for i in seq:
                #For safety and sanity I only allows roman alphabet characters in fasta sequences.
                if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'):
                    raise RuntimeError("Invalid FASTQ character, ASCII code = \'%d\', found in input sequence %s" % (ord(i), name))
            yield name, seq
        line = fileHandle.readline()


def call_muscle(name1, seq1, name2, seq2):
	p = Popen(["muscle", "-quiet"], stdout=PIPE, stdin=PIPE).communicate(">{}\n{}\n>{}\n{}\n".format(name1, seq1, name2, seq2))[0].split("\n")
	return (p[1], p[3])

def calculate_identity(align):
	total, match = 0.0, 0.0
	for ref_base, read_base in align:
		if ref_base == read_base and ref_base != "-" and ref_base != "N":
			total += 1; match += 1
		else:
			total += 1
	return match/total


def main(args):
	args = parse_args(args)
	#load all seqs as name, seq dicts into ram because why not
	#this assumes that the naming is set up properly so each file has no duplicates
	template_fastq_dict = dict([x for x in fastqReadNoQual(args.template)])
	complement_fastq_dict = dict([x for x in fastqReadNoQual(args.complement)])
	template_record_dict = {x.qname: x for x in pysam.Samfile(args.template_sam)}
	complement_record_dict = {x.qname: x for x in pysam.Samfile(args.complement_sam)}

	reference = dict([[x[0].split()[0], x[1]] for x in fastaRead(args.reference)])
	samfile = pysam.Samfile(args.twoD_sam)
	aligned = dict()
	for twoD_record in samfile:
		if twoD_record.is_unmapped is False and twoD_record.qname not in template_record_dict and twoD_record.qname not in complement_record_dict:
			#only look at read sets mapped in 2d but not in either complement or template
			#pull down reference region that the 2d mapped too
			ref_stop, ref_start, ref_name = int(twoD_record.aend), int(twoD_record.aend - twoD_record.alen), samfile.getrname(twoD_record.tid)
			ref_seq = reference[ref_name][ref_start:ref_stop]
			template_seq = template_fastq_dict[twoD_record.qname]
			complement_seq = template_fastq_dict[twoD_record.qname]
			aligned_template = call_muscle(ref_name, ref_seq, "template_ " + twoD_record.qname, template_seq)
			aligned_complement = call_muscle(ref_name, ref_seq, "template_ " + twoD_record.qname, complement_seq)
			args.muscle.write("\t".join(call_muscle(ref_name, ref_seq, "template_ " + twoD_record.qname, template_seq))); args.muscle.write("\n")
			args.muscle.write("\t".join(call_muscle(ref_name, ref_seq, "complement_" + twoD_record.qname, complement_seq))); args.muscle.write("\n")
			aligned[twoD_record.qname] = (aligned_template, aligned_complement)
	args.output.write("ReadName\tTemplateIdentity\tComplementIdentity\n")
	for name, (aligned_template, aligned_complement) in aligned.iteritems():
		args.output.write("{}\t{}\t{}\n".format(name, calculate_identity(aligned_template), calculate_identiy(aligned_complement)))



if __name__ == "__main__":
    sys.exit(main(sys.argv))
