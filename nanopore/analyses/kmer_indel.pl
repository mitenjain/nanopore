#! /usr/bin/perl -w
use strict;
die "usage: kmer_indel.pl <last2fa> <ref.fa> <del_out.txt> <ins_out.txt> <kmer_size>\n" unless($ARGV[0] && $ARGV[1] && $ARGV[2] && $ARGV[3] && $ARGV[4]);

open(FILE,"$ARGV[0]");
open(REF,"$ARGV[1]");
open(DEL,">$ARGV[2]");
open(INS,">$ARGV[3]");

my $kmerSize=$ARGV[4];
##########
# Read in all possible 5-mers (f/r orientation) found in a given reference and determine frequency
##########
my $refSeq;
while(<REF>){
    chomp;
    next if($_=~ />\S+/);
    $refSeq.=$_;
}
close REF;
my $refFreq=kmer($refSeq);
my %refSeq_kmerFreq=%$refFreq;
# resulting hash contains primary keys = a given 5-mer, secondary key ='c' (count) or 'f' (frequency), and values are the count or the normalized frequency 
##########



##########
# Read in the alignment data for the reads
##########

my %del;
my %ins;

my $totGapDel=0;
my $totGapIns=0;

while(<FILE>){
    chomp;
 
    my @l=split(/\t+/,$_);
    my $l1=length($l[0]);
    my $l2=length($l[1]);

    my $start = 0;
    my $slide = 1;
    my $windowSize=$kmerSize;
    my $end = $start + $windowSize;
    my $o=0;

    my $s1=$l[0];
    my $s2=$l[1];
    while ($end <= $l1) {
        my $oligo1 = substr($s1,$start,$windowSize);
	my $oligo2 = substr($s2,$start,$windowSize);
	if(($oligo2=~ /\-/)&&($oligo1!~ /\-/)){ #gap in the read and not in the reference/deletion
	    $totGapDel+=2;
	    my $rc=reverse($oligo1);
	    $rc=~ tr/ATCGN/TAGCN/;
	    $del{$oligo1}=1 unless($del{$oligo1}++);
	    $del{$rc}=1 unless ($del{$rc}++);
	}
	if(($oligo2!~ /\-/)&&($oligo1=~ /\-/)){ #gap in the reference and not in the read/insertion
            $totGapIns+=2;
            my $rc=reverse($oligo2);
            $rc=~ tr/ATCGN/TAGCN/;
            $ins{$oligo2}=1 unless($ins{$oligo2}++);
            $ins{$rc}=1 unless ($ins{$rc}++);
        }
        $start += $slide;
        $end = $start + $slide;
    }
}
close FILE;
print DEL "5-mer\tdelCount\tdelFreq\ttotalCount\ttotalFreq\tlogTratio\n";
foreach my $d(sort{$del{$b}<=>$del{$a}} keys %del){
    my $dFreq=$del{$d}/$totGapDel;
    if(exists $refSeq_kmerFreq{$d}){
	my $logT=log($dFreq/$refSeq_kmerFreq{$d}{'f'});
	print DEL "$d\t$del{$d}\t$dFreq\t$refSeq_kmerFreq{$d}{'c'}\t$refSeq_kmerFreq{$d}{'f'}\t$logT\n";
    }
}
close DEL;

print INS "5-mer\tinsCount\tinsFreq\ttotalCount\ttotalFreq\tlogTratio\n";
foreach my $i(sort{$ins{$b}<=>$ins{$a}}keys %ins){
    my $iFreq=$ins{$i}/$totGapIns;
    if(exists $refSeq_kmerFreq{$i}){
        my $logT=log($iFreq/$refSeq_kmerFreq{$i}{'f'});
        print INS "$i\t$ins{$i}\t$iFreq\t$refSeq_kmerFreq{$i}{'c'}\t$refSeq_kmerFreq{$i}{'f'}\t$logT\n";
    }
}
close INS;

sub kmer{
    my ($s)=@_;
    my $len=length($s);
    my $start = 0;
    my $slide = 1;
    my $windowSize=5;
    my $end = $start + $windowSize;
    my $o=0;
    my %freq;
    my %count;
    my $tot=0;

    while ($end <= $len) {
	my $oligo = substr($s,$start,$windowSize);
        $start += $slide;
        $end = $start + $slide;
        if($oligo!~ /N/){
            my $oL=length($oligo);
            if($oL ==$windowSize){
		$tot+=2;
                my $rc=reverse($oligo);
                $rc=~ tr/ATCGN/TAGCN/;
                $count{$oligo}=1 unless($count{$oligo}++);
		$count{$rc}=1 unless($count{$rc}++);
            }
        }
    }
    foreach my $cKmer(keys %count){
	my $freq=$count{$cKmer}/$tot;
	$freq{$cKmer}{'f'}=$freq;
	$freq{$cKmer}{'c'}=$count{$cKmer};
    }
    my $ref=\%freq;
    return($ref);
}
