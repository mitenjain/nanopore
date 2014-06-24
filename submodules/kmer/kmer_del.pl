#! /usr/bin/perl -w
use strict;

die "must provide last2fasta file and cmpKmer.out file\nusage kmer_del.pl last2faFile cmpKmer.out\ncmpKmerout is needed to obtain the list of all possible kmers for randomized null\n" unless($ARGV[0]);
die "must provide last2fasta file and cmpKmer.out file\nusage kmer_del.pl last2faFile cmpKmer.out\ncmpKmerout is needed to obtain the list of all possible kmers for randomized null\n" unless($ARGV[1]);

open(FILE,"$ARGV[0]");
open(RAND,"$ARGV[1]");
my @rand;

open(DEL,">cmpKmer.del.txt");
my $kmerSize;
while(<RAND>){
    chomp;
    my @l=split(/\t+/,$_);
    $kmerSize=length($l[0]);
    push(@rand,$l[0]);
}
close RAND;

my %rand;
my %del;
my %ins;
my $tot=0;
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
	if(($oligo2=~ /\-/)&&($oligo1!~ /\-/)){
	    my $randVal=rand(@rand);
	    $rand{$rand[$randVal]}=1 unless ($rand{$rand[$randVal]}++);
	    $del{$oligo1}=1 unless ($del{$oligo1}++);
	    $tot++;
	}
        $start += $slide;
        $end = $start + $slide;
	#print\t$start\t$oligo\t$rc\n";
    }
}
close FILE;
foreach my $d(keys %del){
    my $dFreq=$del{$d}/$tot;
    if(exists $rand{$d}){
	my $rFreq=$rand{$d}/$tot;
	my $logT=log($dFreq/$rFreq);
	print DEL "$d\t$del{$d}\t$dFreq\t$rand{$d}\t$rFreq\t$logT\n";
    }
}
close DEL;

