#! /usr/bin/perl -w
use strict;
my $t=0;
my $seq;
my $c;
my $h;

die "must provide input file\nusage: kmer.pl fastaFile k-merSize\n" unless ($ARGV[0]);
die "must provide kmer size\nusage: kmer.pl fastaFile k-merSize\n" unless($ARGV[1]);

my $outFile=$ARGV[2];
open(OUT,">$outFile");
open(FA,"$ARGV[0]");
while(<FA>){
    chomp;
    if($_=~ />(\S+)/){
	my $tmpH=$1;
	if($t==0){
	    $t=1;
	    $h=$tmpH;
	}
	else{
	    kmer($seq,$h);
	    $seq='';
	    $h=$tmpH;
	}
    }
    else{
	$seq.=uc($_);
    }
}
kmer($seq,$h);

sub kmer{
    my ($s,$header)=@_;
    my $len=length($s); 
    my $start = 0;
    my $slide = 1; 
    my $windowSize=12;
    my $end = $start + $windowSize;
    my $o=0;    
    
    while ($end <= $len) {
	my $oligo = substr($s,$start,$windowSize);
	$start += $slide;
	$end = $start + $slide;
	if($oligo!~ /N/){
	    my $oL=length($oligo);
	    if($oL ==$windowSize){
		my $rc=reverse($oligo);
		$rc=~ tr/ATCGN/TAGCN/;
		print OUT "$header\t$start\t$oligo\t$rc\n";
	    }
	}
    }
}
close OUT;
