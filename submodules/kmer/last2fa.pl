#! /usr/bin/perl -w
use strict;

die "provide LAST tab-delimited output file\nusage: last2fa.pl LASToutput\n" unless($ARGV[0]);
open(FILE,"$ARGV[0]");
my $outFile=$ARGV[0] . '.last2fa.txt';
open(CMP,">$outFile");
my $c=0;
my $bi;
my $ch;
while(<FILE>){
    chomp;
    next if($_=~ /^\#/);
    if($_=~ /^s/){
	if($_=~ /channel/){
	    my @l=split(/\s+/,$_);
	    $ch=$l[6];
	    print CMP "$bi\t$ch\n";
	}
	else{
	    my @l=split(/\s+/,$_);
	    $bi=$l[6];
	}
    }
}
close FILE;
close CMP;
