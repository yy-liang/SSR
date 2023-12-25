#!/usr/bin/perl
# Usage: pelr $0 *.misa.mapping_results.pl
use warnings;
use strict;

open my $IN, "<", $ARGV[0] or die;
open my $OUT1, ">", "p1.mapping_results.txt" or die;
open my $OUT2, ">", "p2.mapping_results.txt" or die;
open my $OUT3, ">", "p3.mapping_results.txt" or die;
open my $OUT4, ">", "p4.mapping_results.txt" or die;
open my $OUT5, ">", "p5.mapping_results.txt" or die;
open my $OUT6, ">", "p6.mapping_results.txt" or die;
open my $OUT7, ">", "motif_statistics.$ARGV[0].txt" or die;

my ($p1_CDS, $p1_5UTR, $p1_3UTR, $p1_5UTR_CDS_B, $p1_CDS_3UTR_B) = qw/0 0 0 0 0/;
my ($p2_CDS, $p2_5UTR, $p2_3UTR, $p2_5UTR_CDS_B, $p2_CDS_3UTR_B) = qw/0 0 0 0 0/;
my ($p3_CDS, $p3_5UTR, $p3_3UTR, $p3_5UTR_CDS_B, $p3_CDS_3UTR_B) = qw/0 0 0 0 0/;
my ($p4_CDS, $p4_5UTR, $p4_3UTR, $p4_5UTR_CDS_B, $p4_CDS_3UTR_B) = qw/0 0 0 0 0/;
my ($p5_CDS, $p5_5UTR, $p5_3UTR, $p5_5UTR_CDS_B, $p5_CDS_3UTR_B) = qw/0 0 0 0 0/;
my ($p6_CDS, $p6_5UTR, $p6_3UTR, $p6_5UTR_CDS_B, $p6_CDS_3UTR_B) = qw/0 0 0 0 0/;
while (<$IN>) {
  chomp;
  my ($id, $nr, $type, $motif, $size, $start, $end, $location) = split("\t", $_);
##### p1 #####
  if ($type =~ /p1/){
    print $OUT1 "$_\n";
    if ($location =~ /CDS$/) {
      $p1_CDS ++;
    }elsif ($location =~ /5UTR$/) {
      $p1_5UTR ++;
    }elsif ($location =~ /3UTR$/) {
      $p1_3UTR ++;
    }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
      $p1_5UTR_CDS_B ++;
    }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
      $p1_CDS_3UTR_B ++;
    }
  }
##### p2 #####
 if ($type =~ /p2/){
   print $OUT2 "$_\n";
   if ($location =~ /CDS$/) {
     $p2_CDS ++;
   }elsif ($location =~ /5UTR$/) {
     $p2_5UTR ++;
   }elsif ($location =~ /3UTR$/) {
     $p2_3UTR ++;
   }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
     $p2_5UTR_CDS_B ++;
   }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
     $p2_CDS_3UTR_B ++;
   }
 }
##### p3 #####
 if ($type =~ /p3/){
   print $OUT3 "$_\n";
   if ($location =~ /CDS$/) {
     $p3_CDS ++;
   }elsif ($location =~ /5UTR$/) {
     $p3_5UTR ++;
   }elsif ($location =~ /3UTR$/) {
     $p3_3UTR ++;
   }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
     $p3_5UTR_CDS_B ++;
   }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
     $p3_CDS_3UTR_B ++;
   }
 }
##### p4 #####
 if ($type =~ /p4/){
   print $OUT4 "$_\n";
   if ($location =~ /CDS$/) {
     $p4_CDS ++;
   }elsif ($location =~ /5UTR$/) {
     $p4_5UTR ++;
   }elsif ($location =~ /3UTR$/) {
     $p4_3UTR ++;
   }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
     $p4_5UTR_CDS_B ++;
   }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
     $p4_CDS_3UTR_B ++;
   }
 }
##### p5 #####
 if ($type =~ /p5/){
   print $OUT5 "$_\n";
   if ($location =~ /CDS$/) {
     $p5_CDS ++;
   }elsif ($location =~ /5UTR$/) {
     $p5_5UTR ++;
   }elsif ($location =~ /3UTR$/) {
     $p5_3UTR ++;
   }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
     $p5_5UTR_CDS_B ++;
   }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
     $p5_CDS_3UTR_B ++;
   }
 }
##### p6#####
 if ($type =~ /p6/){
   print $OUT6 "$_\n";
   if ($location =~ /CDS$/) {
     $p6_CDS ++;
   }elsif ($location =~ /5UTR$/) {
     $p6_5UTR ++;
   }elsif ($location =~ /3UTR$/) {
     $p6_3UTR ++;
   }elsif ($location =~ /5UTR\_CDS\_BOUNDARY/) {
     $p6_5UTR_CDS_B ++;
   }elsif ($location =~ /CDS\_3UTR\_BOUNDARY/) {
     $p6_CDS_3UTR_B ++;
   }
 }
}

##### print ######
print $OUT7 "\n############### statistics results ###############\n\n";
print $OUT7 " \t5UTR\tCDS\t3UTR\t5UTR_CDS_BOUNDARY\tCDS_3UTR_BOUNDARY\n";
print $OUT7 "p1\t$p1_5UTR\t$p1_CDS\t$p1_3UTR\t$p1_5UTR_CDS_B\t$p1_CDS_3UTR_B\n";
print $OUT7 "p2\t$p2_5UTR\t$p2_CDS\t$p2_3UTR\t$p2_5UTR_CDS_B\t$p2_CDS_3UTR_B\n";
print $OUT7 "p3\t$p3_5UTR\t$p3_CDS\t$p3_3UTR\t$p3_5UTR_CDS_B\t$p3_CDS_3UTR_B\n";
print $OUT7 "p4\t$p4_5UTR\t$p4_CDS\t$p4_3UTR\t$p4_5UTR_CDS_B\t$p4_CDS_3UTR_B\n";
print $OUT7 "p5\t$p5_5UTR\t$p5_CDS\t$p5_3UTR\t$p5_5UTR_CDS_B\t$p5_CDS_3UTR_B\n";
print $OUT7 "p6\t$p6_5UTR\t$p6_CDS\t$p6_3UTR\t$p6_5UTR_CDS_B\t$p6_CDS_3UTR_B\n";
