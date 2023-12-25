#!/usr/bin/perl
# Usage: pelr $0 *.misa
use warnings;
use strict;

open my $IN, "<", $ARGV[0] or die;
open my $OUT1, ">", "p1_for_residual.txt" or die;
open my $OUT2, ">", "p2_for_residual.txt" or die;
open my $OUT3, ">", "p3_for_residual.txt" or die;
open my $OUT4, ">", "p4_for_residual.txt" or die;
open my $OUT5, ">", "p5_for_residual.txt" or die;
open my $OUT6, ">", "p6_for_residual.txt" or die;

while (<$IN>) {
  chomp;
  unless ($_ =~ /ID\t/i || $_ =~ /BOUNDARY/i){
    my ($id, $nr, $type, $motif, $size, $start, $end, $location) = split("\t", $_);
    $motif =~ /\((.*)\)/;
    my ($actual_motif, $actual_motif_a) = ($1, $1);
    my $reverse_motif = $1;
    $reverse_motif =~ tr/ATGC/TACG/;
    my $reverse_motif_a = $reverse_motif = reverse $reverse_motif;
    my $red_rev;
    my $motif_len = length $actual_motif;
    for (my $i = 0; $i < $motif_len; $i ++) {
      $actual_motif =~ s/(.)(.*)/$2$1/;
      $reverse_motif =~ s/(.)(.*)/$2$1/;
      $actual_motif_a = $actual_motif if ($actual_motif lt $actual_motif_a);
      $reverse_motif_a = $reverse_motif if ($reverse_motif lt $reverse_motif_a);
      }
    if ($actual_motif_a lt $reverse_motif_a) {
      $red_rev = "$actual_motif_a/$reverse_motif_a";
    }else{
        $red_rev = "$reverse_motif_a/$actual_motif_a";
    }
    if ($motif_len == 1) {
      print $OUT1 "$_\t$actual_motif\t$red_rev\n";
    }
    if ($motif_len == 2) {
      print $OUT2 "$_\t$actual_motif\t$red_rev\n";
    }
    if ($motif_len == 3) {
      print $OUT3 "$_\t$actual_motif\t$red_rev\n";
    }
    if ($motif_len == 4) {
      print $OUT4 "$_\t$actual_motif\t$red_rev\n";
    }
    if ($motif_len == 5) {
      print $OUT5 "$_\t$actual_motif\t$red_rev\n";
    }
    if ($motif_len == 6) {
      print $OUT6 "$_\t$actual_motif\t$red_rev\n";
    }
  }elsif($_ =~ /ID\t/i){
	  print $OUT1 "$_\tactual motif\tstandardized motif\n";
	  print $OUT2 "$_\tactual motif\tstandardized motif\n";
	  print $OUT3 "$_\tactual motif\tstandardized motif\n";
	  print $OUT4 "$_\tactual motif\tstandardized motif\n";
	  print $OUT5 "$_\tactual motif\tstandardized motif\n";
	  print $OUT6 "$_\tactual motif\tstandardized motif\n";
  }
}
close $IN;
close $OUT1;
close $OUT2;
close $OUT3;
close $OUT4;
close $OUT5;
close $OUT6;
