#!/usr/bin/perl
# Usage: pelr $0 *.misa
use warnings;
use strict;

open my $IN, "<", $ARGV[0] or die;
open my $OUT, ">", "standardized_".$ARGV[0] or die;

print $OUT "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\tLocation\tactual motif\tstandardized motif 1\tstandardized motif 2\trepeat time\n";

while (<$IN>) {
  chomp;
  unless ($_ =~ /ID\t/i || $_ =~ /BOUNDARY/i){
    my ($id, $nr, $type, $motif, $size, $start, $end, $location) = split("\t", $_);
	$motif =~ /\((.*)\)([0-9]*)/;
    my ($actual_motif, $actual_motif_a) = ($1, $1);
	my $repeat = $2;
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
	my $red_rev2 = $red_rev;
	$red_rev2 =~ s/\/.*//i;
	print $OUT "$_\t$actual_motif\t$red_rev2\t$red_rev\t$repeat\n"
  }
}
