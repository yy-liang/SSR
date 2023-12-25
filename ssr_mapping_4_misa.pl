#!/usr/bin/perl
# Usage: pelr $0 annot.coding.txt *.misa lncRNA.list
use warnings;
use strict;

my (%anno_start_of, %anno_end_of, %len_of, %strand_of); # key为id_type的hash
open my $IN_ANNO, "<", $ARGV[0] or die;
my $both_3UTR_5UTR = 0; # 存储同时注释到3UTR和5UTR的序列数
while (<$IN_ANNO>) {
  chomp;
  $_ =~ s/:/\t/g;
  $_ =~ s/([0-9])-/$1\t/g;
  my ($id, undef, $type,undef, undef, undef, $strand, undef, $anno_start, $anno_end, undef) = split("\t", $_, 11);
  $id =~ s/\|m\..*//;
  my $ha_key = "$id\_$type";
  $strand_of{$id} = $strand unless exists $strand_of{$id};
  #print $anno_end."\n";
  $len_of{$ha_key} = $anno_end - $anno_start + 1;
  if ($strand eq "+"){
	$anno_start_of{$ha_key} = $anno_start;
	$anno_end_of{$ha_key} = $anno_end;  
  }else{
	$anno_start_of{$ha_key} = $anno_end;
	$anno_end_of{$ha_key} = $anno_start;  
  }

  my $key_3UTR = "$id\_3UTR";
  my $key_5UTR = "$id\_5UTR";
##### both 5UTR and 3UTR #####
  if (exists $anno_start_of{$key_3UTR} && exists $anno_start_of{$key_5UTR}) {
    $both_3UTR_5UTR ++;
  }
}
close $IN_ANNO;

my @lnc = `cat $ARGV[2]`;
chomp(@lnc);

open my $IN_MISA, "<", $ARGV[1] or die;
open my $OUT1, ">", "$ARGV[1].mapping_result.txt" or die;
open my $OUT6, ">", "$ARGV[1].Coding.mapping_result.txt" or die;
open my $OUT7, ">", "$ARGV[1].lncRNA.mapping_result.txt" or die;
open my $OUT8, ">", "$ARGV[1].Others.mapping_result.txt" or die;
while (<$IN_MISA>) {
  chomp;
##### print title #####
  if ( $_ =~ m/ID\tSSR/ ) {
    print $OUT1 "#$_\tLocation\tStrand\n";
	print $OUT6 "#$_\tLocation\tStrand\n";
	print $OUT7 "#$_\tLocation\tStrand\n";
	print $OUT8 "#$_\tLocation\tStrand\n";
  }else{

##### splite misa_data #####
    my ($id, $nr, $type, $motif, $size, $start, $end) = split("\t", $_);

##### mapping #####
	if (exists $strand_of{$id} && $strand_of{$id} =~ /\+/){
	  my $id_key = "$id\_CDS";
	  if ($start >= $anno_start_of{$id_key} && $start <= $anno_end_of{$id_key} && $end <= $anno_end_of{$id_key}) {
	  print $OUT6 "$_\tCDS\t\+\n";
	  }elsif ( $start >= $anno_start_of{$id_key} &&  $start <= $anno_end_of{$id_key} && $end > $anno_end_of{$id_key}) {
	  print $OUT6 "$_\tCDS_3UTR_BOUNDARY\t\+\n";
	  }elsif ($start > $anno_end_of{$id_key}) {
	  print $OUT6 "$_\t3UTR\t\+\n";
	  }elsif ($end < $anno_start_of{$id_key}){
	  print $OUT6 "$_\t5UTR\t\+\n";
	  }elsif ($start < $anno_start_of{$id_key} && $end >= $anno_start_of{$id_key}){
	  print $OUT6 "$_\t5UTR_CDS_BOUNDARY\t\+\n";
	  }		
	}elsif(exists $strand_of{$id} && $strand_of{$id} =~ /\-/){
	  my $id_key = "$id\_CDS";
	  if ($start <= $anno_start_of{$id_key} && $start >= $anno_end_of{$id_key} && $end >= $anno_end_of{$id_key}) {
	  print $OUT6 "$_\tCDS\t\-\n";
	  }elsif ( $end <= $anno_start_of{$id_key} &&  $end >= $anno_end_of{$id_key} && $start < $anno_end_of{$id_key}) {
	  print $OUT6 "$_\tCDS_3UTR_BOUNDARY\t\-\n";
	  }elsif ($start < $anno_end_of{$id_key}) {
	  print $OUT6 "$_\t3UTR\t\-\n";
	  }elsif ($start > $anno_start_of{$id_key}){
	  print $OUT6 "$_\t5UTR\t\-\n";
	  }elsif ($start <= $anno_start_of{$id_key} && $end > $anno_start_of{$id_key}){
	  print $OUT6 "$_\t5UTR_CDS_BOUNDARY\t\-\n";
	  }		
	}elsif($id ~~ @lnc){
	  print $OUT7 "$_\tlncRNA\tNA\n";
	}else{
	  print $OUT8 "$_\tOthers\tNA\n";
	}
  }
} 
close $IN_MISA;
close $OUT1;
close $OUT6;
close $OUT7;
close $OUT8;
`cat "$ARGV[1].Coding.mapping_result.txt" "$ARGV[1].lncRNA.mapping_result.txt" "$ARGV[1].Others.mapping_result.txt" |grep -v "#" >> "$ARGV[1].mapping_result.txt"`;

##### statistics #####
open my $IN_MAPPING, "<", "$ARGV[1].mapping_result.txt" or die;
open my $OUT2, ">", "$ARGV[1].mapping_statistics.txt" or die;
#open my $OUT3, ">", "$ARGV[1].without_mononucleotied_mapping.txt" or die;
open my $OUT4, ">", "$ARGV[1].without_mononucleotied.misa" or die;
#open my $OUT5, ">", "$ARGV[1].mapping_statistics.lncRNA.txt" or die;
my ($both_3UTR_5UTR_2, $Nr_of_CDS, $Nr_of_5UTR, $Nr_of_3UTR, $Nr_of_5UTR_CDS_BOUNDARY, $Nr_of_CDS_3UTR_BOUNDARY, $Nr_of_lncRNA, $Nr_of_others) = map {0} 1..8;
my ($Nr2_of_CDS, $Nr2_of_5UTR, $Nr2_of_3UTR, $Nr2_of_5UTR_CDS_BOUNDARY, $Nr2_of_CDS_3UTR_BOUNDARY, $Nr2_of_lncRNA, $Nr2_of_others) = map {0} 1..8; # without_mononucleotied
my ($Nr_of_polyA, $Nr_of_polyT, $Nr_of_polyC, $Nr_of_polyG) = (0) x 4; # coding
my ($Nr3_of_polyA, $Nr3_of_polyT, $Nr3_of_polyC, $Nr3_of_polyG) = (0) x 4; #lncRNA
my ($Nr4_of_polyA, $Nr4_of_polyT, $Nr4_of_polyC, $Nr4_of_polyG) = (0) x 4; #others
while (<$IN_MAPPING>) {
  my @array = split("\t", $_);
  my $id = $array[0];
  my $location = $array[-2];
  my $key_3UTR = "$id\_3UTR";
  my $key_5UTR = "$id\_5UTR";
##### both 5UTR and 3UTR #####
  if (exists $anno_start_of{$key_3UTR} && exists $anno_start_of{$key_5UTR}) {
    $both_3UTR_5UTR_2 ++;
  }
##### statistics - all ssr #####
  if ($location eq "CDS"){
    $Nr_of_CDS ++;
  }
  if ($location eq "3UTR"){
    $Nr_of_3UTR ++;
  }
  if ($location eq "5UTR"){
    $Nr_of_5UTR ++;
  }
  if ($location eq "CDS_3UTR_BOUNDARY"){
    $Nr_of_CDS_3UTR_BOUNDARY ++;
  }
  if ($location eq "5UTR_CDS_BOUNDARY"){
    $Nr_of_5UTR_CDS_BOUNDARY ++;
  } 
  if ($location eq "lncRNA"){
    $Nr_of_lncRNA ++;
  }
  if ($location eq "Others"){
    $Nr_of_others ++;
  }

  if  ($_ =~ /#/)  {
#    print $OUT3 "$_";
    print $OUT4 "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\n";
  }else{
##### statistics - mononucleotide #####
    if ($_ =~ /\t\([A]\)[0-9]*\t/){
		if($location =~ /lncRNA/){
			$Nr3_of_polyA ++;
		}elsif($location =~ /Others/i){
			$Nr4_of_polyA ++;
		}else{
			$Nr_of_polyA ++;
		}    
    }elsif ($_ =~ /\t\([T]\)[0-9]*\t/){
		if($location =~ /lncRNA/){
			$Nr3_of_polyT ++;
		}elsif($location =~ /Others/i){
			$Nr4_of_polyT ++;
		}else{
			$Nr_of_polyT ++;
		}    
    }elsif ($_ =~ /\t\([G]\)[0-9]*\t/){
		if($location =~ /lncRNA/){
			$Nr3_of_polyG ++;
		}elsif($location =~ /Others/i){
			$Nr4_of_polyG ++;
		}else{
			$Nr_of_polyG ++;
		}    
    }elsif ($_ =~ /\t\([C]\)[0-9]*\t/){
		if($location =~ /lncRNA/){
			$Nr3_of_polyC ++;
		}elsif($location =~ /Others/i){
			$Nr4_of_polyC ++;
		}else{
			$Nr_of_polyC ++;
		}  
    }else{
##### statistics - without mononucleotide #####
#      print $OUT3 "$_";
      print $OUT4 "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\n";
      if ($location eq "CDS"){
        $Nr2_of_CDS ++;
      }
      if ($location eq "3UTR"){
        $Nr2_of_3UTR ++;
      }
      if ($location eq "5UTR"){
        $Nr2_of_5UTR ++;
      }
      if ($location eq "CDS_3UTR_BOUNDARY"){
        $Nr2_of_CDS_3UTR_BOUNDARY ++;
      }
      if ($location eq "5UTR_CDS_BOUNDARY"){
        $Nr2_of_5UTR_CDS_BOUNDARY ++;
      }
	  if ($location eq "lncRNA"){
        $Nr2_of_lncRNA ++;
      }
	  if ($location eq "Others"){
        $Nr2_of_others ++;
  }
    }
  }
}

my $total1 = $Nr_of_5UTR + $Nr_of_5UTR_CDS_BOUNDARY + $Nr_of_CDS + $Nr_of_CDS_3UTR_BOUNDARY + $Nr_of_3UTR;
my $total2 = $Nr2_of_5UTR + $Nr2_of_5UTR_CDS_BOUNDARY + $Nr2_of_CDS + $Nr2_of_CDS_3UTR_BOUNDARY + $Nr2_of_3UTR;
my $total3 = $Nr_of_polyA + $Nr_of_polyT + $Nr_of_polyC + $Nr_of_polyG ;

my $prop_of_5UTR = $Nr_of_5UTR / $total1 * 100;
$prop_of_5UTR = sprintf "%.2f", $prop_of_5UTR;
my $prop_of_3UTR = $Nr_of_3UTR / $total1 * 100;
$prop_of_3UTR = sprintf "%.2f", $prop_of_3UTR;
my $prop_of_CDS = $Nr_of_CDS / $total1 * 100;
$prop_of_CDS = sprintf "%.2f", $prop_of_CDS;
my $prop_of_5UTR_CDS_BOUNDARY = $Nr_of_5UTR_CDS_BOUNDARY / $total1 * 100;
$prop_of_5UTR_CDS_BOUNDARY = sprintf "%.2f", $prop_of_5UTR_CDS_BOUNDARY;
my $prop_of_CDS_3UTR_BOUNDARY = $Nr_of_CDS_3UTR_BOUNDARY / $total1 * 100;
$prop_of_CDS_3UTR_BOUNDARY = sprintf "%.2f", $prop_of_CDS_3UTR_BOUNDARY;
print $OUT2 "\n############## Coding Sequences ##############\n";
print $OUT2 "\nNo. of sequence with both 3UTR & 5UTR predicted by ANGEL: $both_3UTR_5UTR\n";
print $OUT2 "\n######## Statistics - all ########\n\n";
print $OUT2 "No. of SSR containing sequence with both 3UTR & 5UTR: $both_3UTR_5UTR_2\n\n";
print $OUT2 "No. of SSR in 5UTR: $Nr_of_5UTR ($prop_of_5UTR%)\n";
print $OUT2 "No. of SSR in 5UTR-CDS Boundary: $Nr_of_5UTR_CDS_BOUNDARY ($prop_of_5UTR_CDS_BOUNDARY%)\n";
print $OUT2 "No. of SSR in CDS: $Nr_of_CDS ($prop_of_CDS%)\n";
print $OUT2 "No. of SSR in CDS-3UTR Boundary: $Nr_of_CDS_3UTR_BOUNDARY ($prop_of_CDS_3UTR_BOUNDARY%)\n";
print $OUT2 "No. of SSR in 3UTR: $Nr_of_3UTR ($prop_of_3UTR%)\n";
print $OUT2 "Total: $total1\n";

my $prop_of_polyA = $Nr_of_polyA / $total3 * 100;
$prop_of_polyA = sprintf "%.2f", $prop_of_polyA;
my $prop_of_polyT = $Nr_of_polyT / $total3 * 100;
$prop_of_polyT = sprintf "%.2f", $prop_of_polyT;
my $prop_of_polyG = $Nr_of_polyG / $total3 * 100;
$prop_of_polyG = sprintf "%.2f", $prop_of_polyG;
my $prop_of_polyC = $Nr_of_polyC / $total3 * 100;
$prop_of_polyC = sprintf "%.2f", $prop_of_polyC;
print $OUT2 "\n######## Statistics - mononucleotide ########\n\n";
print $OUT2 "No. of poly-A: $Nr_of_polyA ($prop_of_polyA%)\n";
print $OUT2 "No. of poly-T: $Nr_of_polyT ($prop_of_polyT%)\n";
print $OUT2 "No. of poly-G: $Nr_of_polyG ($prop_of_polyG%)\n";
print $OUT2 "No. of poly-C: $Nr_of_polyC ($prop_of_polyC%)\n";
print $OUT2 "Total: $total3\n";

my $prop2_of_5UTR = $Nr2_of_5UTR / $total2 * 100;
$prop2_of_5UTR = sprintf "%.2f", $prop2_of_5UTR;
my $prop2_of_3UTR = $Nr2_of_3UTR / $total2 * 100;
$prop2_of_3UTR = sprintf "%.2f", $prop2_of_3UTR;
my $prop2_of_CDS = $Nr2_of_CDS / $total2 * 100;
$prop2_of_CDS = sprintf "%.2f", $prop2_of_CDS;
my $prop2_of_5UTR_CDS_BOUNDARY = $Nr2_of_5UTR_CDS_BOUNDARY / $total2 * 100;
$prop2_of_5UTR_CDS_BOUNDARY = sprintf "%.2f", $prop2_of_5UTR_CDS_BOUNDARY;
my $prop2_of_CDS_3UTR_BOUNDARY = $Nr2_of_CDS_3UTR_BOUNDARY / $total2 * 100;
$prop2_of_CDS_3UTR_BOUNDARY = sprintf "%.2f", $prop2_of_CDS_3UTR_BOUNDARY;
print $OUT2 "\n######## Statistics - without mononucleotide ########\n\n";
print $OUT2 "No. of SSR in 5UTR: $Nr2_of_5UTR ($prop2_of_5UTR%)\n";
print $OUT2 "No. of SSR in 5UTR-CDS Boundary: $Nr2_of_5UTR_CDS_BOUNDARY ($prop2_of_5UTR_CDS_BOUNDARY%)\n";
print $OUT2 "No. of SSR in CDS: $Nr2_of_CDS ($prop2_of_CDS%)\n";
print $OUT2 "No. of SSR in CDS-3UTR Boundary: $Nr2_of_CDS_3UTR_BOUNDARY ($prop2_of_CDS_3UTR_BOUNDARY%)\n";
print $OUT2 "No. of SSR in 3UTR: $Nr2_of_3UTR ($prop2_of_3UTR%)\n";
print $OUT2 "Total: $total2\n";

my $total4 = $Nr3_of_polyA + $Nr3_of_polyT + $Nr3_of_polyG + $Nr3_of_polyC;
my $prop3_of_polyA = $Nr3_of_polyA / $total4 * 100;
$prop3_of_polyA = sprintf "%.2f", $prop3_of_polyA;
my $prop3_of_polyT = $Nr3_of_polyT / $total4 * 100;
$prop3_of_polyT = sprintf "%.2f", $prop3_of_polyT;
my $prop3_of_polyG = $Nr3_of_polyG / $total4 * 100;
$prop3_of_polyG = sprintf "%.2f", $prop3_of_polyG;
my $prop3_of_polyC = $Nr3_of_polyC / $total4 * 100;
$prop3_of_polyC = sprintf "%.2f", $prop3_of_polyC;
print $OUT2 "\n\n############## lncRNA ##############\n";
print $OUT2 "\nTotal No. of SSR in lncRNA: $Nr_of_lncRNA\n";
print $OUT2 "\n######## Statistics - mononucleotide ########\n\n";
print $OUT2 "No. of poly-A: $Nr3_of_polyA ($prop3_of_polyA%)\n";
print $OUT2 "No. of poly-T: $Nr3_of_polyT ($prop3_of_polyT%)\n";
print $OUT2 "No. of poly-G: $Nr3_of_polyG ($prop3_of_polyG%)\n";
print $OUT2 "No. of poly-C: $Nr3_of_polyC ($prop3_of_polyC%)\n";
print $OUT2 "Total: $total4\n";
print $OUT2 "\n######## Statistics - without mononucleotide ########\n\n";
print $OUT2 "Total: $Nr2_of_lncRNA\n";

my $total5 = $Nr4_of_polyA + $Nr4_of_polyT + $Nr4_of_polyG + $Nr4_of_polyC;
my $prop4_of_polyA = $Nr4_of_polyA / $total5 * 100;
$prop4_of_polyA = sprintf "%.2f", $prop4_of_polyA;
my $prop4_of_polyT = $Nr4_of_polyT / $total5 * 100;
$prop4_of_polyT = sprintf "%.2f", $prop4_of_polyT;
my $prop4_of_polyG = $Nr4_of_polyG / $total5 * 100;
$prop4_of_polyG = sprintf "%.2f", $prop4_of_polyG;
my $prop4_of_polyC = $Nr4_of_polyC / $total5 * 100;
$prop4_of_polyC = sprintf "%.2f", $prop4_of_polyC;
print $OUT2 "\n\n############## Other Sequences ##############\n";
print $OUT2 "\nTotal No. of SSR in Other Sequences: $Nr_of_others\n";
print $OUT2 "\n######## Statistics - mononucleotide ########\n\n";
print $OUT2 "No. of poly-A: $Nr4_of_polyA ($prop4_of_polyA%)\n";
print $OUT2 "No. of poly-T: $Nr4_of_polyT ($prop4_of_polyT%)\n";
print $OUT2 "No. of poly-G: $Nr4_of_polyG ($prop4_of_polyG%)\n";
print $OUT2 "No. of poly-C: $Nr4_of_polyC ($prop4_of_polyC%)\n";
print $OUT2 "Total: $total5\n";
print $OUT2 "\n######## Statistics - without mononucleotide ########\n\n";
print $OUT2 "Total: $Nr2_of_others\n";

close $IN_MAPPING;
close $OUT2;
#close $OUT3;
close $OUT4;
#close $OUT5;
