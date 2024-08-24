#! perl
# Usage: perl $0 *.misa
use strict;
use warnings;
use List::Util qw /max sum/;

open my $IN, "<", $ARGV[0] or die;
my(%mono, %di, %tri, %tetra, %penta, %hexa);
my($n, $sum_len, $sum_rep, $mean_len, $mean_rep);
open my $OUT1, ">", $ARGV[0].".distribution.xls" or die;
open my $OUT2, ">", $ARGV[0].".p1.distribution.xls" or die;
open my $OUT4, ">", $ARGV[0].".p2.distribution.xls" or die;
open my $OUT5, ">", $ARGV[0].".p3.distribution.xls" or die;
open my $OUT6, ">", $ARGV[0].".p4.distribution.xls" or die;
open my $OUT7, ">", $ARGV[0].".p5.distribution.xls" or die;
open my $OUT8, ">", $ARGV[0].".p6.distribution.xls" or die;

while(<$IN>){
	unless($_ =~ /^ID/ or $_ =~ /^#/){
		chomp;
		$n ++;
		my (undef, undef, $type, $motif, undef) = split("\t", $_, 5);
		#print "$type\n";
		my ($length, $size);
		if ($motif =~ m/\((.*)\)(.*)$/){
			$length = length($1) * $2;
			$size = $2;
			#print "$length\t$size\n";
		}
		print $OUT1 "$_\t$length\t$size\n";
		$sum_len += $length;
		$sum_rep += $size;
		if($type =~ m/p1/){
			#print "mono\n";
			$mono{$size} ++;
			print $OUT2 "$_\t$length\t$size\n";
		}elsif($type =~ m/p2/){
			$di{$size} ++;
			print $OUT4 "$_\t$length\t$size\n";
		}elsif($type =~ m/p3/){
			$tri{$size} ++;
			print $OUT5 "$_\t$length\t$size\n";
		}elsif($type =~ m/p4/){
			$tetra{$size} ++;
			print $OUT6 "$_\t$length\t$size\n";
		}elsif($type =~ m/p5/){
			$penta{$size} ++;
			print $OUT7 "$_\t$length\t$size\n";
		}elsif($type =~ m/p6/){
			$hexa{$size} ++;
			print $OUT8 "$_\t$length\t$size\n";
		}
	}else{
		chomp;
		print $OUT1 "$_\tlength\ttandem_repeat\n";
		print $OUT2 "$_\tlength\ttandem_repeat\n";
		print $OUT4 "$_\tlength\ttandem_repeat\n";
		print $OUT5 "$_\tlength\ttandem_repeat\n";
		print $OUT6 "$_\tlength\ttandem_repeat\n";
		print $OUT7 "$_\tlength\ttandem_repeat\n";
		print $OUT8 "$_\tlength\ttandem_repeat\n";
	}
}
close $IN;
close $OUT1;
print "$n\n";
$mean_len = $sum_len / $n;
$mean_rep = $sum_rep / $n;

open my $OUT3, ">", $ARGV[0].".TandemRepeat" or die;
print $OUT3 "Average of SSR length: $mean_len\n";
print $OUT3 "Average of Tandem Repeats: $mean_rep\n";

print $OUT3 "\n\tMono-\tDi-\tTri-\tTetra-\tPenta-\tHexa-\n";

my $max = max(keys %mono, keys %di, keys %tri, keys %tetra, keys %penta, keys %hexa);
for(my $i = 3; $i <= $max; $i ++){
	$mono{$i} = 0 unless $mono{$i};
	$di{$i} = 0 unless $di{$i};
	$tri{$i} = 0 unless $tri{$i};
	$tetra{$i} = 0 unless $tetra{$i};
	$penta{$i} = 0 unless $penta{$i};
	$hexa{$i} = 0 unless $hexa{$i};
	print $OUT3 "$i\t$mono{$i}\t$di{$i}\t$tri{$i}\t$tetra{$i}\t$penta{$i}\t$hexa{$i}\n";
}

my ($mono, $di, $tri, $tetra, $penta, $hexa);
for(keys %mono){
	$mono += $mono{$_} if $_ > 29;
}
for(keys %di){
	$di += $di{$_} if $_ > 29;
}
for(keys %tri){
	$tri += $tri{$_} if $_ > 29;
}
for(keys %tetra){
	$tetra += $tetra{$_} if $_ > 29;
}
for(keys %penta){
	$penta += $penta{$_} if $_ > 29;
}
for(keys %hexa){
	$hexa += $hexa{$_} if $_ > 29;
}
$mono = 0 unless $mono;
$di = 0 unless $di;
$tri = 0 unless $tri;
$tetra = 0 unless $tetra;
$penta = 0 unless $penta;
$hexa = 0 unless $hexa;
print $OUT3 ">=30\t$mono\t$di\t$tri\t$tetra\t$penta\t$hexa\n";
close $OUT3;