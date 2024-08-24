#!usr/bin/perl -w
use strict;
if(@ARGV!=1)    ## $ARGV[0] only  misa.pl result  file
{
	print "\nUsage:\n\tperl ssr_statistic.pl Gourd-Unigene.fa.cut.statistics\n\n";
	print "\nInput file is program(misa.pl) SSR statistic result:\n\tlike this:\n\tGourd-Unigene.fa.cut.statistics\n\n";
	print "The result file will like this: \n\tGourd-Unigene.fa.cut.statistics.totality.txt\n\tGourd-Unigene.fa.cut.statistics.classify.txt\n\tGourd-Unigene.fa.cut.statistics.drawSVG.txt\n\n\n";
	die;
}
#open(IN,"Gourd-Unigene.fa.statistics")||die"cannot open:$!";
my $misa_statistics_file=$ARGV[0];

open(IN,$misa_statistics_file)||die"cannot open:$!";
open(OUT1,">$misa_statistics_file.totality.txt")||die"cannot open:$!";
open(OUT2,">$misa_statistics_file.classify.txt")||die"cannot open:$!";
open(OUT3,">$misa_statistics_file.drawSVG.txt")||die"cannot open:$!";

my %hash=();
my %most15=();    my %most15_num=();
my @type_order=();
my $mark=0;
my $ssr_total=0;
my $min_ssr_repeat=undef;
my $max_ssr_repeat=undef;
my @array=();
while(<IN>)
{
	chomp;
####  totality.txt
	print OUT1 "$1\t$2\n" if(/(^Total number of sequences examined:)\s+(\d+)/);
	print OUT1 "$1\t$2\n" if(/(^Total size of examined sequences \(bp\):)\s+(\d+)/);
	print OUT1 "$1\t$2\n" if(/(^Total number of identified SSRs:)\s+(\d+)/);
	print OUT1 "$1\t$2\n" if(/(^Number of SSR containing sequences:)\s+(\d+)/);
	print OUT1 "$1\t$2\n" if(/(^Number of sequences containing more than 1 SSR:)\s+(\d+)/);
	print OUT1 "$1\t$2\n" if(/(^Number of SSRs present in compound formation:)\s+(\d+)/);
	if(/^Unit size\tNumber of SSRs/)
	{
		$mark=1;
		print OUT2 "Number of repeat unit\t";
	}
	if(/^1\t(\d+)/ && $mark==1)
	{
		print OUT1 "Mono-nucleotide\t$1\n";
		print OUT2 "Mono-\t";
	}
	if(/^2\t(\d+)/ && $mark==1)
	{
		print OUT1 "Di-nucleotide\t$1\n";
		print OUT2 "Di-\t";
	}
	if(/^3\t(\d+)/ && $mark==1)
	{
		print OUT1 "Tri-nucleotide\t$1\n";
		print OUT2 "Tri-\t";
	}
	if(/^4\t(\d+)/ && $mark==1)
	{
		print OUT1 "Tetra-nucleotide\t$1\n";
		print OUT2 "Tetra-\t";
	}
	if(/^5\t(\d+)/ && $mark==1)
	{
		print OUT1 "Penta-nucleotide\t$1\n";
		print OUT2 "Penta-\t";
	}
	if(/^6\t(\d+)/ && $mark==1)
	{
		print OUT1 "Hexa-nucleotide\t$1\n";
		print OUT2 "Hexa\n";
		$mark=0;
	}
####  classity.txt
	$mark=2 if(/^Frequency of classified repeat types/);
	if(/^Repeats/ && $mark==2)
	{
		$mark=3;  @array=split(/\s+/,$_);
		$min_ssr_repeat=$array[1];   $max_ssr_repeat=$min_ssr_repeat+10;
		next;
	}
	if($mark==3)
	{
		my @a=split(/\t/,$_);    my @b=split(/\//,$a[0]);   my $rep=length $b[0];
		push @type_order,$a[0];    $most15{$a[0]}=$a[-1];    $ssr_total+=$a[-1];
		if(exists $most15_num{$a[-1]})
		{
			$most15_num{$a[-1]}.="\n".$a[0];
		}
		else
		{
			$most15_num{$a[-1]}=$a[0];
		}
		foreach my $out (1 .. $#a-1)
		{
			if($array[$out] <= $max_ssr_repeat && $array[$out] >= $min_ssr_repeat)
			{
				if($a[$out]=~/\d+/)
				{
					$hash{$array[$out]}{$rep}+=$a[$out];
				}
				else
				{
					$hash{$array[$out]}{$rep}+=0;
				}
			}
			if($array[$out] > $max_ssr_repeat)
			{
				if($a[$out]=~/\d+/)
				{
					$hash{$max_ssr_repeat+1}{$rep}+=$a[$out];
				}
				else
				{
					$hash{$max_ssr_repeat+1}{$rep}+=0;
				}
			}
			if($array[$out] < $min_ssr_repeat)
			{
				if($a[$out]=~/\d+/)
				{
					$hash{$min_ssr_repeat-1}{$rep}+=$a[$out];
				}
				else
				{
					$hash{$min_ssr_repeat-1}{$rep}+=0
				}
			}
		}
	}
}
close IN;
foreach my $out(sort {$a<=>$b} keys %hash)
{
	if($out <= $max_ssr_repeat && $out >= $min_ssr_repeat)
	{
		print OUT2 "$out\t";
	}
	if($out > $max_ssr_repeat)
	{
		print OUT2 ">=$out\t";
	}
	if($out < $min_ssr_repeat)
	{
		print OUT2 "<=$out\t";
	}
	foreach my $in(sort {$a<=>$b} keys %{$hash{$out}})
	{
		print OUT2 "$hash{$out}{$in}\t";
	}
	print OUT2 "\n";
}

my %store=();    my $count=0;
foreach my $out (sort {$b<=>$a} keys %most15_num)
{
	my @a=split(/\n/,$most15_num{$out});
	foreach my $in (@a)
	{
		$count++;
		if($count<=15)
		{
			$store{$in}=$out;
		}
		else
		{
			last;
		}
	}
	if($count>15)
	{
		last;
	}
}
my $most15_total=0;
foreach my $out (@type_order)
{
	if(exists $store{$out})
	{
		my $ratio=$most15{$out}/$ssr_total;
		print OUT3 "$out\t$most15{$out}\t$ratio\n";
		$most15_total+=$most15{$out};
	}
}
print OUT3 "others\t",$ssr_total-$most15_total,"\t",($ssr_total-$most15_total)/$ssr_total,"\n";

close OUT1;
close OUT2;
close OUT3;
