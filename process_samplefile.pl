#!/usr/bin/perl -w

use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];

open OUT,">$out";
open IN,"<$in";
while(<IN>)
{
	chomp;
	if($_ =~ /^sample.+/)
	{
		print OUT "$_\n";
	}
	else
	{
		my @arr = split "\t",$_;
		print OUT "$arr[0]\t";
		foreach my $i(1..$#arr)
		{
			if($i != $#arr)
			{
				if($arr[$i] == 1)
				{
					print OUT "$arr[0]\t";
				}
				else
				{
					print OUT "0\t";
				}
				#print "\n";
			}
			else
			{
				if($arr[$i] == 1)
                                {
                                        print OUT "$arr[0]";
                                }
                                else
                                {
                                        print OUT "0";
                                }
                                #print "\n";

			}
		}
		print OUT "\n"
			
	}
}
close IN;
