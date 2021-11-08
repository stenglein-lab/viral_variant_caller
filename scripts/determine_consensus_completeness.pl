#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $print_usage = 0;

my $usage = <<USAGE;

  This script calculates the % completness of provided fasta format sequence(s), 
  defined as the fraction of the sequences that are not Ns, summed over all provided
  sequences.

  Reads fasta-formatted sequences from standard input

  Writes the fraction of non-N bases in sequences to standard output


  Mark Stenglein,  11/4/2021

  Usage: $0 [-h] fasta_formatted_sequences

   [-h]             print this message

USAGE

if ((scalar @ARGV == 0) and -t STDIN) { print $usage and exit; }

GetOptions ("h" => \$print_usage);

my $total_count = 0;
my $N_count = 0;

while (<>)
{
   chomp;

	# ignore fasta header lines
	if (/^>/) { next; }

	# all other lines assumed to be sequence data
	# convert to upper case in case not
	my $seq = uc($_);

	for (my $i = 0; $i < length ($seq); $i++)
	{
	   my $c = substr($seq, $i, 1);
	   if ($c eq "N")
		{
			# warn "$c == N\n";
		   $N_count += 1;
		}
		$total_count += 1;
	}
}

if ($total_count > 0)
{
   my $fraction_n = ($total_count - $N_count) / $total_count;
   my $fraction_n_output = sprintf("%0.3f", $fraction_n);
   print "$fraction_n_output";
}
else
{
   print "0";
}


