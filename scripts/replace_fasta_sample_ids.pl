#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $print_usage = 0;

my $usage = <<USAGE;

  This script replaces fasta sample IDs (fasta headers) with an alternative ID
  provided in a separate tsv file that maps original IDS -> replacement IDs

  Writes the new fasta to stdout


  Mark Stenglein,  3/2/2022 

  Usage: $0 [-h] fasta_file tsv_map_file

   [-h]             print this message

USAGE

if ((scalar @ARGV != 2)) { print $usage and exit; }

GetOptions ("h" => \$print_usage);

my $fasta_file = shift or die("error:  no fasta file specified.\n$usage\n");
my $map_file = shift or die("error: no tsv_map_file specified.\n$usage\n");

open (my $fasta_fh, "<", $fasta_file) or die("error: could not open fasta file $fasta_file.\n$usage\n");
open (my $map_fh, "<", $map_file) or die("error: could not open tsv_map_file $map_file.\n$usage\n");

# create a hash of replacement IDS
my %ids = ();

while (<$map_fh>)
{
  chomp;
  my @fields = split "\t";
  if (scalar (@fields) != 2) 
  {
     die ("error: unexpected input in map file.  Expecting 2-column tab-delimited format.\nLine: $_\n");
  }
  # map of original ID -> replacement
  $ids{$fields[0]} = $fields[1];
}

while (<$fasta_fh>)
{
   chomp;
   # if a fasta header line
   if (/^>(.*)/)
   {
      my $id = $1;
      # if there is a replacement specified
      if ($ids{$id})
      {
         print ">$ids{$id}\n";
      }
      else
      {
	 warn ("no replacement for $id\n");
	 # no replacement - just output original
         print "$_\n";
      }
   }
   else 
   {
      # just print out all other lines
      print "$_\n";
   }
}

