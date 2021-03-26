#!/usr/bin/env perl
# rankedlineage_keep_non_eukaryota_taxa.pl
# 2020_07_08
use warnings;
use strict;

# This script takes rankedlineage.dmp and splits in into Eukaryota and non-Eukaryota
 my $out_non = "rankedlineage_non_euk.txt";
 my $out_euk = "rankedlineage_euk.txt";

 # Read in taxonomy and keep only those that are Archaea, Bacteria or Viruses
open(IN, "rankedlineage.dmp") or die "Couldn't find the taxonomy file $!";
open (OUT, ">$out_non") or die "Couldn't create a file for your non-Eukaryota taxonomy $!";
open (OUT_2, ">$out_euk") or die "Couldn't create a file for your Eukaryota taxonomy $!";
while (my $line = <IN>) {
  chomp($line);
  my($tax_id, $b1, $tax_name, $b2, $Species, $b3, $Genus, $b4, $Family, $b5, $Order, $b6, $Class, $b7, $Phylum, $b8, $Kingdom, $b9, $Superkingdom, $b10) = split(/\t/, $line, 20);
  if (($Superkingdom eq "Archaea") || ($Superkingdom eq "Bacteria") || ($Superkingdom eq "Viruses")) {
    print OUT "$line\n";
  } elsif ($Superkingdom eq "Eukaryota") {
    print OUT_2 "$line\n";
  }
}

close (IN);
close(OUT);
close(OUT_2);

print "Patientia comes est sapientiae. Patience is the companion of wisdom. - St. Augustine\n Your outfiles are $out_non and $out_euk\n";
