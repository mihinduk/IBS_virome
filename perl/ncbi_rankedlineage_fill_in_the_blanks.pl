#!/usr/bin/env perl
# ncbi_rankedlineage_fill_in_the_blanks.pl
# 2020_07_01
use warnings;
use strict;

# This script fills in blanks in the taxonomy files derived from rankedlineage.dmp

# $ARGV[0] = file derived from NCBI rankedlineage.dmp to be parsed

my $outfile = $ARGV[0];
$outfile =~ s/\.txt/_fib.txt/;

# Pre-script
# grep "Viruses" rankedlineage.dmp > rankedlineage_viruses.txt
# echo -e "tax_id\ttax_name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > rankedlineage_viruses_fin.txt
# more rankedlineage_viruses.txt | cut -f1,3,19,15,13,11,9,7,5 >> rankedlineage_viruses_fin.txt


# Read in taxonomy and replace blanks with the previous category _unknown
open(IN, $ARGV[0]) or die "Couldn't find the taxonomy file $!";
open (OUT, ">$outfile") or die "Couldn't create a file for your blankless taxonomy $!";
while (my $line = <IN>) {
  chomp ($line);
  if ($line =~/^tax_id/) {
    print OUT "$line\n"; } else {
      my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
      if ($tax_name eq "NA") {$tax_name = "tax-name_undefined"};
      if ($Kingdom eq "NA") {$Kingdom = "Kingdom_undefined"};
      if ($Phylum eq "NA") {$Phylum = "$Kingdom"."_undefined_Phylum"};
      if ($Class eq "NA") {$Class = "$Phylum"."_undefined_Class"};
      if ($Order eq "NA") {$Order = "$Class"."_undefined_Order"};
      if ($Family eq "NA") {$Family = "$Order"."_undefined_Family"};
      if ($Genus eq "NA") {$Genus = "$Family"."_undefined_Genus"};
      if ($Species eq "NA") {$Species = "$Genus"."_undefined_Species"};
      print OUT "$tax_id\t$tax_name\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\n";
    }
}

close(IN);
close(OUT);

print "Ingenio maximus, arte rudis. Maximum ingenuity, raw technique. - Ovid\nYour outfile is $outfile\n";
