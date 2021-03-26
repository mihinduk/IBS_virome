#!/usr/bin/env perl
# cat_taxonomy_supplement.pl
# 2020_07_08
use warnings;
use strict;

# This script takes the output from CAT and parses it to generate an outfile modeled on parsed Cenote-taker2 output

# $ARGV[0] = output from CAT (contig.taxonomy)
# $ARGV[1] = flye assembler outfile (assembly_info.txt)
# $ARGV[2] = CAT.scores from CAT


my $infile = $ARGV[1];
my $infile_2 = $ARGV[2];

my $outfile = $ARGV[0];
$outfile = "$ARGV[0]"."_clean_tax.txt";

my %ncbi = ();
my %taxid = ();
my %contig_length = ();
my %Topology = ();
my %Coverage = ();
my %Repetitive = ();
my @line = ();

my ($IS, $Completeness, $CTcn, $ocn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTP, $BLASTN, $taxonomy, $tax_name, $Kingdom_now, $Phylum_now, $Class_now, $Order_now, $Family_now, $Genus_now, $Species_now);

my $cat_number = 999;
my $dark_number = 999;

# Read the taxonomy files into a hash
# non_euk
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_non_euk_fin_fib.txt") or die "Couldn't find the non_euk file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $taxonomy = "$tax_name;$Kingdom;$Phylum;$Class;$Order;$Family;$Genus;$Species";
	$ncbi{$tax_id} = $taxonomy; #print "\$tax_name $tax_name linked to $key\n";
}
close(IN);

# euk
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_euk_fin_fib.txt") or die "Couldn't find the euk file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $taxonomy = "$tax_name;$Kingdom;$Phylum;$Class;$Order;$Family;$Genus;$Species";
	$ncbi{$tax_id} = $taxonomy; #print "\$tax_name $tax_name linked to $key\n";
}
close(IN);

# Use the flye assembler assembly_info.txt file to get the length and linearity of all input contigs
#seq_name	length	cov.	circ.	repeat	mult.	alt_group	graph_path
open(IN, "$infile") or die "Couldn't find the flye out file $!";
while (my $line = <IN>) {
	chomp ($line);
  next if ($line =~/^#/);
  my ($seq_name, $length, $cov, $circ, $repeat, $mult, $alt_group, $graph_path) = split(/\t/, $line, 8);
  #contig_5	60791	10	N	N	8	*	*,-57,5,19,19,*
  if ($circ eq "N") {
    $circ = "linear";
  } elsif ($circ eq "Y") {
    $circ = "circular";
  } else {
      $circ = "unknown";
  }
  $contig_length{$seq_name} = $length;
  $Topology{$seq_name} = $circ;
  $Coverage{$seq_name} = $cov;
  $Repetitive{$seq_name} = $repeat;
}
close(IN);

# Get taxids from CAT.scores file
open(IN, "$infile_2") or die "Couldn't find the CAT.scores file $!";
while (my $line = <IN>) {
	chomp ($line);
  next if ($line =~/^#/);
  my ($contig, $classification, $reason, $lineage, $lineage_scores) = split(/\t/, $line, 5);
  if (!defined($lineage)) {
    $lineage = "taxid_undefined";
    push (@line, $lineage);
  } elsif ($lineage =~/\;/) {
    $lineage=~ s/\*$//;
    @line = split(/;/, $lineage);
  } else {
    push (@line, $lineage);
  }
  my $lineback = pop(@line);
  $taxid{$contig} = $lineback;
}

# Read in CAT taxonomy, remove asterisks (if present) and create parsed CT2 columns
# contig	superkingdom	phylum	class	order	family	genus	species
open(IN, "$ARGV[0]") or die "Couldn't find the CAT taxonomy to parse $!";
open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
print OUT "original contig name\tIsolation source\tCompleteness\tCenote-taker contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\tCoverage\tRepetitive\n";
while (my $line = <IN>) {
	chomp ($line);
  next if ($line =~/^#/);
  my ($contig, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 8);
  print OUT "$contig\tunknown\tunknown\tNA\t";
  if (defined($Kingdom)) {
    $Kingdom=~ s/\*//; $Phylum=~ s/\*//; $Class=~ s/\*//; $Order=~ s/\*//; $Family=~ s/\*//; $Genus=~ s/\*//; $Species=~ s/\*//;
    my $length = $contig_length{$contig};
    my $circ = $Topology{$contig};
    print OUT "$length\t";
    my $tax_id = $taxid{$contig}; #print "\$tax_id = $tax_id\n";
    if (!defined($tax_id) || !defined($ncbi{$tax_id})) {$tax_id = "taxid_undefined";} #print "\$tax_id = $tax_id\n";
    if ($tax_id eq "taxid_undefined") {
        undef($taxonomy);
      } else {
        $taxonomy = $ncbi{$tax_id}; #print "\$taxonomy = $taxonomy\n";
      }
    if (defined($taxonomy)) {
      ($tax_name, $Kingdom_now, $Phylum_now, $Class_now, $Order_now, $Family_now, $Genus_now, $Species_now) = split(/;/, $taxonomy, 8);
      if ($Phylum eq "NA") { $Phylum = $Phylum_now; }
      if ($Class eq "NA") { $Class = $Class_now; }
      if ($Order eq "NA") { $Order = $Order_now; }
      if ($Family eq "NA") { $Family = $Family_now; }
      if ($Genus eq "NA") { $Genus = $Genus_now; }
      $cat_number++; $Element = "$tax_name"." sp. cat"."$cat_number";
    } else {
      $cat_number++; $Element = "$Species"." sp. cat"."$cat_number";
    }
    print OUT "$Element\t$circ\tNA\tProdigal (meta)\tno blastn\tGI_undefined\tref_undefined\t";
    my $current_taxid = $taxid{$contig};
    print OUT "$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t";
    my $cov = $Coverage{$contig};
    my $repeat = $Repetitive{$contig};
    print OUT "no protein\t$cov\t$repeat\n";
  } else {
    my $length = $contig_length{$contig}; #print "\$contig = $contig\t\$length = $length\n";
    $dark_number ++; my $dark_element = "dark_matter sp. cat"."$dark_number"; #print "\$dark_element = $dark_element\n";
    my $circ = $Topology{$contig};  #print "\$circ = $circ\n";
    print OUT "$length\t$dark_element\t$circ\tNA\tProdigal (meta)\tno blastn\tGI_undefined\tref_undefined\t";
    my $current_taxid = $taxid{$contig};  #print "\$current_taxid = $current_taxid\n";
    my $cov = $Coverage{$contig};  #print "\$cov = $cov\n";
    my $repeat = $Repetitive{$contig};  #print "\$repeat = $repeat\n";
    my $Species = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_"."$dark_element";
    print OUT "$current_taxid\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\t";
    print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family\t";
    print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus\t$Species\tno protein\t$cov\t$repeat\n";
  }
}

close(IN);
close(OUT);

print "\nQuando omni flunkus moritatus - When all else fails play dead\nYour parsed contig taxonomy file is: $outfile\n";
