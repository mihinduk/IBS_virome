#!/usr/bin/env perl
# cenote-taker2_cat_taxonomy_merger.pl
# 2020_07_08
use warnings;
use strict;

# This script takes the output from CAT and parses it to generate an outfile modeled on parsed Cenote-taker2 output

# $ARGV[0] = parsed DNA output from Cenote-taker2
# $ARGV[1] = parsed RNA output from Cenote-taker2
# $ARGV[2] = parsed output from CAT
# $ARGV[3] = base name for outfiles

my $base = $ARGV[3];
my $comp = "$base". "_CT2_CAT_comp.txt";
my $final_taxonomy = "$base". "_CT2_CAT_contig_taxonomy.txt";

my %CT2 = ();
my %CT_taxonomy = (); # contig linked to ; separated taxonomy
my %CT_taxid = (); # contig linked to taxid
my %CT_out = (); # contig linked to ; separated entire line
my ($taxonomy, $final_tax);

my $CAT_count = 0;
my $ct_count = 0;
my $rule_1 = 0;
my $rule_2 = 0;
my $rule_3 = 0;
my $rule_4 = 0;
my $rule_5 = 0;
my $rule_6 = 0;
my $rule_7 = 0;
my $rule_8 = 0;

# Merge Cenote-taker2 DNA and RNA taxonomy
open(IN, $ARGV[0]) or die "Couldn't find the cenote-taker2 DNA virus file $!";
while (my $line = <IN>) {
  next if ($line =~/^original contig name/);
  chomp ($line);
  my ($ocn, $IS, $Completeness, $CTcn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTN, $GI, $ref, $taxid, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species, $protein) = split (/\t/, $line, 21);
  $taxonomy = "$Kingdom".";$Phylum".";$Class".";$Order".";$Family".";$Genus".";$Species"; #print "\$taxonomy = $taxonomy\n";
  $CT_taxonomy{$ocn} = $taxonomy; #print "$CT_taxonomy{$ocn}\n";
  $CT_out{$ocn} = $line;
  $CT_taxid{$ocn} = $taxid;
}

close(IN);

open(IN, $ARGV[1]) or die "Couldn't find the cenote-taker2 RNA virus file $!";
while (my $line = <IN>) {
  next if ($line =~/^original contig name/);
  chomp ($line);
  my ($ocn, $IS, $Completeness, $CTcn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTN, $GI, $ref, $taxid, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species, $protein) = split (/\t/, $line, 21);
  $taxonomy = "$Kingdom".";$Phylum".";$Class".";$Order".";$Family".";$Genus".";$Species";
  $CT_taxonomy{$ocn} = $taxonomy;
  $CT_out{$ocn} = $line;
  $CT_taxid{$ocn} = $taxid;
}

close(IN);

# Read in CAT results and generate 2 outfiles; 1 to track contigs in both and an outfile for R
open (OUT, ">$final_taxonomy") or die "Couldn't create a file for your merged taxonomy $!";
open (OUT_2, ">$comp") or die "Couldn't create a file for your taxonomy comparison $!";
print OUT "original contig name\tIsolation source\tCompleteness\tCenote-taker contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\tCoverage\tRepetitive\n";
print OUT_2 "original contig name\tIsolation source\tCompleteness\tCenote-taker contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\tCoverage\tRepetitive\n";
open(IN, $ARGV[2]) or die "Couldn't find the parsed CAT file $!";
while (my $line = <IN>) {
	chomp ($line);
  next if (($line =~/^original contig name/)  || ($line =~/^\n/));
  my ($ocn, $IS, $Completeness, $CTcn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTN, $GI, $ref, $taxid, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species, $protein, $Coverage, $Repetitive) = split (/\t/, $line, 23);
  if (defined($CT_taxid{$ocn})) {
    $taxonomy = $CT_taxonomy{$ocn};  #print "\$ocn = $ocn\n"; print "CT taxonomy = $taxonomy\n";
    my($ct_Kingdom, $ct_Phylum, $ct_Class, $ct_Order, $ct_Family, $ct_Genus, $ct_Species) = split(/\;/, $taxonomy, 7);  #print "\$ct_Kingdom = $ct_Kingdom\n";
    if (($Kingdom eq "not classified") || ($Element =~/^dark_matter/)) { # If CAT is not classified, but CT2 is, keep the CT2 taxonomy
      $final_tax = $CT_out{$ocn};
      print OUT "$final_tax\t$Coverage\t$Repetitive\n"; $ct_count++; $rule_1++;
      #print OUT_2 "$final_tax\n$line\n"
    } elsif (($ct_Kingdom eq "Viruses") && ($Kingdom eq "Viruses")) {
      if ($ct_Family eq $Family) { # If CT2 and CAT agree at the family level, keep the CT2 taxonomy
        $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
        print OUT "$final_tax\t$Coverage\t$Repetitive\n"; $ct_count++;  $rule_2++;
        #print OUT_2 "$final_tax\n$line\n";
    } elsif ($ct_Family ne $Family) { # If CT2 and CAT do NOT agree at the family level
        #print "\$ocn = $ocn\n"; #print "CT taxonomy = $taxonomy\n";
        $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
        if (($Element =~/^not classified/) || ($Element =~/^dark_matter/)) { # keep the CT2 taxonomy
          $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
          print OUT "$final_tax\t$Coverage\t$Repetitive\n"; $ct_count++;  $rule_3++;
          #print OUT_2 "$final_tax\tFLAG\n$line\n";
        } else { # keep the CAT taxonomy
          print OUT "$line\n"; $CAT_count++; $rule_4++;
          #print OUT_2 "$final_tax\tFLAG\n$line\n";
        }
      }
    } elsif (($ct_Kingdom eq "Viruses") && ($Kingdom ne "Viruses")) {
      #print "\$Kingdom = $Kingdom\n";
      #print OUT_2 "$final_tax\tFLAG\n$line\n";
      if (($Element =~/^not classified/) || ($Element =~/^dark_matter/)) {
        $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
        print OUT "$final_tax\t$Coverage\t$Repetitive\n"; $ct_count++;  $rule_5++;
        #print OUT_2 "$final_tax\tFLAG\n$line\n";
      } else { # Keep CAT taxonomy
      $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
      print OUT "$line\n"; $CAT_count++;  $rule_6++;
      #print OUT_2 "$final_tax\tFLAG\n$line\n";
      }
    } elsif (($ct_Kingdom ne "Viruses")) {
      $final_tax = $CT_out{$ocn}; #print "\$final_tax = $final_tax\n";
      print OUT "$line\n"; $CAT_count++;  $rule_7++;
      print OUT_2 "$final_tax\tFLAG\n$line\n";
    }
  } elsif(!defined($CT_taxid{$ocn})) { # if only CAT has taxonomy, keep it
    print OUT "$line\n"; $CAT_count++;  $rule_8++;
    #print OUT_2 "$line\tFLAG - CAT ONLY\n";
  }
}

close(IN);
close(OUT);
close(OUT_2);

my $total = $ct_count + $CAT_count;
print "All: CT = $ct_count; CAT = $CAT_count; Total = $total\n";
print "Aut viam inveniam aut faciam. I shall either find a way or make one. - Hannibal.\nYour R taxonomy infile is $final_taxonomy\n";
print "If you chose to compare sequences classified by both classifiers, see $comp\n";

print "Rule count: 1: $rule_1\n2: $rule_2\n3: $rule_3\n4: $rule_4\n5: $rule_5\n6: $rule_6\n7: $rule_7\n8: $rule_8\n";
