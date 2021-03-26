#!/usr/bin/env perl
# cenote-taker2_parser_v3.pl
# 2020_07_07
use warnings;
use strict;

# This script takes the output from Cenote-taker2 and parses it to generate an outfile with taxonomy for R

# $ARGV[0] = output from Cenote-taker2
# $ARGV[1] = Cenote-taker2 mode: DNA (default) or RNA

# Make $ARGV[1] lc
my $na = lc($ARGV[1]); #print "\$na = $na\n";

my $outfile = $ARGV[0];
$outfile =~ s/\.tsv/_clean_tax.txt/; #print "\$outfile = $outfile\n";

my %non_euk_taxid=();
my %non_euk = ();
my ($unneeded, $to_be_split, $GI, $ref, $to_be_split_again, $Species_to_find, $ne_key, $ne_val, $protein);
my $x = 0;
my @non_euk_match =();
my @merger = ();

# Read the taxonomy files into a hash
# non_euk
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_non_euk_fin_fib.txt") or die "Couldn't find the non_euk file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $taxonomy = "$tax_name;$Kingdom;$Phylum;$Class;$Order;$Family;$Genus;$Species";
	$non_euk {$tax_id} = $taxonomy; #print "\$tax_name $tax_name linked to $key\n";
	$non_euk_taxid{$tax_id} = $tax_name;
}
close(IN);

# Supplemental
open(IN, "/Users/handley_lab/Handley\ Lab\ Dropbox/virome/resources/viral_taxonomy/rankedlineage_non_euk_supp.txt") or die "Couldn't find the non_euk file $!";
while (my $line = <IN>) {
	chomp ($line);
	my ($tax_id, $tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/\t/, $line, 9);
	my $key = "$tax_id;$Kingdom;$Phylum;$Class;$Order;$Family;$Genus;$Species";
	$non_euk {$key} = $tax_name; #print "\$tax_name $tax_name linked to $key\n";
	$non_euk_taxid{$tax_id} = $tax_name;
}
close(IN);

# Parse the DNA file
if ($na eq "dna") {
	open(IN, $ARGV[0]) or die "Couldn't find the cenote-taker2 file $!";
	open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
	print OUT "original contig name\tIsolation source\tCompleteness\tCenote-taker contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\n";
	while (my $line = <IN>) {
		next if ($line =~/^Isolation source/);
		next if ($line =~/^\t/);
	  chomp ($line);
		my ($IS, $Completeness, $CTcn, $ocn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTP, $BLASTN) = split(/\t/, $line, 11);
		$IS=~ s/^\s+|\s+$//g; $Completeness=~ s/^\s+|\s+$//g; $CTcn=~ s/^\s+|\s+$//g; $ocn=~ s/^\s+|\s+$//g; $Length=~ s/^\s+|\s+$//g;
		$Element=~ s/^\s+|\s+$//g; $Topology=~ s/^\s+|\s+$//g; $CVD=~ s/^\s+|\s+$//g; $ORFc=~ s/^\s+|\s+$//g; $BLASTP=~ s/^\s+|\s+$//g;
		$BLASTN=~ s/^\s+|\s+$//g;
		print OUT "$ocn\t$IS\t$Completeness\t$CTcn\t$Length\t$Element\t$Topology\t$CVD\t$ORFc\t$BLASTN\t";
		if ($BLASTP !~ /;/) {
			#print "Here's something odd\n"; print "\$BLASTP = $BLASTP\n";
			$x =0;
			my ($taxa, $extra) = split (/ sp\./, $Element, 2); #print "\$taxa = $taxa\t\$extra = $extra\n";
			@non_euk_match =();
			foreach my $name (keys %non_euk_taxid){
    		if ($non_euk_taxid{$name} eq "$taxa") {
					push (@non_euk_match, $name); #print "@non_euk_match\n"; print "\$name = $name\n";
				}
			}
				my $odd_length = scalar(@non_euk_match); #print "\$odd_length = $odd_length\n";
				if ($odd_length == 1) {
					my $current_taxid = shift(@non_euk_match);
					my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
					my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
					print OUT "GI_undefined\tref_undefined\t$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\tno protein\n";
				} elsif ($odd_length > 1) {
					while ($odd_length > $x) { #print "Found > 1 tax_id for this taxa\n";
						my $current_taxid = shift(@non_euk_match);
						my $values = $non_euk{$current_taxid};
						my ($new_tax_name, $new_Kingdom, $new_Phylum, $new_Class, $new_Order, $new_Family, $new_Genus, $new_Species) = split(/;/, $values, 8);
						if ((scalar @merger == 0)) {
							push(@merger, $new_tax_name, $new_Kingdom, $new_Phylum, $new_Class, $new_Order, $new_Family, $new_Genus, $new_Species);
							$x++;
						} else {
							my $current_tax_name = $merger[0]; my $current_Kingdom = $merger[1]; my $current_Phylum = $merger[2]; my $current_Class= $merger[3];
							my $current_Order = $merger[4]; my $current_Family = $merger[5]; my $current_Genus = $merger[6]; my $current_Species = $merger[7];
							if ($current_tax_name eq $new_tax_name) {$merger[0] = $current_tax_name;} else {$merger[0] = $current_tax_name."_multi";}
							if ($current_Kingdom eq $new_Kingdom) {$merger[1] = $current_Kingdom;} else {$merger[1] = "Kingdom_undefined";}
							if ($current_Phylum eq $new_Phylum) {$merger[2] = $current_Phylum;} else {
								my $m1 =$merger[1]; $merger[2] = "$m1"."_Phylum_undefined";}
							if ($current_Class eq $new_Class) {$merger[3] = $current_Class} else {
								my $m2 = $merger[2]; $merger[3] = "$m2"."_Class_undefined";}
							if ($current_Order eq $new_Order) {$merger[4] = $current_Class} else {
								my $m3 = $merger[3]; $merger[4] = "$m3"."_Order_undefined";}
							if ($current_Family eq $new_Family) {$merger[5] = $current_Family} else {
								my $m4 = $merger[4]; $merger[5] = "$m4"."_Family_undefined";}
							if ($current_Genus eq $new_Genus) {$merger[6] = $current_Genus} else {
								my $m5 = $merger[5]; $merger[6] = "$m5"."_Genus_undefined";}
							if ($current_Species eq $new_Species) {$merger[7] = $current_Species} else {
								my $m6 = $merger[6]; $merger[7] = "$m6"."_Species_undefined_multi";}
							$x++; #print "\$x = $x\n";
						}
					}
					my $current_tax_name = $merger[0]; my $current_Kingdom = $merger[1]; my $current_Phylum = $merger[2]; my $current_Class= $merger[3];
					my $current_Order = $merger[4]; my $current_Family = $merger[5]; my $current_Genus = $merger[6]; my $current_Species = $merger[7];
					print OUT "GI_undefined\tref_undefined\ttaxid_undefined\t$current_Kingdom\t$current_Phylum\t$current_Class\t$current_Order\t$current_Family\t$current_Genus\t$current_Species\tno protein\n";
				} elsif ($odd_length == 0) {
					#print "This should fix $taxa\n";
						my $family_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$taxa";
						my $genus_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$taxa"."_Genus_undefined";
						print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\t$family_out\t$genus_out\t$Element\tno protein\n";
					}
				} else {
				#print "\$BLASTP = $BLASTP\n";
				($unneeded, $to_be_split) = split(/gi\|/, $BLASTP, 2); #print "\$to_be_split = $to_be_split\n";
				($GI, $unneeded, $ref, $to_be_split_again) = split(/\|/, $to_be_split, 4); #print "\$to_be_split_again = $to_be_split_again\t\$GI = $GI\t\$ref = $ref\n";
				print OUT "$GI\t$ref\t";
				($protein, $to_be_split) = split(/\[/, $to_be_split_again, 2); #print "\$protein = $protein\t\$to_be_split = $to_be_split\n";
				$protein=~ s/^\s+|\s+$//g;
				($Species_to_find, $unneeded) = split(/\]/, $to_be_split, 2);  #print "\$Species_to_find = $Species_to_find\n";
				@non_euk_match =();
				foreach my $name (keys %non_euk_taxid){
	    		if ($non_euk_taxid{$name} eq "$Species_to_find") {
						push (@non_euk_match, $name); #print @non_euk_match; print "\n"; #print "\$Species_to_find = $Species_to_find\n";
					}
				}
			my $length = scalar(@non_euk_match); #print "\$length = $length\n";
			if ($length == 1) { #print "Full taxa 1\n";
				my $current_taxid = shift(@non_euk_match);
				my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
				my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
				print OUT "$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n"; #print "OUT: $current_taxid\t$Kingdom\n";
			} elsif ($length >1) { #print "Full taxa >1\n";
				my $current_taxid = shift(@non_euk_match);
				my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
				my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
				print OUT "$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n";
			} else{
				print "Nothing matches $Species_to_find\n";
						print OUT "taxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\t";
						print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\t";
						print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Species_to_find\t$protein\n";
					}
				}
			}
} elsif ($na eq "rna") {
	open(IN, $ARGV[0]) or die "Couldn't find the cenote-taker2 file $!";
	open (OUT, ">$outfile") or die "Couldn't create a file for your clean taxonomy $!";
	print OUT "original contig name\tIsolation source\tCompleteness\tCenote-taker contig name\tLength\tElement Name\tTopology\tCommon Viral Domains\tORF caller used\tBLASTN result (if any)\tGI\tref\ttaxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tProtein\n";
	while (my $line = <IN>) {
		next if ($line =~/^Isolation source/);
		next if ($line =~/^\t/);
	  chomp ($line);
		my ($IS, $Completeness, $CTcn, $ocn, $Length, $Element, $Topology, $CVD, $ORFc, $BLASTP, $BLASTN) = split(/\t/, $line, 11);
		$IS=~ s/^\s+|\s+$//g; $Completeness=~ s/^\s+|\s+$//g; $CTcn=~ s/^\s+|\s+$//g; $ocn=~ s/^\s+|\s+$//g; $Length=~ s/^\s+|\s+$//g;
		$Element=~ s/^\s+|\s+$//g; $Topology=~ s/^\s+|\s+$//g; $CVD=~ s/^\s+|\s+$//g; $ORFc=~ s/^\s+|\s+$//g; $BLASTP=~ s/^\s+|\s+$//g;
		$BLASTN=~ s/^\s+|\s+$//g;
		my ($not_wanted, $true_contig) = split(/\@/, $ocn, 2); #print "\$true_contig = $true_contig\n";
		print OUT "$true_contig\t$IS\t$Completeness\t$CTcn\t$Length\t$Element\t$Topology\t$CVD\t$ORFc\t$BLASTN\t";
		if ($BLASTP !~ /;/) {
			#print "Here's something odd\n"; print "\$BLASTP = $BLASTP\n";
			$x =0;
			my ($taxa, $extra) = split (/ sp\./, $Element, 2); #print "\$taxa = $taxa\t\$extra = $extra\n";
			@non_euk_match =();
			foreach my $name (keys %non_euk_taxid){
    		if ($non_euk_taxid{$name} eq "$taxa") {
					push (@non_euk_match, $name); #print "@non_euk_match\n"; print "\$name = $name\n";
				}
			}
				my $odd_length = scalar(@non_euk_match); #print "\$odd_length = $odd_length\n";
				if ($odd_length == 1) {
					my $current_taxid = shift(@non_euk_match);
					my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
					my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
					print OUT "GI_undefined\tref_undefined\t$current_taxid\tKingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\tno protein\n";
				} elsif ($odd_length > 1) {
					while ($odd_length > $x) { #print "Found > 1 tax_id for this taxa\n";
						my $current_taxid = shift(@non_euk_match);
						my $values = $non_euk{$current_taxid};
						my ($new_tax_name, $new_Kingdom, $new_Phylum, $new_Class, $new_Order, $new_Family, $new_Genus, $new_Species) = split(/;/, $values, 8);
						if ((scalar @merger == 0)) {
							push(@merger, $new_tax_name, $new_Kingdom, $new_Phylum, $new_Class, $new_Order, $new_Family, $new_Genus, $new_Species);
							$x++;
						} else {
							my $current_tax_name = $merger[0]; my $current_Kingdom = $merger[1]; my $current_Phylum = $merger[2]; my $current_Class= $merger[3];
							my $current_Order = $merger[4]; my $current_Family = $merger[5]; my $current_Genus = $merger[6]; my $current_Species = $merger[7];
							if ($current_tax_name eq $new_tax_name) {$merger[0] = $current_tax_name;} else {$merger[0] = $current_tax_name."_multi";}
							if ($current_Kingdom eq $new_Kingdom) {$merger[1] = $current_Kingdom;} else {$merger[1] = "Kingdom_undefined";}
							if ($current_Phylum eq $new_Phylum) {$merger[2] = $current_Phylum;} else {
								my $m1 =$merger[1]; $merger[2] = "$m1"."_Phylum_undefined";}
							if ($current_Class eq $new_Class) {$merger[3] = $current_Class} else {
								my $m2 = $merger[2]; $merger[3] = "$m2"."_Class_undefined";}
							if ($current_Order eq $new_Order) {$merger[4] = $current_Class} else {
								my $m3 = $merger[3]; $merger[4] = "$m3"."_Order_undefined";}
							if ($current_Family eq $new_Family) {$merger[5] = $current_Family} else {
								my $m4 = $merger[4]; $merger[5] = "$m4"."_Family_undefined";}
							if ($current_Genus eq $new_Genus) {$merger[6] = $current_Genus} else {
								my $m5 = $merger[5]; $merger[6] = "$m5"."_Genus_undefined";}
							if ($current_Species eq $new_Species) {$merger[7] = $current_Species} else {
								my $m6 = $merger[6]; $merger[7] = "$m6"."_Species_undefined_multi";}
							$x++; #print "\$x = $x\n";
						}
						my $current_tax_name = $merger[0]; my $current_Kingdom = $merger[1]; my $current_Phylum = $merger[2]; my $current_Class= $merger[3];
						my $current_Order = $merger[4]; my $current_Family = $merger[5]; my $current_Genus = $merger[6]; my $current_Species = $merger[7];
						print OUT "GI_undefined\tref_undefined\ttaxid_undefined\t$current_Kingdom\t$current_Phylum\t$current_Class\t$current_Order\t$current_Family\t$current_Genus\t$current_Species\tno protein\n";
						}
					} elsif ($odd_length == 0) {
						#print "This should fix $taxa\n";
						my $family_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$taxa";
						my $genus_out = "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined"."_$taxa"."_Genus_undefined";
						print OUT "GI_undefined\tref_undefined\ttaxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\t$family_out\t$genus_out\t$Element\tno protein\n";
					}
			} else {
				#print "\$BLASTP = $BLASTP\n";
				($unneeded, $to_be_split) = split(/gi\|/, $BLASTP, 2); #print "\$to_be_split = $to_be_split\n";
				($GI, $unneeded, $ref, $to_be_split_again) = split(/\|/, $to_be_split, 4); #print "\$to_be_split_again = $to_be_split_again\t\$GI = $GI\t\$ref = $ref\n";
				print OUT "$GI\t$ref\t";
				($protein, $to_be_split) = split(/\[/, $to_be_split_again, 2); #print "\$protein = $protein\t\$to_be_split = $to_be_split\n";
				$protein=~ s/^\s+|\s+$//g;
				($Species_to_find, $unneeded) = split(/\]/, $to_be_split, 2);  #print "\$Species_to_find = $Species_to_find\n";
				@non_euk_match =();
				foreach my $name (keys %non_euk_taxid){
	    		if ($non_euk_taxid{$name} eq "$Species_to_find") {
						push (@non_euk_match, $name); #print @non_euk_match; print "\n"; #print "\$Species_to_find = $Species_to_find\n";
					}
			}
			my $length = scalar(@non_euk_match); #print "\$length = $length\n";
			if ($length == 1) { #print "Full taxa 1\n";
				my $current_taxid = shift(@non_euk_match);
				my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
				my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
				print OUT "$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n"; #print "OUT: $current_taxid\t$Kingdom\n";
			} elsif ($length >1) { #print "Full taxa >1\n";
				my $current_taxid = shift(@non_euk_match);
				my $values = $non_euk{$current_taxid}; #print "\$values\t$values\n";
				my ($tax_name, $Kingdom, $Phylum, $Class, $Order, $Family, $Genus, $Species) = split(/;/, $values, 8);
				print OUT "$current_taxid\t$Kingdom\t$Phylum\t$Class\t$Order\t$Family\t$Genus\t$Species\t$protein\n";
			} else{
				print "Nothing matches $Species_to_find\n";
						print OUT "taxid_undefined\tKingdom_undefined\tKingdom_undefined_Phylum_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined\t";
						print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined\tKingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined\t";
						print OUT "Kingdom_undefined_Phylum_undefined_Class_undefined_Order_undefined_Family_undefined_Genus_undefined\t$Species_to_find\t$protein\n";
					}
				}
			}
		}

close (IN);
close(OUT);

print "Ingenio maximus, arte rudis. Maximum ingenuity, raw technique. - Ovid\nYour outfile is $outfile\n";
