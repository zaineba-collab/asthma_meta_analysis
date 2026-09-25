#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../tools/ensembl-vep";
use Bio::DB::HTS::Faidx;

@ARGV == 2 or die "Usage: $0 GRCh38.fa corrected_exposure.tsv\n";
my ($fasta, $table) = @ARGV;
my $fa = Bio::DB::HTS::Faidx->new($fasta);
open my $fh, '<', $table or die "Cannot open $table: $!\n";
my $header = <$fh>;
chomp $header;
my @h = split /\t/, $header, -1;
my %i; @i{@h} = (0 .. $#h);
for my $required (qw(SNP chromosome.grch38 position.grch38)) {
  exists $i{$required} or die "Missing column $required in $table\n";
}
print "SNP\tREF_FASTA\n";
while (my $line = <$fh>) {
  chomp $line;
  my @x = split /\t/, $line, -1;
  my ($snp, $chr, $pos) = @x[@i{qw(SNP chromosome.grch38 position.grch38)}];
  next if !defined($pos) || $pos eq '' || $pos eq 'NA';
  my ($base, $length) = $fa->get_sequence("$chr:$pos-$pos");
  $base = defined($base) && length($base) ? uc($base) : 'N';
  print "$snp\t$base\n";
}
