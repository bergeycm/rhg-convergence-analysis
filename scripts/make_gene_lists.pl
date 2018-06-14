#!/usr/bin/perl

# ========================================================================================
# --- Make gene lists given a list of per SNP stats
# --- This also creates an annotated BED file for the inputted stats file
# ========================================================================================

use strict;
use warnings;

use List::Util qw( max sum );

my $stat_in = shift;
chomp $stat_in;

my ($chr_col, $pos_col, $stat_col, $has_header);

my ($is_fst, $is_ihs, $is_pbs, $is_pi, $is_bayenv) = 0;

# Has user passed already-annotated BED file
my $is_anno_bed = 0;

if ($stat_in =~ /\.anno.*\.bed/) {
    $is_anno_bed = 1;
}

if ($stat_in =~ /\.weir\.fst/) {
    $chr_col  = 0;
    $pos_col  = 1;
    $stat_col = 2;
    $has_header = 1;
    $is_fst = 1;
} elsif ($stat_in =~ /ihs\./) {
    $chr_col  = 0;
    $pos_col  = 1;
    $stat_col = 6;
    $has_header = 0;
    $is_ihs = 1;
} elsif ($stat_in =~ /pbs\..*\.bed/) {
    $chr_col  = 0;
    $pos_col  = 1;
    $stat_col = 3;
    $is_pbs = 1;
    $has_header = 0;
    if (!$is_anno_bed) {
        $stat_in =~ s/\.bed$//;
    }
} elsif ($stat_in =~ /pi_ratio\./) {
    $chr_col  = 0;
    $pos_col  = 1;
    $stat_col = 3;
    $is_pi = 1;
    $has_header = 0;
    if (!$is_anno_bed) {
        $stat_in =~ s/\.bed$//;
    }
} elsif ($stat_in =~ /bayenv_BF_pass\.bed/) {
    $chr_col  = 0;
    $pos_col  = 1;
    $stat_col = 3;
    $is_bayenv = 1;
    $has_header = 0;
    if (!$is_anno_bed) {
        $stat_in =~ s/\.bed$//;
    }
}

# ----------------------------------------------------------------------------------------
# --- Convert to BED format
# ----------------------------------------------------------------------------------------

my $out_bed;

# Skip if input is annotated BED
if (!$is_anno_bed) {

    print STDERR "- Converting to BED format...\n";

    $out_bed = "${stat_in}.bed";

    if ($is_fst || $is_ihs) {

        my $convert_cmd = "awk -v OFS='\\t' '{if (\$" . ($stat_col + 1) . " != \"-nan\") ";
        $convert_cmd .= "print \$" . ($chr_col + 1). ",\$" . ($pos_col + 1);
        $convert_cmd .= ",\$" . ($pos_col + 1) . "+1,\$" . ($stat_col + 1) . "}' ";
        $convert_cmd .= "$stat_in ";
        # Remove header if it's there
        if ($has_header) {
            $convert_cmd .= "| tail -n +2 ";
        }
        # Reduce first column (chr:pos) into just chr in iHS files
        if ($is_ihs) {
            $convert_cmd .= '| sed -e "s/^\([0-9]\+\):[0-9]\+/chr\1/" ';
        }
        $convert_cmd .= "> $out_bed";

        print STDERR "CMD: [$convert_cmd]\n";

        system($convert_cmd);

    }
}

# ----------------------------------------------------------------------------------------
# --- Annotate BED with genes
# ----------------------------------------------------------------------------------------

my $anno_bed;

# Skip if input is annotated BED
if (!$is_anno_bed) {

    print STDERR "- Annotating BED file...\n";

    $anno_bed = $out_bed;
    $anno_bed =~ s/\.bed/\.anno.bed/;

    my $mapbed_cmd = "module load bedtools/2.26.0;";
    $mapbed_cmd .= "mapBed -a $out_bed ";
    $mapbed_cmd .= "-b refGene/refGene.sort.simple.justGenes.gtf ";
    $mapbed_cmd .= "-c 9 -o collapse ";
    $mapbed_cmd .= "> $anno_bed";

    print STDERR "CMD: [$mapbed_cmd]\n";

    system($mapbed_cmd);

} else {
    $anno_bed = $stat_in;
}

# ----------------------------------------------------------------------------------------
# --- Read in annotated stats
# ----------------------------------------------------------------------------------------

print STDERR "- Parsing data for [$stat_in]...\n";

open (STATS, "<$anno_bed")
    or die "ERROR: Could not open annotated BED file. $!\n";

# Loop through annotated BED file, adding data to hash keyed by gene

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

my %snp_info;

while (<STATS>) {

    my @info = split;

    my $val = $info[3];
    my $genes_str = $info[4];

    # Split genes if there are multiple ones and remove redundant genes
    my @genes = uniq(split(/,/, $genes_str));

    foreach my $gene (@genes) {
        if (exists $snp_info{$gene}) {
            push @{ $snp_info{$gene} }, $val;
        } else {
            $snp_info{$gene} = [$val];
        }
    }
}

close STATS;

# ----------------------------------------------------------------------------------------
# --- Output maximum value by gene
# ----------------------------------------------------------------------------------------

print STDERR "- Outputting maximum value by gene...\n";

my $max_out_file = "${stat_in}.max.txt";
open (MAX, ">$max_out_file")
    or die "ERROR: Could not open output file for maximum values. $!\n";

foreach my $gene (keys %snp_info) {

    my $max = max @{ $snp_info{$gene} };

    print MAX "$gene\t$max\n";
}

close MAX;

# ----------------------------------------------------------------------------------------
# --- Output mean by gene
# ----------------------------------------------------------------------------------------

print STDERR "- Outputting mean value by gene...\n";

my $mean_out_file = "${stat_in}.mean.txt";
open (MEAN, ">$mean_out_file")
    or die "ERROR: Could not open output file for mean values. $!\n";

foreach my $gene (keys %snp_info) {

    my $total = sum @{ $snp_info{$gene} };
    my $num_genes = scalar @{ $snp_info{$gene} };

    my $mean = $total / $num_genes;

    print MEAN "$gene\t$mean\n";
}

close MEAN;

# ----------------------------------------------------------------------------------------
# --- Output high count by gene (|iHS| >= 2)
# ----------------------------------------------------------------------------------------

if ($is_ihs) {

    print STDERR "- Outputting high value count by gene...\n";

    my $high_out_file = "${stat_in}.high.txt";
    open (HIGH, ">$high_out_file")
        or die "ERROR: Could not open output file for high value count. $!\n";

    foreach my $gene (keys %snp_info) {

        my $high_count = 0;

        for my $this_ihs (@{ $snp_info{$gene} }) {
            if ($this_ihs >= 2) {
                $high_count++;
            }
        }
        print HIGH "$gene\t$high_count\n";
    }

    close HIGH;
}

# Clean up a bit
# Skip if input is annotated BED
if (!$is_anno_bed) {
    unlink $out_bed if (!$is_pbs && !$is_bayenv);
}

exit;
