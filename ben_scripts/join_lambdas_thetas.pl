#!/usr/bin/perl
use strict;
use warnings;

# Check for command line arguments
if (@ARGV != 2) {
    die "Usage: $0 <lambda_file> <theta_file>\n";
}

my ($lambda_file, $theta_file) = @ARGV;

# Read the lambda file into a hash based on the first 3 columns
my %lambda_data;
open my $lambda, '<', $lambda_file or die "Cannot open $lambda_file: $!";
LINE: while (my $line=<$lambda>) {
    chomp $line;
    my @fields = split(/\t/,$line);
    # Lambda has a header but theta doesn't
    if ($fields[0] =~ /^#\s*chrom/) {
        $fields[0] = "chrom";
        unshift(@fields, qw/chr start end theta_mean_h theta_mean_l silhouette best_cell_type 2nd_best_cell_type best_lambda 2nd_best_lambda/);
        print join("\t", @fields), "\n";
    }
    else
    {
        # Cell type fields are $field[6..$#fields]
        my @cell_type_flds = @fields[6..$#fields];
        my ($min_index, $second_min_index, $min_val, $second_min_val, $sil) = silhouette(@cell_type_flds);

        my $key = join("\t", @fields[0..2]);
        $lambda_data{$key} = join("\t", $sil, $min_index, $second_min_index, $min_val, $second_min_val, @fields);
    }
}
close $lambda;

# Read the theta file and compare with the lambda hash
open my $theta, '<', $theta_file or die "Cannot open $theta_file: $!";
while (<$theta>) {
    chomp;
    my @fields = split /\t/;

    my $key = join("\t", @fields[0..2]);
    if (exists $lambda_data{$key}) {
        # Remove brackets and split columns 7 and 8
        my $column7_mean = mean_list($fields[6]);
        my $column8_mean = mean_list($fields[7]);
        
        # Replace columns 7 and 8 with their means
        $fields[6] = $column7_mean;
        $fields[7] = $column8_mean;

        my $new_line = join("\t", @fields[3..7]);
        print "$new_line\t$lambda_data{$key}\n";
    }
    else
    {
        die "Key $key not found in lambda hash\n";
    }
}
close $theta;

sub silhouette
{
    my (@vals) = @_;

    my $min_val = 1E5;
    my $min_index = 1;

    # Find the minimum value and its index
    for my $i (0..$#vals) {
        if ($vals[$i] < $min_val) {
            $min_val = $vals[$i];
            $min_index = $i + 1;
        }
    }

    my $second_min_val = 1E5;
    my $second_min_index = 1;

    # Find the second lowest value and its index, excluding the minimum value
    for my $i (0..$#vals) {
        if ($vals[$i] < $second_min_val && $vals[$i] != $min_val) {
            $second_min_val = $vals[$i];
            $second_min_index = $i + 1;
        }
    }

    my $sil = $second_min_val - $min_val;

    return ($min_index, $second_min_index, $min_val, $second_min_val, $sil);
}


sub mean_list {
    my ($str) = @_;

    $str =~ s/^\s*\[\s*//g;
    $str =~ s/\s*\]\s*$//g;
    my @values = split(/,/, $str);
    return mean(@values);
}

# Function to calculate mean
sub mean {
    my @values = @_;
    my $sum = 0;
    foreach my $value (@values) {
        $sum += $value;
    }
    return $sum / @values;
}
