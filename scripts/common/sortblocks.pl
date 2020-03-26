#!/usr/bin/perl
$|=1;				# auto-flush

# sort a haploblocks file by chromosome

use strict;
no strict "refs";

my $haploblocksFileFull = $ARGV[0];
my $sortedBlocksFile = $ARGV[1];

# comparison function for chromosome names.
sub chr_cmp {
    my ($chra, $chrb) = @_;
    # numerics < non-numerics
   my $c1 = ( ($chra =~ /\D+/) <=> ($chrb =~ /\D+/) );
   if ($c1 != 0) {
       return $c1;
   }
   else {
       return ($chra <=> $chrb);
   }
}

# compare chromosomes and blocks
sub chr_blk_cmp {
    my $c1 = chr_cmp($a->[0], $b->[0]);
    if ($c1 != 0) {
	return $c1;
    }
    else {
	return ($a->[1] <=> $b->[1]);
    }
}

open(HAPLOBLOCKS, "<$haploblocksFileFull") || die("<pre>Internal error: results file $haploblocksFileFull open failed: $! </pre>\n");
my @allblocks = ();
while (<HAPLOBLOCKS>) {
    chomp($_);
    my @blockline = split(/\t/, $_);
    push @allblocks, \@blockline;
}
close(HAPLOBLOCKS);

my @sblocks = sort { chr_blk_cmp($a, $b) } @allblocks;

open(SORTEDBLOCKS, ">$sortedBlocksFile") || die("<pre>Internal error: results file $sortedBlocksFile open failed: $! </pre>\n");

foreach my $blockref (@sblocks) {
    print SORTEDBLOCKS join("\t", @{$blockref}), "\n";
}
close(SORTEDBLOCKS);


