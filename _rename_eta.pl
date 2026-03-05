#!/usr/bin/env perl
# Bulk rename Unicode η → nrna_frac / nRNA fraction
use strict;
use warnings;

# Read all files from ARGV
while (<ARGV>) {
    # catch-all: replace any remaining η with nrna_frac
    s/η/nrna_frac/g;
    print;
}
