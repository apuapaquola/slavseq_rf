#!/usr/bin/env perl
$^W=1;
use ExtUtils::testlib;
use Gapafim;

my $s=shift;
my $t=shift;


my ($sa, $ma, $ta, $sc)=Gapafim::sw($s,$t);

print join("\n", $sa, $ma, $ta, $sc),"\n";
