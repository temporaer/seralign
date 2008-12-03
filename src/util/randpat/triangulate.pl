#!/usr/bin/perl -w

# This script triangulates a given number of points in the plane
# It uses TriangleWrap, which in turn utilizes Triangle
# available at http://www.cs.cmu.edu/~quake/triangle.html

# Input:
# A list of nodes of the form:
#   NodeName X-Coordinate Y-Coordinate
#   NodeName X-Coordinate Y-Coordinate
#   NodeName X-Coordinate Y-Coordinate

# Output:
# A (redundant!) list of the form:
#   NodeName1 NodeName2
#   NodeName1 NodeName2
#   NodeName1 NodeName2
# which describes the nodes which are connected by the triangulation.

# Author: Hannes Schulz <mail at hannes-schulz dot de>

BEGIN{push @INC, "@SERALIGN_SCRIPT_DIR@"}

use TriangleWrap;
use Data::Dumper;

my %p;
while(<>){
	chomp;
	s/^\s*//g;
	next if /^$/;
	next if /^#/;
	my($name,$x,$y) = split /\s+/,$_;
	$p{$name} = [$x,$y];
}
@ts = TriangleWrap::getTriangles(\%p);
foreach $t(@ts){
	print "$$t[0] $$t[1]\n";
	print "$$t[1] $$t[2]\n";
	print "$$t[2] $$t[0]\n";
}

__DATA__
A 0 0
B 0 1
C 1 0
D 1 1
E 0.8 0.6
