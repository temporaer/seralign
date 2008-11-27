package TriangleWrap;
use strict;
use warnings;

sub getTriangles{
	my $vertices = shift;
	my $fn = "tmp_support";
	open NODE, ">$fn.node" or die $!;
	print NODE join(' ', scalar(keys %$vertices), 2, 0, 0), "\n";
	my $i = 0;
	my %nr2point;
	while (my ($k,$v) = each %$vertices){
		$nr2point{$i} = $k;
		print NODE join (' ', $i++,$$v[0],$$v[1]), "\n";
	}
	close NODE;
	system("../../third_party/triangle/triangle","-Q","$fn.node");
	open ELE, "<$fn.1.ele" or die $!;
	$_=<ELE>;  # ignore 1st line
	my @triangles;
	while(<ELE>){
		# <leading space><line-number> <pointid1> <pointid2> <pointid3>
		next if /^\s*#/;   # ignore comments
		chomp;             # ignore space at end...
		s/^\s+//g;         # ... and beginning
		@_ = split /\s+/;  # get tokens
 		shift;             # ignore 1st token
		push @triangles, [@nr2point{@_}];  # save triangle
	}
	close ELE;
	return @triangles;
}

1;
