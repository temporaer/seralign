#!/usr/bin/perl -w

foreach $f( <s*.pl> ){
	open F, "<$f" or die $!;
	@lines = <F>;
	@pos   = @neg   = ();
	push @pos, map { /\((.*)\)/; $1 } grep {  /^active/} @lines;
	push @neg, map { /\((.*)\)/; $1 } grep { !/^active/} @lines;
	close F;

	$f =~ /^(.*)\.pl/;
	open OUT, ">$1.py" or die $!;
	print OUT "pos = {}\n";
	map {print OUT "pos[$_] = True\n"}@pos;
	map {print OUT "pos[$_] = False\n"}@neg;
	close OUT;
}
