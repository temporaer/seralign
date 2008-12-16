#!/usr/bin/perl -w

@pos   = @neg   = ();
foreach $f( <s*.pl> ){
	open F, "<$f" or die $!;
	@lines = <F>;
	push @pos, map { /\((.*)\)/; $1 } grep {  /^active/} @lines;
	push @neg, map { /\((.*)\)/; $1 } grep { !/^active/} @lines;
	close F;
}
open OUT, ">is_active.txt" or die $!;
map {print OUT qq|$_ 1\n|}@pos;
map {print OUT qq|$_ 0\n|}@neg;
close OUT;
