#!/usr/bin/perl -w

use Smart::Comments;
use Data::Dumper;

my @curModel;
my $curModelName = "";
while(<>){
	chomp;
	next if /^\s*$/;
	if(/^begin\((.*)\)\.$/){
		$curModelName = $1;
		next;
	}
	if(/^end\((.*)\)\.$/){
		$mn = $1;
		unless($mn eq $curModelName){
			die "Syntax Error: Sec Start `$curModelName' != Sec End `$mn'\n";
		}
		readModule($curModelName, @curModel);
		@curModel=();
		next;
	}
	push @curModel, $_;
}
print "EOF\n";

sub readModule{
	my $name  = shift();
	my @lines = @_;
	my @objs = grep{$_}map{/\b(o\d+)/ and $1}@lines;
	my %tmp;
	map{$tmp{$_}++}@objs;
	@objs = keys %tmp;
	@objs = sort{ $a cmp $b } @objs;
	#print "$name: @objs\n";
	%properties = ();
	%relations  = ();
	foreach(@lines){
		if(/^(\w+)\((\w+)\)\.$/){ # one-arg
			$properties{$2} = [] unless exists $properties{$2};
			push @{$properties{$2}}, $1;
		}
		if(/^(\w+)\((\w+),(\w+)\)\.$/){ # two-arg
			if($1 eq "config"){
				push @{$properties{$2}}, $3;
				next;
			}
			#$relations{$1} = {} unless exists $properties{$1};
			$relations{$1}{"$2-$3"}++;
		}
	}
	print scalar @objs, "\n";
	foreach my $i (0..$#objs){
		next unless $properties{$objs[$i]};
		print "$i ", join(" ", @{$properties{$objs[$i]}}), "\n";
	}
	foreach my $r (keys %relations){
		my @pairs = ();
		foreach my $o1 (0..$#objs){
			foreach my $o2 (0..$#objs){
				next unless exists ${$relations{$r}}{$objs[$o1] . "-" . $objs[$o2]};
				push @pairs, "$o1 $o2";
			}
		}
		print "$r ", join(" ", @pairs), "\n";
	}
	print "END\n";
}
