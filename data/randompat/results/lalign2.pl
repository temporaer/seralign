use Algorithm::NeedlemanWunsch;
use Statistics::Descriptive;

sub score_sub {
	if (!@_) {
		return -1; # gap penalty
	}

	return ($_[0] eq $_[1]) ? 1 : -1;
}
my $align = Algorithm::NeedlemanWunsch->new(\&score_sub);
$align->local(1);
$align->gap_open_penalty(-6);
$align->gap_extend_penalty(-1);


foreach $f(@ARGV){

	@ar1 = getAsArr($f);

	$stats = Statistics::Descriptive::Full->new();

	$x = [ qw(X X X X X X X) ];

	foreach $a (@ar1){
		$score = $align->align($a,$x);
		$stats->add_data($score);
	}
	print "$f : ", $stats->mean(), "\n";
}




sub getAsArr{
	open(IN, "<".shift) or die $!;
	my @arr = ();
	while(<IN>){
		my ($l) = /^\w+\s*:\s*(\w+)\s*$/;
		push @arr, [split //,$l];
	}
	return @arr;
}
