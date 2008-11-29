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

$fn1 = shift @ARGV || die "need 2 files";
$fn2 = shift @ARGV || die "need 2 files";

@ar1 = getAsArr($fn1);
@ar2 = getAsArr($fn2);

$stats_AB = Statistics::Descriptive::Full->new();
$stats_AA = Statistics::Descriptive::Full->new();
$stats_BB = Statistics::Descriptive::Full->new();

$x = [ qw(X X X X X X X) ];

foreach $a (@ar1){
	$score = $align->align($a,$x);
	$stats_AA->add_data($score);
}
foreach $a (@ar2){
	$score = $align->align($a,$x);
	$stats_BB->add_data($score);
}


$mAA = $stats_AA->mean();
$mBB = $stats_BB->mean();

#print "AB=$mAB AA=$mAA BB=$mBB\n";
print "AX=$mAA BX=$mBB\n";


sub getAsArr{
	open(IN, "<".shift) or die $!;
	my @arr = ();
	while(<IN>){
		my ($l) = /^\w+\s*:\s*(\w+)\s*$/;
		push @arr, [split //,$l];
	}
	return @arr;
}
