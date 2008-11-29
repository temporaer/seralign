use Lingua::POSAlign;
use Statistics::Descriptive;
my $align = new Lingua::POSAlign;
%Lingua::POSAlign::penalty = qw(X -1 _ 0);

$fn1 = shift @ARGV || die "need 2 files";
$fn2 = shift @ARGV || die "need 2 files";

@ar1 = getAsArr($fn1);
@ar2 = getAsArr($fn2);

$stats_AB = Statistics::Descriptive::Full->new();
$stats_AA = Statistics::Descriptive::Full->new();
$stats_BB = Statistics::Descriptive::Full->new();

foreach $a (@ar1){
	foreach $b(@ar2){
		$align->align($a,$b);
		$stats_AB->add_data($align->total_score);
	}
}
foreach $a (@ar1){
	foreach $b(@ar1){
		next if($a==$b);
		$align->align($a,$b);
		$stats_AA->add_data($align->total_score);
	}
}
foreach $a (@ar2){
	foreach $b(@ar2){
		next if($a==$b);
		$align->align($a,$b);
		$stats_BB->add_data($align->total_score);
	}
}


$mAB = $stats_AB->mean();
$mAA = $stats_AA->mean();
$mBB = $stats_BB->mean();

print "AB=$mAB AA=$mAA BB=$mBB\n";


sub getAsArr{
	open(IN, "<".shift) or die $!;
	my @arr = ();
	while(<IN>){
		my ($l) = /^\w+\s*:\s*(\w+)\s*$/;
		push @arr, [split //,$l];
	}
	return @arr;
}
