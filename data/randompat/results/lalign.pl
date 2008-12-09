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

@sim_mat = ();
foreach $a (@ar1, @ar2){
	@sim_row = ();
	foreach $b(@ar1, @ar2){
		$align->align($a,$b);
		$stats_BB->add_data($align->total_score);
		push @sim_row, $align->total_score;
	}
	push @sim_mat, [@sim_row];
}

print join(',', map{"w$_"} (0..(scalar(@ar1)+scalar(@ar2)))), "\n";
foreach $r (@sim_mat){
	@r = @$r;
	print join(',',@r),"\n";
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
