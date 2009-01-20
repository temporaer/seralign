#!/usr/bin/perl -w
use Data::Dumper;
use Algorithm::SVM;
use Algorithm::SVM::DataSet;
use Statistics::Descriptive;

my $fixedAtomNum = shift @ARGV || die "Usage: $0 fixedAtomNum";

my @c1 = readClass("/tmp/c1.dmp");
my @c0 = readClass("/tmp/c0.dmp");


shift @c1 while(scalar(@c1) > scalar(@c0));
shift @c0 while(scalar(@c0) > scalar(@c1));

foreach (0..9){
	crossval(10,$_);
}

sub crossval{
	my $cross = shift;
	my $iter  = shift;

	my $size1 = scalar(@c1)/$cross;
	my $size0 = scalar(@c0)/$cross;

	my @testrange1  = ($size1*$iter, $size1*($iter+1));
	my @testrange0  = ($size0*$iter, $size0*($iter+1));
	my @trainrange1 = grep{ $a=$_; not grep{$_==$a}@testrange1}(0..$#c1);
	my @trainrange0 = grep{ $a=$_; not grep{$_==$a}@testrange0}(0..$#c0);

	my @trainset1 = @c1[@trainrange1];
	my @trainset0 = @c0[@trainrange0];

	my @testset1 = @c1[$testrange1[0]..$testrange1[1]];
	my @testset0 = @c0[$testrange0[0]..$testrange0[1]];

	$svm1 = new Algorithm::SVM(Type   => 'one-class');
	$svm0 = new Algorithm::SVM(Type   => 'one-class');


	$svm1->train(flatten(@trainset1));
	$svm0->train(flatten(@trainset0));

	$stat = Statistics::Descriptive::Full->new();

	foreach $ds (@testset1){
		my ($sum1, $sum0) = (0,0);
		foreach $atom ( @$ds ){
			$v   = $svm1->predict_value($atom);
			$sum1 += $v if $v>0;
			$sum0 -= $v if $v<0;

			$v   = $svm0->predict_value($atom);
			$sum1 -= $v if $v<0;
			$sum0 += $v if $v>0;
		}
		my $decision = ($sum1>$sum0)?1:0;
		$stat->add_data(($decision==1)?1:0);
		print "1: Predicted Class: ", $decision, " sure:", ($sum1-$sum0),"\n";
	}

	foreach $ds (@testset0){
		my ($sum1, $sum0) = (0,0);
		foreach $atom ( @$ds ){
			$v   = $svm1->predict_value($atom);
			$sum1 += $v if $v>0;
			$sum0 -= $v if $v<0;

			$v   = $svm0->predict_value($atom);
			$sum1 -= $v if $v<0;
			$sum0 += $v if $v>0;
		}
		my $decision = ($sum1>$sum0)?1:0;
		$stat->add_data(($decision==0)?1:0);
		print "0: Predicted Class: ", $decision, " sure:", ($sum1-$sum0), "\n";
	}
	print "Accuracy: ", $stat->mean(), "\n";
}


sub readClass{
	my $file  = shift;
	open F, "<$file" or die $!;
	$_ = <F>;
	$_ = <F>;
	my @res;
	my @molekuele;
	my $cnt=0;
	while(<F>){
		last unless /^\d+\s/;
		chomp;
		my @a = split/\s/;
		shift @a;
		if($a[0]!=0){
			push @res, new  Algorithm::SVM::DataSet(Label => 1, Data => \@a); 
		}
		$cnt++;
		if($cnt%40 == 0){
		  push @molekuele, [@res];
		  @res = ();
		}
	}
	close F;
	return @molekuele;
}

sub flatten { map @$_, @_ }
