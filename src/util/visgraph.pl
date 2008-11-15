use GraphViz;
use List::Util qw/reduce max sum/;
use Data::Dumper;

$dir = shift @ARGV;

$adj   = `cat ../../build/$dir/adjmat.dat`;
@paths = split /\n/, `cat ../../build/$dir/paths.dat`;

sub getmat{
	$_ = shift();
	s/\](\s*;\s*)+$//g;
	s/^.*\[/[/g;
	s/(\d)\s+(\d)/\1,\2/g;
	s/(\d)\s+(\d)/\1,\2/g;
	s/(\d)\s*;\s*(\d)/\1],[\2/g;
	s/;/]/;
	return $_;
}

#create your adjacency matrix, nodes do whatever else you need.
my @admat = (
	eval(getmat($adj))
);

@nodes = map{ $_ } (0..$#admat);

#create a new graph
my $graph = GraphViz->new( 
	layout => 'circo',
	directed => 0);

#add all your nodes to the graph
foreach (@nodes) {
	#print "adding node $_\n";
    $graph->add_node($_);
}

my %edgesDone = ();
my %pathstrength = ();
print "Create paths...\n";
@colors = qw(red orange green blue);
foreach $p (0..$#paths){
	addPathStrength(split(/\s+/, $paths[$p]));
	#drawpath($colors[$p],  split(/\s+/, $paths[$p]));
}
$maxpathstrength = max(values(%pathstrength));
print "maxpathstrength = $maxpathstrength\n";
drawstrengths();

#create edges as defined by your adjacecncy matrix
for($i=0; $i<@nodes; $i++){
    for($j=$i+1; $j<@nodes; $j++){
		next if($edgesDone{"$i-$j"});
        if($admat[$i][$j] eq 1){
            #print "adding edge from $i to $j\n";
			$graph->add_edge($nodes[$i] => $nodes[$j], color => 'gray');
		}
    }
}



sub drawstrengths{
	foreach $k (keys(%pathstrength)){
		($i,$j) = split /-/, $k;
		$c = $pathstrength{$k} / $maxpathstrength;
		$m = sum(values(%pathstrength)) / scalar(values(%pathstrength)) / $maxpathstrength ;
		if($c > 3*$m/2) {
		#if($pathstrength{$k} == $maxpathstrength) {
			$graph->add_edge($i => $j, 	color => "0.0,$c,0.99", 
										label => $pathstrength{$k},
			                           	fontcolor =>"0.0,$c,0.0");
			$edgesDone{"$i-$j"}++;
			$edgesDone{"$j-$i"}++;
		}
	}
}

sub addPathStrength{
	while(@_){
		$a = shift;
		$b = shift;
		($b,$a) = ($a,$b) if($a > $b);
		$pathstrength{"$a-$b"} ++;
	}
}

sub drawpath{
	$color = shift;
	print "  drawing $color path...\n";
	while(@_){
		$a = shift;
		$b = shift;
		die "unsymm" unless ($admat[$a][$b] == $admat[$b][$a]);
		if( $admat[$a][$b] != 0){
			$graph->add_edge($a => $b, color => $color, arrowhead => 'normal'); 
		}else{
			print "    jumping from ".($a)." to ".($b)."\n";
		};
	};
};


#render the graph to a png-file
$graph->as_png("pretty.png");
