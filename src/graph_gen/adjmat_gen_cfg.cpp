#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class AdjMatGenCfg{
	public:
		AdjMatGenCfg();
};

AdjMatGenCfg::AdjMatGenCfg(){
	options_description od("========== Adjacency Matrix Generation ==========");
	od.add_options()
		("adjmat_gen.verbose", value<bool>()->default_value(false), "Whether to be verbose")
		;

	options_description mutagenesis("  Mutagenesis");
	mutagenesis.add_options()
		("mutagenesis.kernel",                   value<string>()->default_value("kernelNull"), "Kernel to use")
		("mutagenesis.in-file",                  value<string>(), "file with mutagenesis CSV data")
		("mutagenesis.num-seriation-neighbours", value<int>()->default_value(0), "how many neighbours to each side in seriation to include in output")
		("mutagenesis.include-graph-neighbours", value<bool>()->default_value(true), "whether to include graph neighbours in output")
		;
	options_description rand_amg("  Random AdjMatGen");
	rand_amg.add_options()
		("rand_adj_mat_gen.size", value<int>()->default_value(20), "Size of the Matrix to be generated")
		("rand_adj_mat_gen.out-degree", value<float>()->default_value(3.0f), "Average out-degree of a vertex")
		("rand_adj_mat_gen.seed", value<float>()->default_value(1.0f), "Seed for random number generator")
		("rand_adj_mat_gen.weighted", value<bool>()->default_value(false), "Whether edges should be weighted")
		;
	options_description randpat_amg("  Random Pattern AdjMatGen");
	randpat_amg.add_options()
		("rand_pat_adjmat_gen.size",    value<int>()->default_value(30), "Size of random pattern matrix to be generated")
		("rand_pat_adjmat_gen.pattern_size", value<int>()->default_value(10), "Size of pattern to be generated")
		("rand_pat_adjmat_gen.edge_inclusion_prob", value<float>()->default_value(0.5f), "Probability of an edge not part of pattern to be included")
		("rand_pat_adjmat_gen.gen_pat", value<bool>()->default_value(true), "If false, randomly assign nodes the part_of_pattern-property")
		;
	options_description fullconn_amg("  FullConn AdjMatGen");
	fullconn_amg.add_options()
		("fullconn_adjmat_gen.size", value<int>()->default_value(15), "Size of Graph to generate")
		("fullconn_adjmat_gen.patsize",  value<int>()->default_value(5), "Size of Pattern in Graph to generate")
		("fullconn_adjmat_gen.jumble",  value<bool>()->default_value(false), "Whether to jumble the matrices")
		;
	options_description sdf_amg("  SDF AdjMatGen");
	sdf_amg.add_options()
		("SDFAdjmatGen.files", value<string>()->default_value(""), "sdf-files, format: file:class,file:class")
		("SDFAdjmatGen.fixed_size", value<int>()->default_value(0), "fixed size of molecules")
		;
	od.add(rand_amg);
	od.add(randpat_amg);
	od.add(fullconn_amg);
	od.add(mutagenesis);
	od.add(sdf_amg);
	gCfg().addModuleOptions(od);
}

namespace {
	AdjMatGenCfg _adjmatgencfg;
}
