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
	options_description od("Adjacency Matrix Generation");

	options_description rand_amg("Random AdjMatGen");
	rand_amg.add_options()
		("rand_adj_mat_gen.size,n", value<int>()->default_value(20), "Size of the Matrix to be generated")
		("rand_adj_mat_gen.out-degree,p", value<float>()->default_value(0.4f), "Average out-degree of a vertex")
		("rand_adj_mat_gen.seed,s", value<float>()->default_value(1.0f), "Seed for random number generator")
		("rand_adj_mat_gen.weighted", value<bool>()->default_value(false), "Whether edges should be weighted")
		;
	od.add(rand_amg);
	gCfg().addModuleOptions(od);
}

namespace {
	AdjMatGenCfg _adjmatgencfg;
}
