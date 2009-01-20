#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class GraphEmbedderCfg{
	public:
		GraphEmbedderCfg();
};

GraphEmbedderCfg::GraphEmbedderCfg(){
	options_description od("========== Action Options ==========");

	options_description hkge("  HeatkernelGraphEmbedder Options");
	hkge.add_options()
		("HeatkernelGraphEmbedder.t",    value<float>()->default_value(1.0), "locality parameter")
		;
	od.add(hkge);
	gCfg().addModuleOptions(od);
}

namespace {
	GraphEmbedderCfg _cfg;
}
