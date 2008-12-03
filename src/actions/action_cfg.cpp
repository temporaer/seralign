#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

class ActionCfg{
	public:
		ActionCfg();
};

ActionCfg::ActionCfg(){
	options_description od("Action Options");

	options_description action("Serialize Options");
	action.add_options()
		("serialize.adjmat_gen",    value<string>()->default_value("RandomAdjMatGen"), "Adjacency Matrix Generator")
		("serialize.seriation_gen", value<string>()->default_value("SDPSeriationGen"), "Seriation Generator")
		("serialize.max_num", value<int>()->default_value(0), "Maximum number of graphs to serialize")
		;
	od.add(action);
	gCfg().addModuleOptions(od);
}

namespace {
	ActionCfg _actioncfg;
}
