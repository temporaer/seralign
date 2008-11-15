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
		;
	od.add(action);
	gCfg().addModuleOptions(od);
}

namespace {
	ActionCfg _actioncfg;
}
