/*       Created   :  10/07/2008 09:04:23 PM
 *       Last Change: Sat Dec 20 01:00 PM 2008 CET
 */

#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

using namespace boost::program_options;
using namespace std;

struct SeriationGenCfg{
	SeriationGenCfg();
};

SeriationGenCfg::SeriationGenCfg(){
	options_description od("========== Seriation Generator Options ==========");
	od.add_options()
		("ser-gen.verbose",value<bool>()->default_value(false),"verbose seriation generation output")
		("ser-gen.sdp-wrapper,w",value<string>()->default_value("SDPAWrapper"),"Which SDP-Solver to use")
		("ser-gen.sdp-rand-plane-tries,t",value<int>()->default_value(10000),"How many random-hyperplane tries")
		("ser-gen.sdp-rounding-method,m",value<int>()->default_value(2),"How to round pos semdef Matrix to Vec")
		;
	options_description rsg("  Repetative Seriation Generator");
	rsg.add_options()
		("repetative-ser-gen.serializer", value<string>()->default_value("SDPSeriationGen"), "Which Serializer to call repeatedly")
		("repetative-ser-gen.repetition_num", value<int>()->default_value(3), "How often to call the Serializer")
		;
	options_description gsg("  GDist Seriation Generator");
	gsg.add_options()
		("gdist-ser-gen.seriation_gen", value<string>()->default_value("LEVSeriationGen"), "Serializer used to determine start vertex in graph")
		("gdist-ser-gen.normalize", value<bool>()->default_value(true), "Whether all distances should be confined btw [0,1]")
		;

	gCfg().addModuleOptions(od);
	gCfg().addModuleOptions(rsg);
	gCfg().addModuleOptions(gsg);
}
namespace{
	SeriationGenCfg _tmp;
}

