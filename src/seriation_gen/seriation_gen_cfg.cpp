/*       Created   :  10/07/2008 09:04:23 PM
 *       Last Change: Sat Nov 15 09:00 PM 2008 CET
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
	options_description od("Seriation Generator Options");
	od.add_options()
		("ser-gen.method,g", value<string>()->default_value("SDP"), "How to generate seriation")
		("ser-gen.sdp-wrapper,w",value<string>()->default_value("SDPAWrapper"),"Which SDP-Solver to use")
		("ser-gen.sdp-rand-plane-tries,t",value<int>()->default_value(10000),"How many random-hyperplane tries")
		("ser-gen.sdp-rounding-method,m",value<int>()->default_value(2),"How to round pos semdef Matrix to Vec")
		;
	gCfg().addModuleOptions(od);
	gCfg().dependent_options("ser-gen.method","ser-gen.sdp-wrapper");
}
namespace{
	SeriationGenCfg _tmp;
}

