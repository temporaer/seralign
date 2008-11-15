/*       Created   :  10/07/2008 12:08:01 PM
 *       Last Change: Mon Oct 27 11:00 PM 2008 CET
 */

#include <configuration.hpp>
#include <boost/program_options.hpp>
#include <factory/factory.h>

struct SDPWrapperCfg{
	SDPWrapperCfg();
};

using namespace boost::program_options;
using namespace std;
SDPWrapperCfg::SDPWrapperCfg(){
	options_description OD;
#ifdef HAVE_DSDP
	options_description dsdp("DSDP Options");
	//dsdp.add_options()
		//("f",value<float>(),"f, a float.")
		//;
	OD.add(dsdp);
#endif
#ifdef HAVE_CSDP
	options_description csdp("CSDP Options");
	//csdp.add_options()
		//("g",value<float>(),"g, a float.")
		//;
	OD.add(csdp);
#endif
#ifdef HAVE_SDPA
    options_description sdpa("SDPA Options");
    sdpa.add_options()
        ("sdpa.param-file",value<string>()->default_value(SDPA_PARAM_FILE), "SDPA Parameter File")
        ;
    OD.add(sdpa);
#endif
#ifdef HAVE_SDPLR
    options_description sdplr("SDPLR Options");
    sdplr.add_options()
        ("sdplr.param-file",value<string>()->default_value(SDPLR_PARAM_FILE), "SDPLR Parameter File")
        ;
    OD.add(sdplr);
#endif
	gCfg().addModuleOptions(OD);
}
namespace{
	SDPWrapperCfg _sdpwrappercfg;
}
