#include <dlfcn.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include "configuration.hpp"
#include "action.hpp"

#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C
using namespace boost;
using namespace std;


int main(int argc, char* argv[]){
	void* err;
	err = dlopen("actions/libactions.so",             RTLD_LAZY); I(err);
	err = dlopen("sdp_wrappers/libsdp_wrappers.so",   RTLD_LAZY); I(err);

	gCfg().parsecfg(argc,argv);

	string action_name = gCfg().getString("action");
	auto_ptr<Action> action = genericFactory<Action>::instance().create(action_name);
	if(!action.get())
		throw logic_error(string("Supplied action `") + action_name + "' does not exist");
	(*action)();
}
