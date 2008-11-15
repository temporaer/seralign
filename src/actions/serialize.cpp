#include <dlfcn.h>
#include <boost/shared_ptr.hpp>
#include <factory/factory.h>
#include <cstdio>

#include <configuration.hpp>
#include <sdp_wrapper.hpp>
#include <adjmat_gen.hpp>

#include "serialize.hpp"
#include <nana.h>

using namespace std;
using namespace boost;

void Serialize::operator()()
{
	string adjmat_gen_name               = gCfg().getString("serialize.adjmat_gen");
	auto_ptr<AdjMatGen> adjmat_gen       = genericFactory<AdjMatGen>::instance().create(adjmat_gen_name);
	if(!adjmat_gen.get())
		throw logic_error(string("Supplied AdjMatGen `") + adjmat_gen_name + "' does not exist");
	adjmat_gen->configure();

	shared_ptr<AdjMatGen::AdjMatT> adjmat_ptr((*adjmat_gen)());

	//string seriation_gen_name            = gCfg().getString("serialize.seriation_gen");
	//auto_ptr<SeriationGen> seriation_gen = genericFactory<SeriationGen>::instance().create(seriation_gen_name);
	//if(!seriation_gen.get())
		//throw logic_error(string("Supplied SeriationGen `") + seriation_gen_name + "' does not exist");

	//seriation_gen->configure();

	//SeriationGen::SeriationT randwalk;
	//randwalk = (*seriation_gen)(adjmat_ptr);
}

Serialize::~Serialize()
{
}

namespace{
	registerInFactory<Action, Serialize> registerBase("Serialize");
}
