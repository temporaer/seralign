#include <dlfcn.h>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <factory/factory.h>
#include <cstdio>

#include <configuration.hpp>
#include <sdp_wrapper.hpp>
#include <adjmat_gen.hpp>


#include <SerGenAdj.hpp>

#include <DegreeSort.hpp>

#include <matlab_io.hpp>

#include "serialize.hpp"

#include <postproc.hpp>

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

	string postproc_name                 = gCfg().getString("output-format");
	auto_ptr<PostProc> out_ptr = genericFactory<PostProc>::instance().create(postproc_name);
	if(!out_ptr.get())
		throw logic_error(string("Supplied Postprocessor `") + postproc_name + "' does not exist");
	
	int max_num = gCfg().getInt("serialize.max_num");

	string seriation_gen_name            = gCfg().getString("serialize.seriation_gen");
	auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(seriation_gen_name);
	if(!seriation_gen.get())
		throw logic_error(string("Supplied SerGenAdj `") + seriation_gen_name + "' does not exist");
	seriation_gen->configure();

	out_ptr->atStart();
	int cnt=0;
	while(true){
		L("Action::Serialize %03d: Generating...\n", cnt);
		ProbAdjPerm tmp =  (*adjmat_gen)() ;
		ProbAdjLapPerm prob ( tmp );
		if(!adjmat_gen->hasNext())
			break;

#if 1
		DegreeSort ds;
		ds.sort(prob);
#else
		shared_ptr<PermMat::PermMatT> P_ptr(new PermMat::PermMatT(prob.getAdjMat()->size1(),prob.getAdjMat()->size1()));
		*P_ptr = numeric::ublas::identity_matrix<double>(prob.getAdjMat()->size1(),prob.getAdjMat()->size1());
		prob.setPermMat(P_ptr);
#endif
		prob.calculateLaplacian();

		L("Action::Serialize %03d: Serializing...\n", cnt);
		Serialization randwalk;
		randwalk = (*seriation_gen)(prob);

		L("Action::Serialize %03d: Postprocessing...\n", cnt);
		out_ptr->atSeriation(*adjmat_gen, randwalk, prob);

		cnt++;
		if(cnt == max_num)
			break;
	}
	out_ptr->atEnd();
}

Serialize::~Serialize()
{
}

namespace{ registerInFactory<Action, Serialize> registerBase("Serialize"); }
