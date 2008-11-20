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
#include <nana.h>

using namespace std;
using namespace boost;

void prologhead(ostream& os){ os  << "<sequences>"<<endl; }
void prologfoot(ostream& os){ os  << "</sequences>"<<endl; }
void prologseqhead(ostream& os, AdjMatGen& gen){ os  << "  <sequence id=\"" << gen.getGraphID() <<"\">"<<endl; }
void prologseqfoot(ostream& os){ os  << "  </sequence>"<<endl; }
void prologout(int idx,ostream& os,AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob){
	prologseqhead(os,gen);
	for(int i=0; i<ser.size(); i++)
	{
		int idx = prob.getOriginalIndex(ser[i]);
		os	<<"    <atom>" <<endl
			<<"      <symbol>" << gen.getPrologDescription(idx) << "</symbol>" << endl
			<<"      <label>none</label>"   <<endl
			<<"    </atom>"<<endl;
	}
	prologseqfoot(os);
}

void Serialize::operator()()
{
	ofstream os(gCfg().getString("output").c_str());
	string adjmat_gen_name               = gCfg().getString("serialize.adjmat_gen");
	auto_ptr<AdjMatGen> adjmat_gen       = genericFactory<AdjMatGen>::instance().create(adjmat_gen_name);
	if(!adjmat_gen.get())
		throw logic_error(string("Supplied AdjMatGen `") + adjmat_gen_name + "' does not exist");
	adjmat_gen->configure();

	prologhead(os);
	int cnt=0;
	while(true){
		ProbAdjPerm tmp =  (*adjmat_gen)() ;
		ProbAdjLapPerm prob ( tmp );
		if(!adjmat_gen->hasNext())
			break;

		DegreeSort ds;
		ds.sort(prob);
		prob.calculateLaplacian();

		string seriation_gen_name            = gCfg().getString("serialize.seriation_gen");
		auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(seriation_gen_name);
		if(!seriation_gen.get())
			throw logic_error(string("Supplied SerGenAdj `") + seriation_gen_name + "' does not exist");

		seriation_gen->configure();

		Serialization randwalk;
		randwalk = (*seriation_gen)(prob);

		prologout(cnt, os,*adjmat_gen, randwalk, prob);

		cnt++;
	}
	prologfoot(os);
}

Serialize::~Serialize()
{
}

namespace{ registerInFactory<Action, Serialize> registerBase("Serialize"); }
