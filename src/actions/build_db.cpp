#include <boost/shared_ptr.hpp>
#include <fstream>
#include <factory/factory.h>
#include <cstdio>
#include <cstring>

#include <boost/numeric/ublas/matrix.hpp>

#include <configuration.hpp>
#include <adjmat_gen.hpp>

#include <progressbar.hpp>
#include <SerGenAdj.hpp>

#include <DegreeSort.hpp>
#include <GraphFromAdj.hpp>

#include "build_db.hpp"

#include <postproc.hpp>

#include <nana.h>

using namespace std;
using namespace boost;


void BuildDB::operator()()
{
	mDB.configure();

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

	bool wantDegreeSorting = gCfg().getBool("serialize.want_degree_sorting");

	out_ptr->atStart();
	int  cnt=0;
	bool verbose    = gCfg().getBool("verbose");
	bool nonverbose = !gCfg().getBool("quiet") && !gCfg().getBool("verbose");
	
	ProgressBar pb(max_num==0?200:max_num, "serializing");
	while(true){
		if(nonverbose) { pb.inc();                                          }
		if(verbose)    { L("Action::BuildDB %03d: Generating...\n", cnt); }
		ProbAdjPerm tmp =  (*adjmat_gen)() ;
		ProbAdjLapPerm prob ( tmp );
		if(!adjmat_gen->hasNext())
			break;
		if(wantDegreeSorting){
			DegreeSort ds;
			ds.sort(prob);
		}else{
			shared_ptr<PermMat::PermMatT> P_ptr(new PermMat::PermMatT(prob.getAdjMat()->size1(),prob.getAdjMat()->size1()));
			*P_ptr = numeric::ublas::identity_matrix<double>(prob.getAdjMat()->size1(),prob.getAdjMat()->size1());
			prob.setPermMat(P_ptr);
		}
		prob.calculateLaplacian();

		if(verbose){ L("Action::BuildDB %03d: Building internal representation...\n", cnt);}
		GDistProjectedDB::TCloud cloud = mDB.add(adjmat_gen->getGraphID(), prob);
		out_ptr->atSeriation(*adjmat_gen, cloud, prob);

		cnt++;
		if(cnt == max_num)
			break;
	}
	if(nonverbose) { pb.finish(); }
	out_ptr->atEnd();

	ProgressBar matching(cnt, "Matching");
	ofstream os(gCfg().getOutputFile("output").c_str());
	// output distances as CSV
	for(vector<string>::const_iterator it=mDB.getIDs().begin(); it!=mDB.getIDs().end(); it++){
		os <<","<<(*it);
	}
	os<<endl;
	numeric::ublas::matrix<double> dmat(cnt,cnt);
	int dmat_i=0, dmat_j=0;
	for(GDistProjectedDB::iterator it=mDB.begin(); it!=mDB.end(); it++){
		if(nonverbose){matching.inc();}
		vector<GDistProjectedDB::point_type> query;
		copy(it->begin(),it->end(),back_inserter(query));
		os << "dummy";
		for(GDistProjectedDB::iterator it2=mDB.begin(); it2!=mDB.end(); it2++){
			it2->match(query.begin(),query.end(),
					GDistProjectedDB::point_type::Translator(),
					GDistProjectedDB::point_type::Rotator<GDistProjectedDB::TICP::TQuat>());
			dmat(dmat_i,dmat_j) = it2->getMatchingError(); 
			dmat_i++;
		}
		dmat_j++; dmat_i=0;
	}
	if(nonverbose){matching.finish();}
	for(int i=0;i<cnt;i++){
		os << "dummy";
		for(int j=0;j<cnt;j++){
			os << "," <<min(dmat(i,j),dmat(j,i));
		}
		os << endl;
	}
}

BuildDB::~BuildDB()
{
}

namespace{ registerInFactory<Action, BuildDB> registerBase("BuildDB"); }
