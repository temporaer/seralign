#include <boost/shared_ptr.hpp>
#include <fstream>
#include <factory/factory.h>
#include <cstdio>
#include <cstring>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
namespace bl = boost::lambda;

#include <boost/filesystem/path.hpp>
namespace fs = boost::filesystem;

#include <boost/lambda/casts.hpp>

#include <configuration.hpp>
#include <adjmat_gen.hpp>

#include <progressbar.hpp>
#include <SerGenAdj.hpp>

#include <DegreeSort.hpp>
#include <GraphFromAdj.hpp>

#include "build_db.hpp"

#include <postproc.hpp>
#include <stats.hpp>

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

	//out_ptr->atStart();
	int  cnt=0;
	bool verbose    = gCfg().getBool("verbose");
	bool nonverbose = !gCfg().getBool("quiet") && !gCfg().getBool("verbose");
	
	string ofile("centroids.dat");
	string opath =  gCfg().getString("output-dir");
	ofstream ofp( (fs::path(opath) / fs::path(ofile)).string().c_str());
	
	ProgressBar pb(max_num==0?200:max_num, "serializing");
	while(true){
		if(nonverbose) { pb.inc();                                          }
		if(verbose)    { L("Action::BuildDB %03d: Generating...\n", cnt); }
		ProbAdjPerm tmp =  (*adjmat_gen)() ;
		ProbAdjLapPerm prob ( tmp );
		if(!adjmat_gen->hasNext())
			break;
		// WARNING: Permutation is dangerous, since modules later 
		// on do not have access to the permutation matrix!
		shared_ptr<PermMat::PermMatT> P_ptr(new PermMat::PermMatT(prob.getAdjMat()->size1(),prob.getAdjMat()->size1()));
		*P_ptr = numeric::ublas::identity_matrix<double>(prob.getAdjMat()->size1(),prob.getAdjMat()->size1());
		prob.setPermMat(P_ptr);
		prob.calculateLaplacian();

		if(verbose){ L("Action::BuildDB %03d: Building internal representation...\n", cnt);}
		GDistProjectedDB::TCloud cloud        = mDB.addSpectral(adjmat_gen->getGraphID(), prob);
		//GDistProjectedDB::TCloud cloud        = mDB.add(adjmat_gen->getGraphID(), prob);
		GDistProjectedDB::TICP::TVec centroid = mDB.back().getModelCentroid();
		ofp << centroid[0] <<"\t"<<centroid[1]<<"\t"<<centroid[2]<<"\t"<<adjmat_gen->getClassID()<<endl;
		

		//out_ptr->atSeriation(*adjmat_gen, cloud, prob);

		cnt++;
		if(cnt == max_num)
			break;
	}
	if(nonverbose) { pb.finish(); }
	//out_ptr->atEnd();

#if 0
	printDistMatrix(cnt);
#else
	int qid = gCfg().getInt("BuildDB.query_id");
	ExactDescriptiveStatistics stats("accuracy");
	for(int i=0;i<cnt;i++){
		int klass = knn_classify(cnt, i, *adjmat_gen, *out_ptr, 3);
		stats.notify(klass == adjmat_gen->getClassID(mDB.getIDs()[i]) ? 1 : 0);
	}
	cout <<stats<<endl;
#endif
}

int BuildDB::knn_classify(int cnt, int id, AdjMatGen& gen, PostProc& out, int k){
	//bool verbose    = gCfg().getBool("verbose");
	bool nonverbose = !gCfg().getBool("quiet") && !gCfg().getBool("verbose");

	ProgressBar matching(cnt, "Matching");
	ofstream os(gCfg().getOutputFile("output").c_str());
	std::vector<double>  dists(cnt); // icp distances
	std::vector<double> cdists(cnt); // centroid distances

	for(int i=0; i<cnt; i++){
		if(nonverbose){matching.inc();}
		cdists[i] = norm_2(mDB[i].getModelCentroid() - mDB[id].getModelCentroid());
		dists[i]  = match(mDB[i],mDB[id]);
	}
	if(nonverbose){matching.finish();}

	vector<unsigned int> granks(cnt), cranks(cnt);
	for(int i=0; i<cnt; i++) granks[i]=i;
	for(int i=0; i<cnt; i++) cranks[i]=i;
	sort(granks.begin(), granks.end(), bl::var( dists)[bl::_1] < bl::var( dists)[bl::_2]);
	sort(cranks.begin(), cranks.end(), bl::var(cdists)[bl::_1] < bl::var(cdists)[bl::_2]);
	cout << "Cranks: "<<endl;
	for(int i=0;i<k+1;i++){
		int idx = cranks[i];
		GDistProjectedDB::TICP& icp = mDB[idx];
		int n=icp.size();
		Serialization ser(n);
		int nodecount=0;
		for(GDistProjectedDB::TICP::iterator it=icp.begin(); it!=icp.end(); it++, nodecount++){
			ser.getRanks()[nodecount]   = it->id;
			ser.setPosition(it->id, it->pos);
		}
		cout << "  Match i="<<i<<"\t"<<mDB.getIDs()[idx]<<endl;
		out.atSeriation(gen, ser, mDB.getIDs()[idx]);
	}

	cout << "ICPranks: "<<endl;
	out.atStart();
	map<int,double> classcnt;
	typedef map<int,double>::iterator CCIt;
	typedef pair<int,double> CCItP;
	for(int i=0;i<k+1;i++){
		int idx = granks[i];
		double dist = dists[idx];
		GDistProjectedDB::TICP& icp = mDB[idx];
		int n=icp.size();
		Serialization ser(n);
		int nodecount=0;
		for(GDistProjectedDB::TICP::iterator it=icp.begin(); it!=icp.end(); it++, nodecount++){
			ser.getRanks()[nodecount]   = it->id;
			ser.setPosition(it->id, it->pos);
		}
		cout << "  Match i="<<i<<"\t"<<mDB.getIDs()[idx]<<endl;
		out.atSeriation(gen, ser, mDB.getIDs()[idx]);
		if(i>0) // TODO: need more elegant sol'n. ignore pattern itself!
		{
			int klass = gen.getClassID(mDB.getIDs()[idx]);
			if(classcnt.find(klass) == classcnt.end())
				classcnt[klass]  = 1.0/dist;
			else
				classcnt[klass] += 1.0/dist;
		}
	}
	out.atEnd();

	// return the class which was voted for most
	return max_element(classcnt.begin(), classcnt.end(), 
			bl::bind<const int&>(&CCItP::second, bl::_1) < 
			bl::bind<const int&>(&CCItP::second, bl::_2))
			->first;
}
void BuildDB::printDistMatrix(int cnt){
	//bool verbose    = gCfg().getBool("verbose");
	bool nonverbose = !gCfg().getBool("quiet") && !gCfg().getBool("verbose");

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
			if(dmat_i < dmat_j){
				dmat_i++;
				continue;
			}
			double dist = match(*it, *it2);
			dmat(dmat_i,dmat_j) = dmat(dmat_j,dmat_i) = dist;
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
double
BuildDB::match(GDistProjectedDB::TICP& a,GDistProjectedDB::TICP& b){
	vector<GDistProjectedDB::point_type> query1, query2;
	copy(a.begin(),a.end(),back_inserter(query1));
	copy(b.begin(),b.end(),back_inserter(query2));
	b.match(query1.begin(),query1.end(),
			GDistProjectedDB::point_type::Translator(),
			GDistProjectedDB::point_type::Rotator<GDistProjectedDB::TICP::TQuat>());
	a.match(query2.begin(),query2.end(),
			GDistProjectedDB::point_type::Translator(),
			GDistProjectedDB::point_type::Rotator<GDistProjectedDB::TICP::TQuat>());
	return a.getMatchingError() + b.getMatchingError();
}

BuildDB::~BuildDB()
{
}


namespace{ registerInFactory<Action, BuildDB> registerBase("BuildDB"); }
