#include <boost/shared_ptr.hpp>
#include <fstream>
#include <factory/factory.h>
#include <cstdio>
#include <cstring>

#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;

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

#include <graph_embedder.hpp> 
#include <feature_embedder.hpp>
#include <anndb.hpp>
#include <pointnn.hpp>

#include "build_db.hpp"

#include <postproc.hpp>
#include <stats.hpp>

#include <gnuplot_i.hpp>

#include <nana.h>


using namespace std;
using namespace boost;


void BuildDB::operator()()
{
	//mDB.reset( new ANNDB() );
	mDB.reset( new PointNNDB() );
	mDB->configure();

	string adjmat_gen_name               = gCfg().getString("serialize.adjmat_gen");
	shared_ptr<AdjMatGen> adjmat_gen     ( genericFactory<AdjMatGen>::instance().create(adjmat_gen_name) );
	if(!adjmat_gen.get())
		throw logic_error(string("Supplied AdjMatGen `") + adjmat_gen_name + "' does not exist");
	adjmat_gen->configure();

	string postproc_name                 = gCfg().getString("output-format");
	auto_ptr<PostProc> out_ptr = genericFactory<PostProc>::instance().create(postproc_name);
	if(!out_ptr.get())
		throw logic_error(string("Supplied Postprocessor `") + postproc_name + "' does not exist");

	int max_num = gCfg().getInt("serialize.max_num");

	//out_ptr->atStart();
	bool verbose    = gCfg().getBool("verbose");
	bool nonverbose = !gCfg().getBool("quiet") && !gCfg().getBool("verbose");
	int embed_dim = gCfg().getInt("BuildDB.embed_dim");

	ofstream classAcc("/tmp/acc.dat");

	string gemb_name                     = gCfg().getString("BuildDB.embed_meth");
	auto_ptr<GraphEmbedder> gemb_ptr     = genericFactory<GraphEmbedder>::instance().create(gemb_name);
	if(!gemb_ptr.get())
		throw logic_error(string("Supplied Embedding Method `") + gemb_name + "' does not exist");
	if(gemb_name == "FeatureGraphEmbedder")
		((FeatureGraphEmbedder*)gemb_ptr.get())->setAdjMatGen(adjmat_gen, embed_dim);
	
	string ofile("centroids.dat");
	string opath =  gCfg().getString("output-dir");
	ofstream ofp( (fs::path(opath) / fs::path(ofile)).string().c_str());


	gemb_ptr->configure();

	bool iterate = false;
	do{
		int  cnt=0;
		ExactDescriptiveStatistics nstats("nstats");
		
		ProgressBar pb(max_num==0?600:max_num, "serializing");
		while(true){
			if(nonverbose) { pb.inc();                                          }
			if(verbose)    { L("Action::BuildDB %03d: Generating...\n", cnt); }
			ProbAdjPerm prob =  (*adjmat_gen)() ;
			if(!adjmat_gen->hasNext())
				break;
			nstats.notify(prob.getAdjMat()->size1());
			// WARNING: Permutation is dangerous, since modules later 
			// on do not have access to the permutation matrix!
			shared_ptr<PermMat::PermMatT> P_ptr(new PermMat::PermMatT(prob.getAdjMat()->size1(),prob.getAdjMat()->size1()));
			*P_ptr = numeric::ublas::identity_matrix<double>(prob.getAdjMat()->size1(),prob.getAdjMat()->size1());
			prob.setPermMat(P_ptr);

			if(verbose){ L("Action::BuildDB %03d: Building internal representation...\n", cnt);}
			GraphEmbedder::cloud_type cloud = (*gemb_ptr)(prob,embed_dim);
			//if(adjmat_gen->getClassID() == 1)
			mDB->add(prob, cloud);

			cnt++;
			if(cnt == max_num)
				break;
		}
		if(nonverbose) { pb.finish(); }
		cout <<nstats<<endl;
		mDB->finish();
		//out_ptr->atEnd();
	
		int evalMode = gCfg().getInt("BuildDB.evalmode");
	
		if(0);
		else if(evalMode == 1)
			spatialAnalysis(cnt, *adjmat_gen);
		else if(evalMode == 2){
			mDB->evaluate(*adjmat_gen, 5);
			ExactDescriptiveStatistics acc("Accuracy ");
			while(true){
				ProbAdjPerm pap = (*adjmat_gen)();
				if(pap.getAdjMat().get()==NULL)
					break;
				int tc = adjmat_gen->getClassID(pap.getBackground());
				int cc = mDB->classify((*gemb_ptr)(pap,embed_dim), *adjmat_gen, 10);
				acc += tc==cc;
			}
			cout << acc<<endl;
			classAcc << acc.getMean() <<endl;
			mDB->clear();
			if(gemb_name == "FeatureGraphEmbedder") // reset this, too
				((FeatureGraphEmbedder*)gemb_ptr.get())->setAdjMatGen(adjmat_gen, embed_dim);
			adjmat_gen->rewind();
			iterate = true;
		}else
			cout << "Unknown EvalMode"<<endl;
	}while(iterate);
}


BuildDB::~BuildDB()
{
}

int spatialHash(
		const ublas::vector<int>& pt,
		const ublas::vector<int>& res){
	int n = res.size();
	int hashval = 0;
	for(int i=n-1;i>=0;i--)
	{
		int prevprod=1;
		for(int j=i-1;j>=0;j--)
			prevprod *= res(j);
		hashval += pt(i) * prevprod;
	}
	return hashval;
}
int spatialHash(
		const ublas::vector<double>& minv, 
		const ublas::vector<double>& maxv,
		const ublas::vector<int> &   res,
		const ublas::vector<double> &pt){
	
	int n = pt.size();
	ublas::vector<int> npt(n);
	for(int i=0;i<n;i++){
		npt(i) = (int)(res(i) * (pt(i)-minv(i))/(maxv(i)-minv(i)));
	}

	return spatialHash(npt, res);
}

ublas::vector<double> f(const ublas::vector<double>& v){
	ublas::vector<double> w(v);
	//double d = ublas::norm_2(v);
	//if(d>0.00000000001)
		//w /= d*d;
	return w;
}

void BuildDB::spatialAnalysisCloud(int cnt, AdjMatGen& adjmatgen)
{
	ofstream os("/tmp/spatdb.csv");
	for(int c=0;c<cnt;c++){
		const GraphDB::TCloud& cloud = mDB->getCloud(c);
		string s                     = mDB->getPap(c).getId();
		int klass                    = adjmatgen.getClassID(s);
		if(c==0){
			// header
			for(unsigned int i=0;i<cloud[0].size()*cloud.size();i++)
				os<<"c"<<i<<",";
			os<<"klass"<<endl;
		}
		for(unsigned int i=0;i<cloud.size();i++){
			const ublas::vector<double>& c = cloud[i];
			for(unsigned int j=0;j<c.size();j++){
				double d = c[i];
				if(fabs(d)>1E6) // TODO: woher kommen die groszen zahlen?
					d=0.0;
				if(fabs(d)<1E-20) 
					d=0.0;
				os << d<<",";
			}
		}
		os << "c"<<klass<<endl;
	}
}

void BuildDB::spatialAnalysis(int cnt, AdjMatGen& adjmatgen)
{  
	ofstream os("/tmp/spatdb.csv");
	for(int c=0;c<cnt;c++){
		ublas::vector<int> v  = mDB->getFeatures(c);
		const boost::any& ref = mDB->getPap(c).getBackground();
		int klass             = adjmatgen.getClassID(ref);
#define NUMATOM_ONLY 0
#if NUMATOM_ONLY
		if(c==0)
			os << "numAtom,klass"<<endl;
		os << mDB->getPap(c).getAdjMat(ref)->size1()<<",";
#else
		if(c==0){
			// header
			for(unsigned int i=0;i<v.size();i++)
				os<<"c"<<i<<",";
			os<<"klass"<<endl;
		}
		for(unsigned int i=0;i<v.size();i++)
			os<<v(i)<<",";
#endif
		os<<"c"<<klass<<endl;
	}
}


namespace{ registerInFactory<Action, BuildDB> registerBase("BuildDB"); }
