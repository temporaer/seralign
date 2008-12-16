#include <boost/shared_ptr.hpp>
#include <string>
#include <fstream>
#include <vector>
#include <factory/factory.h>
#include <cstdio>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <utility/fastmap.h>
#include <stats.hpp>
#include <gnuplot_i.hpp>
#include <map>

#include <configuration.hpp>
#include <fastmap_distancemat.hpp>
#include <nana.h>

using namespace std;
using namespace boost;
using namespace utility;
namespace ublas = boost::numeric::ublas;

typedef ublas::matrix<double> DistMatT;
typedef shared_ptr<DistMatT>  DistMatP;

struct DistMatCSV{
	vector<string> headers;
	DistMatP       mat;
};

class D : public abstract_distance_functor<vector<float>, float>{
	public:
	DistMatP     mDistMat;
	double       mMin,mMax;

	D(DistMatP p) : mDistMat(p) {
		ExactDescriptiveStatistics stats("diststats");
		unsigned int n = mDistMat->size1();

		for(unsigned int i=0;i<mDistMat->size1();i++)
			for(unsigned int j=0;j<mDistMat->size2();j++) {
				if(i==j) continue;
				stats += (*mDistMat)(i,j);
			}
		// make sure dist between 1 and 2 (for log())
		(*mDistMat) += ublas::scalar_matrix<double>(n,n,-stats.getMin());
		(*mDistMat) /= stats.getRange();
		//(*mDistMat) += ublas::scalar_matrix<double>(n,n,1.0);
		
		stats.reset();
		for(unsigned int i=0;i<mDistMat->size1();i++)
			for(unsigned int j=0;j<mDistMat->size2();j++)
			{
				if(i==j) continue;
				(*mDistMat)(i,j) = exp(-(*mDistMat)(i,j));
				stats += (*mDistMat)(i,j);
			}
		cout << stats<<endl;
		(*mDistMat) += ublas::scalar_matrix<double>(n,n,-stats.getMin());
		(*mDistMat) /= stats.getRange();
		mMin = 0;
		mMax = 1;
	}
	virtual bool is_symmetric()const {return true;}
	virtual bool zero_means_same() const {return false;}
	double getDist(int i, int j)const{
		if(i==j) return mMin;
		return 1-exp(-(*mDistMat)(i,j));
	}
	virtual score_type operator() (const feature_type &i, const feature_type &j) const{
		unsigned int a = round(i[0]);
		unsigned int b = round(j[0]);
		return   getDist(a,b);
	}
};

DistMatCSV
readDistanceMatrixCSVWithNameCol(const char* fn){
	ifstream is(fn);
	string line;
	getline(is,line);
	trim(line);
	DistMatCSV csv;
	vector<string> headers;
	split( headers, line, boost::is_any_of(",") );
	copy(headers.begin()+1, headers.end(), back_inserter(csv.headers));

	unsigned int n = csv.headers.size();
	DistMatP distmat_ptr( new DistMatT(n,n) );
	DistMatT& distmat = *distmat_ptr;

	vector<string> vals;
	for(unsigned int i=0;i<n; i++){
		getline(is,line); trim(line);
		split( vals, line, boost::is_any_of(",") );
		I(vals.size() == n+1);
		for(unsigned int j=1;j<=n;j++){
			if(vals[j] == "None")
				distmat(i,j-1) = 0;
			else
			{
				//L("Converting %s to double\n", vals[j].c_str());
				distmat(i,j-1) = boost::lexical_cast<double>(vals[j]);
			}
		}
	}
	csv.mat = distmat_ptr;
	is.close();
	return csv;;
}

void FastmapDistancemat::operator()()
{
	Gnuplot gp;
	vector<string> dmfns = gCfg().get<vector<string> >("fastmap.matrices");
	BOOST_FOREACH(string& dmfn, dmfns){
		DistMatCSV dm = readDistanceMatrixCSVWithNameCol(dmfn.c_str());
		D d(dm.mat);
		fastmap<vector<float> > fm;
		fm.set_distance_function(&d);
		vector<vector<float> > dataset;
		for(unsigned int i=0;i<dm.mat->size1(); i++)
			dataset.push_back(vector<float>(1,i));
		fm.make_map(2,dataset);
		vector<vector<float> > mapped = fm.get_map();

		map<string, int> classes;
		ifstream is(gCfg().getString("fastmap.classlabels").c_str());
		while(is.good()){
			string s; int i;
			is >> s >> i;
			classes[s] = i;
		}

		string rewhat = gCfg().getString("fastmap.filename_replace_what");
		string rewith = gCfg().getString("fastmap.filename_replace_with");
		string outfilename = regex_replace(dmfn, boost::regex(rewhat), rewith, boost::match_default | boost::format_sed);
		string plotfile    = regex_replace(dmfn, boost::regex(rewhat), "", boost::match_default | boost::format_sed);
		L("Saving fastmap to %s\n", outfilename.c_str());
		ofstream os(outfilename.c_str());
		for(unsigned int i=0;i<mapped.size(); i++){
			copy(mapped[i].begin(), mapped[i].end(), ostream_iterator<float>(os,"\t"));
			os << classes[dm.headers[i]];
			os << endl;
		}
		os.close();
		gp.cmd("set title \"Fastmap of Distance Matrix\"");
		gp.cmd("set terminal png giant size 1024,768");
		gp.cmd("set view 0,0");
		gp.cmd("unset border");
		gp.cmd("unset colorbox");
		gp.cmd("set pointsize 1 ");
		gp.cmd("set palette defined ( 10 \"#ff0000\", 90 \"#0000ff\" )");
		gp.cmd("set output \"" + plotfile + ".png\"");
		gp.cmd("splot \""+outfilename+"\" using 1:2:3 w p pt 7 palette t \"" + dmfn +"\"");
	}
}

FastmapDistancemat::~FastmapDistancemat()
{
}

namespace{ registerInFactory<Action, FastmapDistancemat> registerBase("FastmapDistancemat"); }
