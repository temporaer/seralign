#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <factory/factory.h>
#include <GraphFromAdj.hpp> 
#include <configuration.hpp>
//#include <matlab_io.hpp>
#include "gdist_seriation_gen.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;

// 2008-12-20 Hannes Schulz <mail at hannes-schulz dot de>

GDistSeriationGen::GDistSeriationGen()
{
}

GDistSeriationGen::~GDistSeriationGen()
{
  // cleanup
}

struct PositionSorter{
	PositionSorter(const Serialization::PosT& p)
		:mPos(p)
	{
	}
	inline bool operator()(int i, int j)const{
		return mPos(i) < mPos(j);
	}
	const Serialization::PosT&  mPos;
};

Serialization GDistSeriationGen::operator()(const ProbAdjPerm& pap)
{
	auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(mSeriationGenName);
	if(!seriation_gen.get())
		throw logic_error(string("Supplied SerGenAdj `") + mSeriationGenName + "' does not exist");
	seriation_gen->configure();

	Serialization ser = (*seriation_gen)(pap);

	// drop out if serialization empty
	int n = ser.getRanks().size();
	if(n==0) return ser;

	// create graph, calculate distances starting at 1st node
	GraphFromAdj g(pap,ser.getRanks()[0], ser.getRanks()[n-1]);

	Serialization::PosT  pos(n);
	Serialization::RankT ranks(n);
	for(int i=0;i<n;i++){
		pos(i)   = g.getDist(i);
#ifndef NDEBUG
		if(pos(i)!=pos(i))
			cerr <<"Warning: Distance is NaN for index "<<i<< endl;
#endif
		ranks(i) = i;
	}

	// order the columns in the adjmat according to their distance to the
	// first column in the original serialization
	sort(ranks.begin(),ranks.end(),PositionSorter(pos));
	// order the positions ascendingly (-->matching the order of the ranks)
	sort(pos.begin(),pos.end());

	if(mNormalize){
		// min is 0 already, all positive --> just div by max
		pos /= *max_element(pos.begin(),pos.end());
	}

	Serialization ret(n);
	ret.setPositions(pos);
	ret.setRanks(ranks);
	return ret;
}

void GDistSeriationGen::configure()
{
	mSeriationGenName  = gCfg().getString("gdist-ser-gen.seriation_gen");
	mNormalize         = gCfg().getBool("gdist-ser-gen.normalize");
}

namespace{ registerInFactory<SerGenAdj, GDistSeriationGen> registerBase("GDistSeriationGen"); }
