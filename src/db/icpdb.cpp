#include <boost/foreach.hpp>
#include <configuration.hpp>
#include "icpdb.hpp"

#include <nana.h>

using namespace std;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

ICPDB::ICPDB()
{
}

void ICPDB::findNearest(const TCloud& cloud)
{
	//TCloud c = cloud;
	//BOOST_FOREACH(TICP& icp, mDB){
		//icp.match(c.begin(), c.end()
			////,icp_point_type::Translator(),
			////,icp_point_type::Rotator<TICP::TQuat>());
			//);
	//}
}


void ICPDB::add(const ProbAdjPerm& pap, const TCloud& cloud)
{
	TICP icp(mICPMaxIters,mICPLambda);
	mPaps.push_back(pap);
	mDB.push_back(icp);
	std::vector<icp_point_type> v;
	int idx = 0;
	BOOST_FOREACH(const point_type& p, cloud){
		v.push_back(icp_point_type());
		v.back().pos = p;
		v.back().id  = idx++;
	}
	mDB.back().registerModel(v.begin(), v.end(), icp_point_type::Translator());
}


void ICPDB::init(int dim)
{
	I(DIM == dim);
}



ICPDB::~ICPDB()
{
  // cleanup
}

void ICPDB::configure()
{
	GraphDB::configure();
	mICPMaxIters      = gCfg().getInt("GDistProjectedDB.maxiter");
	mICPLambda        = gCfg().getFloat("GDistProjectedDB.lambda");
}


