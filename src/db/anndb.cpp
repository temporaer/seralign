#include <ANN/ANN.h>
#include <map>
#include <fstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/foreach.hpp>
#include <configuration.hpp>
#include "anndb.hpp"
#include <nana.h>

using namespace std;
namespace ublas = boost::numeric::ublas;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

ANNDB::ANNDB()
{
}

void ANNDB::configure()
{
	mDim        = gCfg().getInt("BuildDB.embed_dim");
	mBucketSize = gCfg().getInt("ANNDB.bucket_size");
}

void ANNDB::init(int dim)
{
	mDim=dim;
}

void ANNDB::finish()
{
	int ptnum = 0;
	BOOST_FOREACH(const TCloud& c, mClouds){
		ptnum += c.size();
	}
	ANNpointArray dataPts;
	dataPts = annAllocPts(ptnum, mDim);

	int idx=0;
	BOOST_FOREACH(const TCloud& c, mClouds){
		BOOST_FOREACH(const point_type& d, c){
			for(int i=0;i<mDim;i++)
				dataPts[idx][i] = d[i];
			idx++;
		}
	}

	mAllTree.reset(new ANNkd_tree( // build search structure
      dataPts,         
      ptnum, 
      mDim,      // dimension of space
	  mBucketSize
	  ));        
	ofstream os("/tmp/tree.dmp");
	mAllTree->Dump(ANNtrue,os);

	char str[255];
	mLeafCount = 0;
	BOOST_FOREACH(const TCloud& c, mClouds){
		BOOST_FOREACH(const point_type& d, c){
			memset(&str[0],0,255);
			mAllTree->annTrace(const_cast<double*>(&d[0]), &str[0]);
			for(int c=254;c>=0;c--){
				if(str[c]==0) continue;
				if(mLeafs.find(string(&str[0])) == mLeafs.end())
					mLeafs[string(&str[0])] = mLeafCount++;
				str[c]='\0';
			}
		}
	}
	cout << endl<<"Number of leafs: " << mLeafCount<<endl;
	//mClouds.clear();
}

boost::numeric::ublas::vector<int>    
ANNDB::getFeatures(int i)
{
	I(i<(int)mClouds.size());
	return getFeatures(mClouds[i]);
}

const ANNDB::TCloud&
ANNDB::getCloud(int i)
{
	I(i<(int)mClouds.size());
	return mClouds[i];
}


boost::numeric::ublas::vector<int>    
ANNDB::getFeatures(const TCloud& cloud)
{
	char str[255];
	ublas::vector<int> res(mLeafCount,0);
	BOOST_FOREACH(const point_type& p, cloud){
		memset(&str[0],0,255);
		mAllTree->annTrace(const_cast<double*>(&p[0]), &str[0]);
		for(int c=254;c>=0;c--){
			if(str[c]==0) continue;
			res(mLeafs[string(&str[0])]) ++;
			//res(mLeafs[string(&str[0])]) =1;
			str[c]='\0';
		}
	}
	return res;
}


void ANNDB::add(const ProbAdjPerm& pap, const TCloud& cloud)
{
	mPaps.push_back(pap);
	mClouds.push_back(cloud);
}


ANNDB::~ANNDB()
{
  // cleanup
}
