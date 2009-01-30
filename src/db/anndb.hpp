#ifndef __ANNDB_HPP__
#define __ANNDB_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <memory>
#include <map>
#include "db.hpp"

class ANNkd_tree;
class ANNDB : public GraphDB 
{
  public:
    ANNDB();
	virtual void configure();
	virtual void add(const ProbAdjPerm& pap, const TCloud& cloud);
	virtual boost::numeric::ublas::vector<int>    getFeatures(const TCloud& cloud);
	virtual boost::numeric::ublas::vector<int>    getFeatures(int i);
	virtual const TCloud&                         getCloud(int i);
	virtual void finish();
	virtual void init(int dim);
    virtual ~ANNDB();

  private:
	std::auto_ptr<ANNkd_tree> mAllTree;
	std::vector<TCloud>       mClouds;
	std::map<std::string,int>      mLeafs;
	int mBucketSize;
	int mLeafCount;
	int mDim;
};

#endif /* #ifndef __ANNDB_HPP__ */
