#ifndef __POINTNN_HPP__
#define __POINTNN_HPP__
#include <ProbData.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "db.hpp"

class AdjMatGen;
class PointNNDB : public GraphDB
{
  public:
    PointNNDB();
	virtual void configure();
    virtual ~PointNNDB();

	virtual void findNearest(const TCloud& cloud);
	virtual void evaluate(AdjMatGen& gen, int k);
	virtual void add(const ProbAdjPerm& pap, const TCloud& cloud);
	virtual void init(int dim);
	virtual void finish();
	virtual void clear();
	virtual int classify(const TCloud& cloud, AdjMatGen& gen, int k);

  private:
	struct Impl;
	boost::shared_ptr<Impl> mImpl;
	
};

#endif /* #ifndef __POINTNN_HPP__ */
