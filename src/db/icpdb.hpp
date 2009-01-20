#include <ProbData.hpp>
#include <icp.hpp>
#include "db.hpp"
#include "DBPoint.hpp"
#ifndef __ICPDB_HPP__
#define __ICPDB_HPP__

#define DIM 3
class ICPDB : public GraphDB
{
  public:
	  typedef double                                    TPrecision;
	  typedef DBPoint<TPrecision,DIM>                   icp_point_type;
	  typedef util::ICP<DIM,icp_point_type,TPrecision>  TICP;
	  typedef TICP                                      value_type;
	  typedef std::vector<TICP>                         TICPVec;
  public:
    ICPDB();
	virtual void configure();
    virtual ~ICPDB();

	virtual void findNearest(const TCloud& cloud);
	virtual void add(const ProbAdjPerm& pap, const TCloud& cloud);
	virtual void init(int dim);


  private:
	TICPVec mDB;
	int               mICPMaxIters;
	float             mICPLambda;

};

#endif /* #ifndef __ICPDB_HPP__ */
