#ifndef __GRAPHFROMADJ_HPP__
#define __GRAPHFROMADJ_HPP__

#include <boost/shared_ptr.hpp>

/**
 * Make Boost Graph from adjmat
 */
class ProbAdjPerm;
class GraphFromAdj
{
  public:
    GraphFromAdj(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2, int dim=3);
    GraphFromAdj(const ProbAdjPerm&, int dim=3);
    virtual ~GraphFromAdj();

	double getDist(int refid, int i);
	double getTotalDist(int refid);
	double getDistFromA(int refid, int i);
	double getDistFromB(int refid, int i);
	void   setA(int refid, int i);
	void   setB(int refid, int i);
	unsigned int    getFarthestFromA(int refid);
	unsigned int    getFarthestFromB(int refid);

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __GRAPHFROMADJ_HPP__ */
