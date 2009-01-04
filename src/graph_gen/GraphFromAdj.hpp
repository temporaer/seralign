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

	double getProjection(unsigned int k, int i);
	double getDist(unsigned int k, int i, int j);
	double getTotalDist(unsigned int k);
	//double getDistFromA(unsigned int k, int i);
	//double getDistFromB(unsigned int k, int i);
	void   setA(unsigned int k, int i);
	void   setB(unsigned int k, int i);
	unsigned int    getFarthestFromA(unsigned int k);
	unsigned int    getFarthestFromB(unsigned int k);

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __GRAPHFROMADJ_HPP__ */
