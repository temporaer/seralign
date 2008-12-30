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
    GraphFromAdj(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2);
    GraphFromAdj(const ProbAdjPerm&);
    virtual ~GraphFromAdj();

	double getDist(int i);
	double getTotalDist();
	double getDistFromIdx1(int i);
	double getDistFromIdx2(int i);
	void   setIdx1(int i);
	void   setIdx2(int i);
	unsigned int    getFarthestFromIdx1();
	unsigned int    getFarthestFromIdx2();

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __GRAPHFROMADJ_HPP__ */
