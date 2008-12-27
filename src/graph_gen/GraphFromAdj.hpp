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
    virtual ~GraphFromAdj();

	double getDist(int i);

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __GRAPHFROMADJ_HPP__ */
