#ifndef __GRAPHVIZPRINT_HPP__
#define __GRAPHVIZPRINT_HPP__

#include <iostream>
#include <vector>
#include <fstream>

#include "postproc.hpp"

class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief GraphVizPrint print seriation to a .dot file
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-27
 */

class GraphVizPrint : public PostProc
{
  public:

    /**
     * Default constructor
     */
    GraphVizPrint();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atSeriation(AdjMatGen& gen, std::vector<boost::numeric::ublas::vector<double> >& cloud, ProbAdjPerm& prob);
	virtual void atEnd();


    /**
     * Destructor
     */
    virtual ~GraphVizPrint();

  private:
	std::ofstream mOS;
	void seqhead(AdjMatGen& gen, int size);
	void seqfoot();

};

#endif /* #ifndef __GRAPHVIZPRINT_HPP__ */

