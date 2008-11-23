#ifndef __PLAINPRINT_HPP__
#define __PLAINPRINT_HPP__

#include <iostream>
#include <fstream>

#include "postproc.hpp"

class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief PlainPrint print seriation in one line each
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-21
 */

class PlainPrint : public PostProc
{
  public:

    /**
     * Default constructor
     */
    PlainPrint();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atEnd();


    /**
     * Destructor
     */
    virtual ~PlainPrint();

  private:
	std::ofstream mOS;
	void seqhead(AdjMatGen& gen);
	void seqfoot();

};

#endif /* #ifndef __PLAINPRINT_HPP__ */

