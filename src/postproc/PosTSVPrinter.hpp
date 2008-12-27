#ifndef __POSTSVPRINT_HPP__
#define __POSTSVPRINT_HPP__

#include <iostream>
#include <fstream>

#include "postproc.hpp"

class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief PosTSVPrint print seriation in one line each, but only _positions_
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-12-23
 */

class PosTSVPrint : public PostProc
{
  public:

    /**
     * Default constructor
     */
    PosTSVPrint();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atEnd();


    /**
     * Destructor
     */
    virtual ~PosTSVPrint();

  private:
	std::ofstream mOS;
	void seqhead(AdjMatGen& gen);
	void seqfoot();

};

#endif /* #ifndef __POSTSVPRINT_HPP__ */

