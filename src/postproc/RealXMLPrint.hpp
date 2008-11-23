#ifndef __REALXMLPRINT_HPP__
#define __REALXMLPRINT_HPP__

#include <iostream>
#include <fstream>

#include "postproc.hpp"

class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief RealXMLPrint print seriation to REAL xml
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-21
 */

class RealXMLPrint : public PostProc
{
  public:

    /**
     * Default constructor
     */
    RealXMLPrint();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atEnd();


    /**
     * Destructor
     */
    virtual ~RealXMLPrint();

  private:
	std::ofstream mOS;
	void seqhead(AdjMatGen& gen);
	void seqfoot();

};

#endif /* #ifndef __REALXMLPRINT_HPP__ */
