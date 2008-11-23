#ifndef __POSTPROC_HPP__
#define __POSTPROC_HPP__
#include <Serialization.hpp>
class AdjMatGen;
class ProbAdjPerm;

/**
 * @brief PostProc post-processing of seriations
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-21
 */

class PostProc
{
  public:

    /**
     * Default constructor
     */
    PostProc();

    /**
     * Destructor
     */
    virtual ~PostProc();

	virtual void atStart();
	virtual void atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob);
	virtual void atEnd();

  private:

};

#endif /* #ifndef __POSTPROC_HPP__ */
