#ifndef __POSTPROC_HPP__
#define __POSTPROC_HPP__
#include <Serialization.hpp>
#include <vector>
#include <DBPoint.hpp>
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
	typedef std::vector<DBPoint<double,3> > DBCloud;

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
	virtual void atSeriation(AdjMatGen& gen, const DBCloud&, ProbAdjPerm& prob);
	virtual void atEnd();

  private:

};

#endif /* #ifndef __POSTPROC_HPP__ */
