#ifndef __PROBDATA_HPP__
#define __PROBDATA_HPP__

#include	"AdjMat.hpp"
#include    "Laplacian.hpp"
#include	"PermMat.hpp"
#include	"SimMat.hpp"

class ProbAdjPerm
: public AdjMat
, public PermMat
{
  public:

    /**
     * Default constructor
     */
    ProbAdjPerm();

    /**
     * Destructor
     */
    virtual ~ProbAdjPerm();

  private:

};

#endif /* #ifndef __PROBDATA_HPP__ */
