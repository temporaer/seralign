#ifndef __PROBDATA_HPP__
#define __PROBDATA_HPP__

#include	"AdjMat.hpp"
#include    "Laplacian.hpp"
#include	"PermMat.hpp"
#include	"SimMat.hpp"
#include	<string>

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
	
	inline std::string getId()const{return mId;}
	inline void setId(const std::string& s){mId=s;}

  private:

	std::string mId;
};

class ProbAdjLapPerm : public ProbAdjPerm , public Laplacian
{
	public:
	void calculateLaplacian();
	ProbAdjLapPerm(const ProbAdjPerm&);
};

#endif /* #ifndef __PROBDATA_HPP__ */
