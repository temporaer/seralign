#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

/**
 * @brief Serialization A Serialization
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */
#include <vector>
#include <boost/numeric/ublas/vector.hpp>

class Serialization
{
  public:
	  typedef boost::numeric::ublas::vector<int>   RankT;
	  typedef boost::numeric::ublas::vector<float> PosT;
	  typedef std::vector<PosT>                    PosVecT;
  private:
	  RankT   mRanks;
	  PosT    mPositions;
	  PosVecT mPosVecs;
	  bool    mHasPosVecs;
  public:

	  Serialization(int n):mRanks(n),mPositions(n),mPosVecs(n),mHasPosVecs(false){}
	  Serialization(const RankT&r, const PosT&p):mRanks(r),mPositions(p),mPosVecs(r.size()),mHasPosVecs(false){}
	  Serialization(const RankT&r):mRanks(r),mPositions(r.size()),mPosVecs(r.size()),mHasPosVecs(false){}

	  inline void setRanks    (const RankT& r){ mRanks = r; }
	  inline void setPositions(const PosT& p) { mPositions = p; }
	  
	  template<class T> 
	  inline void setPosition(int i, const boost::numeric::ublas::vector_expression<T>& v){
		  mHasPosVecs = true;
		  mPosVecs[i] = v;
	  }

	  inline const RankT&   getRanks    ()const{ return mRanks; }
	  inline       RankT&   getRanks    (){ return mRanks; }
	  inline const PosT&    getPositions()const{ return mPositions; }
	  inline const PosVecT& getPosVecs  ()const{ return mPosVecs;   }

	  virtual ~Serialization();

};

#endif /* #ifndef __SERIALIZATION_HPP__ */
