
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include "adjmat_gen.hpp"


class RandomAdjMatGen : public AdjMatGen{
	private:
		int   mMatrixSize;
		float mConnectionProb;
		bool  mWeightedEdges;
		float mSeed;
	public:
		virtual void configure();
		inline void setMatrixSize(int n){mMatrixSize=n;}
		inline void setConnectionProb(float f){mConnectionProb=f;}
		inline void setWeightedEdges(bool b){mWeightedEdges=b;}
		inline void setSeed(float s){mSeed=s;}

		virtual boost::shared_ptr<AdjMatT> operator()();
		virtual ~RandomAdjMatGen();
};

