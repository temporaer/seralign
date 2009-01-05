#ifndef __FULL_CONN_ADJMAT_HPP__
#define __FULL_CONN_ADJMAT_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adjmat_gen.hpp>
#include <vector>

class FullConnAdjmatGen: public AdjMatGen{
	private:
		int     mSize;             //< Size of AdjMatrix to be generated
		int     mDim;
		int     mRunningID;        //< used to construct unique ID for a graph
		boost::numeric::ublas::matrix<double> mAdj;
		std::vector<boost::numeric::ublas::vector<double> > mVertexPos;
	public:
		virtual void configure();
		FullConnAdjmatGen();
		virtual ProbAdjPerm operator()();
		virtual ~FullConnAdjmatGen();
		virtual bool hasNext();
		virtual std::string getPlainDescription(int ser_idx, const Serialization&, const std::string&);
		virtual std::string getGraphID(const std::string&);
		virtual std::string getGraphVizNodeAttribs(int idx, const std::string&); //< should start with a comma!
		virtual int getClassID(const std::string&);          ///< which class the graph belongs to

		inline void  setSize(int i){mSize=i;}
		inline int   getSize()     {return mSize;}
};

#endif /* __FULL_CONN_ADJMAT_HPP__ */

