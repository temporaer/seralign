#ifndef __RANDOM_PAT_ADJMAT_HPP__
#define __RANDOM_PAT_ADJMAT_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <adjmat_gen.hpp>
#include <vector>

class RandomPatAdjmatGen: public AdjMatGen{
	private:
		int     mSize;             //< Size of AdjMatrix to be generated
		int     mPatSize;          //< Size of Pattern hidden in Random matrix
		float   mEdgeInclusionProb;//< probability of an edge to be included
		int     mRunningID;        //< used to construct unique ID for a graph
		bool    mGeneratePattern;  //< if false, randomly assign nodes the part_of_pattern property
		struct RNode{
			bool part_of_pattern;
			float x;
			float y;
		};
		struct Edge : public std::pair<int, int>{ };
		std::vector<RNode> mNodes;
	public:
		virtual void configure();
		RandomPatAdjmatGen();
		virtual ProbAdjPerm operator()();
		virtual ~RandomPatAdjmatGen();
		virtual bool hasNext();
		virtual std::string getPrologDescription(int ser_idx, const Serialization&, const std::string&);
		virtual std::string getPlainDescription(int ser_idx, const Serialization&, const std::string&);
		virtual std::string getGraphID(const std::string&);
		virtual std::string getGraphVizNodeAttribs(int idx, const std::string&); //< should start with a comma!

		inline void  setSize(int i){mSize=i;}
		inline int   getSize()     {return mSize;}
		inline void  setPatSize(int i){mPatSize=i;}
		inline int   getPatSize()     {return mPatSize;}
		inline void  setEdgeInclusionProb(float f){mEdgeInclusionProb=f;}
		inline float getEdgeInclusionProb()       {return mEdgeInclusionProb;}
		inline void  setGeneratePattern(bool b){mGeneratePattern=b;}
		inline bool  getGeneratePattern()       {return mGeneratePattern;}

	private:
		void fillNodeArrayCirclePattern();
		void triangulateNodes(std::vector<Edge>& edges);
};

#endif /* __RANDOM_PAT_ADJMAT_HPP__ */

