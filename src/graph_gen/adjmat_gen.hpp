#ifndef __ADJMAT_GEN_HPP__
#define __ADJMAT_GEN_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <Serialization.hpp>
#include <ProbData.hpp>



/**
 * Adjacency matrix generator base class.
 */
class AdjMatGen {
	public:
		typedef boost::numeric::ublas::vector<double> feature_t;
		virtual void configure();          ///< call setters/getters using data from config-file
		virtual ProbAdjPerm operator()();  ///< generate a adjacency matrix
		virtual ~AdjMatGen();              ///< virtual destructor
		virtual bool hasNext();            ///< whether there are more adjacency-matrices available.

		///< which class the graph belongs to
		virtual int getClassID(const boost::any&);          

		/// get a description of the node content in the form of a prolog program
		virtual std::string getPrologDescription(int ser_idx, const Serialization& s, const boost::any&);

		/// get a description of the node content suitable for one-line-per serialization presentation
		virtual std::string getPlainDescription(int ser_idx, const Serialization& s, const boost::any&);

		/// get a unique id for the graph
		virtual std::string getGraphID(const boost::any&);

		/// get adjacency matrix of graph
		virtual boost::shared_ptr<AdjMat::AdjMatT> getAdjMat(const boost::any&);

		/// for drawing: return additional node attributes in graphviz-syntax.
		/// should start with a comma if none-empty!
		virtual std::string getGraphVizNodeAttribs(int idx,const boost::any&); 

		/// get number of features of a node
		virtual unsigned int getNumFeatures(){ return 0;}

		/// get features of i-th node in problem referred to by ref
		virtual feature_t getFeatures(int idx, const boost::any&);

		/// set feature weights
		virtual void setFeatureWeights(const boost::numeric::ublas::vector<double>&);

		/// set feature weights
		virtual const boost::numeric::ublas::vector<double>& getFeatureWeights()const;

		/// start over with generating the adj-matrices
		virtual void rewind();

		/// whether to be verbose
		static inline void setVerbose(bool b)  {mVerbose=b;}
		static inline bool  isVerbose()        {return mVerbose;}
	private:
		static bool mVerbose;
};

#endif /* __ADJMAT_GEN_HPP__ */

