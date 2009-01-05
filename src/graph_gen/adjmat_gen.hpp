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
		virtual void configure();          ///< call setters/getters using data from config-file
		virtual ProbAdjPerm operator()();  ///< generate a adjacency matrix
		virtual ~AdjMatGen();              ///< virtual destructor
		virtual bool hasNext();            ///< whether there are more adjacency-matrices available.

		///< which class the graph belongs to
		virtual int getClassID(const std::string&ref="");          

		/// get a description of the node content in the form of a prolog program
		virtual std::string getPrologDescription(int ser_idx, const Serialization& s, const std::string& ref="");

		/// get a description of the node content suitable for one-line-per serialization presentation
		virtual std::string getPlainDescription(int ser_idx, const Serialization& s, const std::string& ref="");

		/// get a unique id for the graph
		virtual std::string getGraphID(const std::string&ref="");

		/// get adjacency matrix of graph
		virtual boost::shared_ptr<AdjMat::AdjMatT> getAdjMat(const std::string&ref="");

		/// for drawing: return additional node attributes in graphviz-syntax.
		/// should start with a comma if none-empty!
		virtual std::string getGraphVizNodeAttribs(int idx,const std::string&ref=""); 

		/// whether to be verbose
		static inline void setVerbose(bool b)  {mVerbose=b;}
		static inline bool  isVerbose()        {return mVerbose;}
	private:
		static bool mVerbose;
};

#endif /* __ADJMAT_GEN_HPP__ */

