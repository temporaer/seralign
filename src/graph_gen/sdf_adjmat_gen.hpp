#ifndef __SDF_ADJMAT_HPP__
#define __SDF_ADJMAT_HPP__

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/filesystem.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
//#include <boost/serialization/hash_set.hpp>
//#include <boost/serialization/hash_map.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include <adjmat_gen.hpp>
namespace bgl = boost;


namespace boost{
	enum edge_bondnum_t {bondnum=112};
	enum graph_classid_t {graph_classid=111};
	BOOST_INSTALL_PROPERTY(graph, classid);
	BOOST_INSTALL_PROPERTY(edge, bondnum);
};

class SDFAdjmatGen: public AdjMatGen{
	private:
		typedef bgl::adjacency_list<
			bgl::listS, // store out-edges of vertex in std::list
			bgl::vecS, // store vertex-set in std::vector
			bgl::undirectedS, // not directed
			bgl::property<bgl::vertex_name_t, std::string>,  // vertex-properties
			bgl::property<bgl::edge_weight_t, double,        // edge-properties
				bgl::property<bgl::edge_bondnum_t, int>
			>,       
			bgl::property<bgl::graph_name_t, std::string,    // graph-properties
			  bgl::property<bgl::graph_classid_t, int>
			>
		>    SDFGraph;
		typedef bgl::property_map<SDFGraph, bgl::graph_name_t>::type SDFGraphNameMap;
		typedef bgl::property_map<SDFGraph, bgl::graph_classid_t>::type SDFGraphClassIDMap;
		typedef bgl::property_map<SDFGraph, bgl::vertex_name_t>::type SDFVertexNameMap;
		typedef bgl::property_map<SDFGraph, bgl::edge_weight_t>::type SDFEdgeWeightMap;
		typedef bgl::graph_traits<SDFGraph>::vertex_descriptor SDFVertex;
		typedef bgl::graph_traits<SDFGraph>::edge_descriptor SDFEdge;
		typedef bgl::graph_traits<SDFGraph>::vertex_iterator SDFVertexIterator;
		typedef bgl::graph_traits<SDFGraph>::adjacency_iterator SDFAdjIterator;
		typedef bgl::graph_traits<SDFGraph>::edge_iterator SDFEdgeIterator;
		typedef bgl::graph_traits<SDFGraph>::out_edge_iterator SDFOutEdgeIterator;
		int          mFixedSize;
	public:
		virtual void configure();
		SDFAdjmatGen();
		virtual ProbAdjPerm operator()();
		virtual ~SDFAdjmatGen();
		virtual bool hasNext();
		virtual std::string getPlainDescription(int ser_idx, const Serialization&, const boost::any&);
		virtual std::string getGraphID(const boost::any&);
		virtual std::string getGraphVizNodeAttribs(int idx, const boost::any&); //< should start with a comma!
		virtual int getClassID(const boost::any&);          ///< which class the graph belongs to
		virtual void rewind();

		virtual unsigned int getNumFeatures();
		virtual feature_t getFeatures(int idx, const boost::any&);
		virtual const boost::numeric::ublas::vector<double>& getFeatureWeights()const;
		virtual void setFeatureWeights(const boost::numeric::ublas::vector<double>&);
		
	private:
		int mOutputCounter;

		struct FileDescriptor{
			std::string name;
			int         classid;
		};
		std::vector<FileDescriptor> mInputFiles; 
		boost::numeric::ublas::vector<double> mFeatureWeights;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& ar, const unsigned int version){
			ar & mDescriptors;
		}

		struct Descriptor{
			SDFGraph          graph;
			int               classid;
			std::string       name;
			template <class Archive>
				void serialize(Archive& ar, const unsigned int version){
					ar & graph;
					ar & classid;
					ar & name;
				}
		};
		std::vector<Descriptor>           mDescriptors;

		void setInputFiles(const std::string& s);
		void readInputFiles();
		bool readMolekule(std::istream& is,int klass, int cnt);
		feature_t edge_features(const SDFEdge& e, const SDFGraph&);
};

#endif /* __SDF_ADJMAT_HPP__ */

