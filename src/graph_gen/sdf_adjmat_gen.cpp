#include <sstream>
#include <fstream>
#include <numeric>
#include <configuration.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/graphviz.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <factory/factory.h>
#include <matlab_io.hpp>
#include "sdf_adjmat_gen.hpp"
//#include "jumbled_adjmat_gen.hpp"
#include <nana.h>

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;
namespace fs = boost::filesystem;
namespace ll = boost::lambda;

bool nextLineMatches(istream& is, string re){
	string s;
	getline(is, s);
	regex e(re);
	bool b = regex_match(s, e);
	if(!b){
		cerr << "Warning: Line `"<<s<<"' does not match expected pattern `"<<re<<"'"<<endl;
		return false;
	}
	return true;
}

void SDFAdjmatGen::configure()
{
	AdjMatGen::configure(); 
	string files = gCfg().getString("SDFAdjmatGen.files");
	mFixedSize   = gCfg().getInt("SDFAdjmatGen.fixed_size");
	setInputFiles(files);
}

SDFAdjmatGen::SDFAdjmatGen()
	:mFixedSize(0)
{
	// randomly between 0.5 and 1
	mFeatureWeights = ublas::vector<double>(getNumFeatures());
	generate(mFeatureWeights.begin(), mFeatureWeights.end(), drand48);
	transform(mFeatureWeights.begin(), mFeatureWeights.end(), mFeatureWeights.begin(), 0.5+ll::_1/2.0);
}
bool SDFAdjmatGen::hasNext()
{
	return mOutputCounter < (int)(0.90*mDescriptors.size());
}
SDFAdjmatGen::~SDFAdjmatGen()
{
}
std::string SDFAdjmatGen::getGraphID(const boost::any& ref)
{
	return any_cast<Descriptor*>(ref)->name;
}
std::string SDFAdjmatGen::getGraphVizNodeAttribs(int idx, const boost::any& ref)
{
	return string("");
}
std::string SDFAdjmatGen::getPlainDescription(int ser_idx, const Serialization& ser, const boost::any& ref)
{
	stringstream str;
	return str.str();
}
int SDFAdjmatGen::getClassID(const boost::any& ref)
{
	return any_cast<Descriptor*>(ref)->classid;
}

inline double my_exp(const double& d){ return exp(d); }

ProbAdjPerm SDFAdjmatGen::operator()()
{
	if(mDescriptors.size()==0){
		if(fs::exists("/tmp/graphs.ser")){
			ifstream ifs("/tmp/graphs.ser");
			archive::binary_iarchive ar(ifs);
			ar >> (*this);
		}else
			readInputFiles();
		mOutputCounter = 0;
	}
	if((unsigned int)mOutputCounter>=mDescriptors.size()){
		ProbAdjPerm pap;
		pap.setAdjMat(boost::shared_ptr<AdjMat::AdjMatT>());
		return pap;
	}

	Descriptor& desc  = mDescriptors[mOutputCounter++];
	string&     gname = desc.name;
	SDFGraph& graph   = desc.graph;
	SDFVertexIterator vi, vi_end;
	tie(vi, vi_end) = bgl::vertices(graph);
	int n = std::distance(vi, vi_end);

	shared_ptr<AdjMat::AdjMatT> adj(new AdjMat::AdjMatT(n, n));
	*adj = ublas::scalar_matrix<double>(n,n,0);
	SDFEdgeIterator ei, ei_end;

	for(tie(ei,ei_end)=bgl::edges(graph);ei!=ei_end;ei++){
		feature_t f = ublas::element_prod(edge_features(*ei, graph), mFeatureWeights);
		double d = 
			//1.0-exp(-accumulate(f.begin(), f.end(), 0.0, ll::_1 + ll::_2));
			accumulate(f.begin(), f.end(), 0.0, ll::_1 + ll::_2);
		d = max(d, 1E-5); // required by embedders
		bgl::put(bgl::edge_weight, graph, *ei, d);
		(*adj)(bgl::source(*ei, graph),bgl::target(*ei,graph)) = 
		(*adj)(bgl::target(*ei, graph),bgl::source(*ei,graph)) = 
		d;
	}

	ProbAdjPerm pap;
	pap.setAdjMat(adj);
	pap.setId(gname);
	pap.setBackground(&desc);

	return pap;
}

void SDFAdjmatGen::setInputFiles(const std::string& inputfileline)
{
	string line = boost::trim_copy(inputfileline);

	vector<string> file_class_pairs;
	boost::split( file_class_pairs, line, boost::is_any_of(",") );

	BOOST_FOREACH(const string& pair, file_class_pairs){
		vector<string> tmp;
		boost::split( tmp, pair, boost::is_any_of(":") );
		if(!tmp.size()==2)
		{
			cerr << "Warning: SDFAdjmatGen: InputFormatError: "<<pair<< " in " << line<< " should have format filename:classid"<<endl;
			continue;
		}
		FileDescriptor fd;
		fd.name = tmp[0];
		fd.classid = boost::lexical_cast<int>(tmp[1]);
		mInputFiles.push_back(fd);
	}
}

void SDFAdjmatGen::readInputFiles()
{
	int graphCount = 0;
	BOOST_FOREACH(const FileDescriptor& fd, mInputFiles){
		L("Reading SDF input file `%s'...", fd.name.c_str());
		ifstream is(fd.name.c_str());
		int fcnt = 0;
		while(readMolekule(is, fd.classid, graphCount++) && fcnt++<200);
		L(". %d molecules read.\n", graphCount);
	}
	ofstream os("/tmp/graphs.ser");
	archive::binary_oarchive ar(os);
	random_shuffle(mDescriptors.begin(), mDescriptors.end());
	ar << const_cast<const SDFAdjmatGen&>(*this);
}

bool SDFAdjmatGen::readMolekule(std::istream& is, int klass, int graphCount)
{
	Descriptor desc;
	desc.classid = klass;

	string line, id;
	if(!nextLineMatches(is, "\\s*"))            return false;
	if(!nextLineMatches(is, ".*OpenBabel.*"))   return false;
	if(!nextLineMatches(is, "\\s*"))            return false;
	if(!nextLineMatches(is, ".*V3000"))         return false;
	if(!nextLineMatches(is, ".*BEGIN	CTAB")) return false;
	getline(is, line); // M V30 COUNTS numberAtoms numberBonds 0 0 0
	boost::trim(line);
	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of("	") );
	if(!strvec.size()==8){
		cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: input description line `"<<line<<"' has unexpected format."<<endl;
		return false;
	}
	if(!nextLineMatches(is, ".*BEGIN	ATOM")) return false;

	// create a matrix with numberAtoms size
	int numberAtoms = boost::lexical_cast<int>(strvec[3]);
	int numberBonds = boost::lexical_cast<int>(strvec[4]);
	int n = (mFixedSize==0) ? numberAtoms : mFixedSize;

	desc.graph      = SDFGraph(n);
	SDFGraph& graph = desc.graph;

	// read all atoms
	SDFVertexIterator vi, vi_end;
	int atomCnt=0;
	for(tie(vi,vi_end) = bgl::vertices(graph); vi!=vi_end; ++vi, ++atomCnt){
		if(atomCnt == numberAtoms) break;
		getline(is, line);
		boost::trim(line);
		boost::split( strvec, line, boost::is_any_of("	") );
		if(strvec.size() < 8){
			cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: atom "<<graphCount<<" description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		string atomType = strvec[3]; // C or N or...
		bgl::get(bgl::vertex_name_t(),graph)[*vi] = atomType;
	}
	if(!nextLineMatches(is,".*END	ATOM"))return false;
	if(!nextLineMatches(is,".*BEGIN	BOND"))return false;
	
	// read all bonds
	for(int i=0;i<numberBonds; i++){
		getline(is, line);
		boost::trim(line);
		boost::split( strvec, line, boost::is_any_of("	") );
		if(strvec.size() != 6){
			cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: bond description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		int atomFrom = boost::lexical_cast<int>(strvec[4])-1;
		int atomTo   = boost::lexical_cast<int>(strvec[5])-1;
		int bondType = boost::lexical_cast<int>(strvec[3]);
		if(atomFrom >= numberAtoms || atomTo >= numberAtoms){
			cerr << "Warning: SDFAdjmatGen: InputFileFormat Error: bond description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		bgl::add_edge(atomFrom, atomTo, graph);
		bgl::get(bgl::get(bgl::bondnum, graph),bgl::edge(atomFrom, atomTo, graph).first) = bondType;
	}
	if(!nextLineMatches(is,".*END	BOND"))return false;
	if(!nextLineMatches(is,".*END	CTAB"))return false;
	if(!nextLineMatches(is,".*END"))return false;
	if(!nextLineMatches(is,"\\$\\$\\$\\$"))return false;

	stringstream sname;
	sname << "SDF_class"<<klass<<"_"<<graphCount;
	string gname = sname.str();
	desc.name = gname;
	bgl::get_property(graph, bgl::graph_name) = gname;
	bgl::get_property(graph, bgl::graph_classid) = klass;

	// debug output
	/*
	stringstream filename;
	filename << "../../data/hiv/results/"<<gname<<".dot";
	ofstream os(filename.str().c_str());
	if(!os.is_open())
		throw runtime_error(string("could not open for writing: ")+sname.str());

	std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
	vertex_attr["label"] = gname;

	bgl::write_graphviz(os, graph 
			,bgl::make_label_writer(bgl::get(bgl::vertex_name,graph))
			,bgl::make_label_writer(bgl::get(bgl::edge_weight,graph))
			//,bgl::make_graph_attributes_writer(graph)
			);
			*/
	
	mDescriptors.push_back(desc);

	return true;
}
void SDFAdjmatGen::rewind()
{
	for(unsigned int i=0;i<mDescriptors.size(); i++){
		Descriptor& d = mDescriptors[i];
		string s = (boost::format("../../data/hiv/results/SDF_class%d_%d.dot")%d.classid%i).str();
		ofstream os(s.c_str());
		bgl::write_graphviz(os, d.graph 
				,bgl::make_label_writer(bgl::get(bgl::vertex_name,d.graph))
				,bgl::make_label_writer(bgl::get(bgl::edge_weight,d.graph))
				);
	}
	mOutputCounter = 0;
	mDescriptors.clear();
}

void SDFAdjmatGen::setFeatureWeights(const boost::numeric::ublas::vector<double>& v)
{
	mFeatureWeights = v;
}


unsigned int SDFAdjmatGen::getNumFeatures() {
	return 17;
}
SDFAdjmatGen::feature_t
SDFAdjmatGen::edge_features(const SDFEdge& e, const SDFGraph& graph)
{
	feature_t f(getNumFeatures(),0);
	SDFVertex src = bgl::source(e, graph);
	SDFVertex tar = bgl::target(e, graph);
	int c=0;
	f(c++) = 
		bgl::get(bgl::bondnum,graph,e);
	f(c++) = 
		bgl::get(bgl::bondnum,graph,e) == 1;
	f(c++) = 
		bgl::get(bgl::bondnum,graph,e) == 2;
	f(c++) = 
		bgl::get(bgl::bondnum,graph,e) == 7;
	f(c++) = 
		bgl::get(bgl::vertex_name,graph,src) ==
		bgl::get(bgl::vertex_name,graph,tar);
	f(c++) = 
		bgl::get(bgl::vertex_name,graph,src) == "C" ||
		bgl::get(bgl::vertex_name,graph,tar) == "C"   ;
	f(c++) = 
		bgl::get(bgl::vertex_name,graph,src) == "C" &&
		bgl::get(bgl::vertex_name,graph,tar) == "C"   ;
	f(c++) = 
		(bgl::get(bgl::vertex_name,graph,src) == "O")  +
		(bgl::get(bgl::vertex_name,graph,tar) == "O")   ;
	f(c++) = 
		(bgl::get(bgl::vertex_name,graph,src) == "S")  +
		(bgl::get(bgl::vertex_name,graph,tar) == "S")   ;
	f(c++) = 
		(bgl::get(bgl::vertex_name,graph,src) == "N")  +
		(bgl::get(bgl::vertex_name,graph,tar) == "N")   ;
	f(c++) = 
		(bgl::get(bgl::vertex_name,graph,src) == "Cl")  +
		(bgl::get(bgl::vertex_name,graph,tar) == "Cl")   ;
	f(c++) = 
		(bgl::get(bgl::vertex_name,graph,src) == "P")  +
		(bgl::get(bgl::vertex_name,graph,tar) == "P")   ;

	int snap=0;
	// 2nd order features
	SDFVertex srcdest[2] = {src,tar};
	SDFAdjIterator it, it_end;
	for(int i = 0; i< 2; i++){
		for(tie(it,it_end) = adjacent_vertices(srcdest[i],graph);
				it!=it_end; it++){
			snap = c;
			f(snap) += bgl::get(bgl::vertex_name,graph,*it) == "N";  snap++;
			f(snap) += bgl::get(bgl::vertex_name,graph,*it) == "C";  snap++;
			f(snap) += bgl::get(bgl::vertex_name,graph,*it) == "S";  snap++;
			f(snap) += bgl::get(bgl::vertex_name,graph,*it) == "P";  snap++;
			f(snap) += bgl::get(bgl::vertex_name,graph,*it) == "Cl"; snap++;
		}
	}

	for(unsigned int i=0;i<f.size();++i)
		if(f(i)==0) f(i) = -1;
	
	return f;
}

SDFAdjmatGen::feature_t 
SDFAdjmatGen::getFeatures(int idx, const boost::any& ref)
{
	SDFAdjmatGen::feature_t f(getNumFeatures(),0);
	Descriptor& desc = * any_cast<Descriptor*>(ref);
	SDFGraph& graph   = desc.graph;
	SDFVertex v = bgl::vertex(idx,graph);

	SDFOutEdgeIterator it, it_end;
	tie(it,it_end) = bgl::out_edges(v, graph);
	int cnt=0;
	for(;it!=it_end;++it, ++cnt){
		SDFEdge   e   = *it;
		f += edge_features(e, graph);
	}
	if(cnt>0)
		f/=cnt;
	
	return f;
}

const boost::numeric::ublas::vector<double>& 
SDFAdjmatGen::getFeatureWeights()const
{
	return mFeatureWeights;
}

namespace{ registerInFactory<AdjMatGen, SDFAdjmatGen> registerBase("SDFAdjmatGen"); }
