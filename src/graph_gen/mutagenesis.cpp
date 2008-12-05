// 2008-11-17 Hannes Schulz <mail at hannes-schulz dot de>

#include	<string>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<list>
#include	<factory/factory.h>
#include	<boost/algorithm/string.hpp>
#include	<boost/lexical_cast.hpp>
#include    <boost/numeric/ublas/matrix.hpp>
#include	<matlab_io.hpp>

#include	<configuration.hpp>
#include	"mutagenesis.hpp"
#include	<nana.h>
#undef A

using namespace std;
namespace ublas = boost::numeric::ublas;

/**********************************************************
 *          Mutagenesis Implementation                *
 **********************************************************/
struct Mutagenesis::Impl{
	string mInputFilename;
	ifstream mInputStream;

	Impl();
	~Impl();

	ProbAdjPerm nextGraph();
	void open();
	bool hasNext();
	string getPrologDescription(int idx,const Serialization&);
	string getPlainDescription(int idx,const Serialization&);

	vector<string> mNames;
	vector<string> mTypes;
	vector<string> mParam1;
	vector<string> mParam2;
	string         mName;
};

string Mutagenesis::Impl::getPlainDescription(int ser_idx,const Serialization&s){
	return mTypes[s[ser_idx]];
}
string Mutagenesis::Impl::getPrologDescription(int ser_idx,const Serialization&s){
	unsigned int idx = s[ser_idx];
	if(idx>=mTypes.size()){
		L("idx = %d, size = %d\n",idx,mTypes.size()); 
		throw runtime_error("Mutagenesis::getPrologDescription(): Index too large");
	}
	ostringstream o;
	if(mTypes[idx] == "b"){
		// bond
		o << "chemBond("<<mParam1[idx] << ")";
	}else{
		// atom
		//o << "chemAtom(" << mName <<","<< mNames[idx] << "," << mTypes[idx] << ")";
		o << "chemAtom(" << mTypes[idx] << ")";
		o << ";weight(" << mParam1[idx] << ")";
	}
	for(int i=-4;i<=4;i++){
		int nidx = idx+i;              // nidx: position of neighbour in serialization
		if(nidx<0)         continue;
		if(nidx>=s.size()) continue;
		int idx = s[nidx];             // idx: original position of neighbour
		if(mTypes[idx]!="b"){
			o << ";neighChemAtom("<<mTypes[idx]<<")";
		}
	}
	return o.str();
}

bool Mutagenesis::Impl::hasNext(){
	return mInputStream.is_open() && mInputStream.good();
}
void Mutagenesis::Impl::open(){
	mInputStream.open(mInputFilename.c_str());
}

ProbAdjPerm Mutagenesis::Impl::nextGraph(){
	string line("      ");
	int cnt=0;
	vector<string> names;
	vector<string> types;
	vector<string> param1;
	vector<string> param2;
	int n = -1;
	boost::shared_ptr<AdjMat::AdjMatT> A_ptr;

	getline(mInputStream, line);
	mName = line;

	while(line.length()>1 && getline(mInputStream, line)){
		list<string> strvec;
		boost::split( strvec, line, boost::is_any_of(",") );
		int addParams=4;
		if(cnt==0){
			for(int i=0;i<addParams;i++) strvec.pop_front(); // ignore 1st 4 "None"s
			list<string>::iterator it = strvec.begin();
			copy(it,strvec.end(),back_inserter(names));
			n = names.size();
			A_ptr.reset(new ublas::matrix<double>(n,n));
			cnt++;
			continue;
		}
		AdjMat::AdjMatT& A = *A_ptr;
		if(A.size1() != strvec.size()-addParams){
			if(cnt-1 == n)
				break;
			L("Breaking at cnt=%d n=%d A.size1()=%d, strvecsize=%d\n",cnt,n,A.size1(),strvec.size());
			throw runtime_error("Mutagenesis::Impl::nextGraph() -- input file format error");
		}
		list<string>::iterator it = strvec.begin();
		string name = *it;
		it++;
		types.push_back( *it ); it++;
		param1.push_back( *it ); it++;
		param2.push_back( *it ); it++;
		for(int i=0;it!=strvec.end();it++,i++){
			int val = boost::lexical_cast<double>(*it);
			A(i,cnt-1) = A(cnt-1,i) = val>0 ? 1:0;
		}
		cnt++;
	}
	ProbAdjPerm prob;
	prob.setAdjMat(A_ptr);
	mNames = names;
	mTypes = types;
	mParam1 = param1;
	mParam2 = param2;
	return prob;
}

// Impl constructor
Mutagenesis::Impl::Impl(){
}

// Impl destructor
Mutagenesis::Impl::~Impl(){
}


/**********************************************************
 *          Mutagenesis Interface                     *
 **********************************************************/
Mutagenesis::Mutagenesis()
	:mImpl(new Impl())
{
}

Mutagenesis::~Mutagenesis()
{
  // cleanup
}
void Mutagenesis::configure()
{
	mImpl->mInputFilename = gCfg().getString("mutagenesis.in-file");
	mImpl->open();
}
ProbAdjPerm Mutagenesis::operator()()
{
	return mImpl->nextGraph();
}
bool Mutagenesis::hasNext()
{
	return mImpl->hasNext();
}
void Mutagenesis::open()
{
	mImpl->open();
}
string Mutagenesis::getPlainDescription(int ser_idx, const Serialization& s)
{
	return mImpl->getPlainDescription(ser_idx,s);
}
string Mutagenesis::getPrologDescription(int ser_idx,const Serialization& s)
{
	return mImpl->getPrologDescription(ser_idx,s);
}
std::string Mutagenesis::getGraphID()
{
	return mImpl->mName;
}


namespace{ registerInFactory<AdjMatGen, Mutagenesis> registerBase("Mutagenesis"); }
