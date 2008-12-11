// 2008-11-17 Hannes Schulz <mail at hannes-schulz dot de>

#include	<string>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<list>
#include    <map>
#include	<factory/factory.h>
#include	<boost/algorithm/string.hpp>
#include	<boost/lexical_cast.hpp>
#include    <boost/numeric/ublas/matrix.hpp>
#include    <boost/numeric/ublas/matrix_proxy.hpp>
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

	double kernelNull(int i, int j);
	double kernelTypeAndWeightEq(int i, int j);
	double kernelTypeSum(int i, int j);
	double kernelNienhuysCheng(int i, int j);

	double (Mutagenesis::Impl:: *mKernel)(int, int);
	
	boost::shared_ptr<AdjMat::AdjMatT> mA_ptr;
	vector<string> mNames;
	vector<string> mTypes;
	vector<string> mParam1;
	vector<string> mParam2;
	string         mName;
};

double Mutagenesis::Impl::kernelNienhuysCheng(int i, int j){
	double d = 0.0;
	int    n = 0;
	if(mTypes[i]  == mTypes[j])  d+=1; n++;
	if(mParam1[i] == mParam1[j]) d+=1; n++;
	return d/(2*n);
}
double Mutagenesis::Impl::kernelTypeSum(int i, int j){
	double d = 0.0;
	map<string, double> m;
	m["n"]  = 0.3;
	m["f"]  = 0.3;
	m["cl"] = 0.3;
	m["o"]  = 0.2;
	m["h"]  = 0.1;
	m["c"]  = 0.4;

	d += m[mTypes[i]] + m[mTypes[j]];

	return d;
}
double Mutagenesis::Impl::kernelTypeAndWeightEq(int i, int j){
	double d = 0.0;
	if(mTypes[i] == mTypes[j]) d+=0.25;
	if(mParam1[i] == mParam1[j]) d+=0.25;
	return d;
}
double Mutagenesis::Impl::kernelNull(int i, int j){
	return 0.0;
}
string Mutagenesis::Impl::getPlainDescription(int ser_idx,const Serialization&s){
	stringstream str;
	AdjMat::AdjMatT& A = *mA_ptr;
	int idx = s[ser_idx];
	str << mTypes[idx];
	bool incGraphNeigh = gCfg().getBool("mutagenesis.include-graph-neighbours");
	if(incGraphNeigh){
		str << "[";
		ublas::vector<double> col = ublas::column(A,idx);
		for(int i=0;i<A.size1();i++){
			if(i==idx)              continue;
			if(col(i) < 0.0001)     continue;
				str << mTypes[i];
		}
		str<< "]";
	}
	return str.str();
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
	int maxNeigh = gCfg().getInt("mutagenesis.num-seriation-neighbours");
	for(int i=-maxNeigh;i<=maxNeigh;i++){
		if(i==0)                continue;
		int nidx = idx+i;              // nidx: position of neighbour in serialization
		if(nidx<0)              continue;
		if(nidx>=(int)s.size()) continue;
		int idx = s[nidx];             // idx: original position of neighbour
		if(mTypes[idx]!="b"){
			o << ";neighChemAtom("<<mTypes[idx]<<")";
		}
	}
	bool incGraphNeigh = gCfg().getBool("mutagenesis.include-graph-neighbours");
	if(incGraphNeigh){
		AdjMat::AdjMatT& A = *mA_ptr;
		ublas::vector<double> col = ublas::column(A,idx);
		for(int i=0;i<A.size1();i++){
			if(i==idx)              continue;
			if(col(i) < 0.0001)     continue;
			o << ";graphNeigh("<<mTypes[i]<<")";
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
			mA_ptr.reset(new ublas::matrix<double>(n,n));
			cnt++;
			continue;
		}
		AdjMat::AdjMatT& A = *mA_ptr;
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
			A(i,cnt-1) = A(cnt-1,i) = (val>0 ? 1:0);
		}
		cnt++;
	}
	mNames = names;
	mTypes = types;
	mParam1 = param1;
	mParam2 = param2;

	if(n>0){
		AdjMat::AdjMatT& A = *mA_ptr;
		for(int i=0;i<n;i++)
			for(int j=i+1;j<n;j++)
				if(A(i,j)>0.001){
					A(i,j) = A(j,i) = A(i, j) - 0.5*(this->*mKernel)(i,j);
				}
	}

	ProbAdjPerm prob;
	prob.setAdjMat(mA_ptr);

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
	AdjMatGen::configure();
	string kernel = gCfg().getString("mutagenesis.kernel");
	if(0);
	else if(kernel == "kernelNull")
		mImpl->mKernel = &Impl::kernelNull;
	else if(kernel == "kernelTypeAndWeightEq")
		mImpl->mKernel = &Impl::kernelTypeAndWeightEq;
	else if(kernel == "kernelTypeSum")
		mImpl->mKernel = &Impl::kernelTypeSum;
	else if(kernel == "kernelNienhuysCheng")
		mImpl->mKernel = &Impl::kernelNienhuysCheng;
	else
		throw runtime_error("Mutagenesis: Supplied kernel unknown");
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
