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
	string getPrologDescription(int idx,const Serialization&,const string& ref);
	string getPlainDescription(int idx,const Serialization&,const string& ref);

	double kernelNull(int i, int j, double b);
	double kernelTypeAndWeightEq(int i, int j, double b);
	double kernelTypeSum(int i, int j, double b);
	double kernelNienhuysCheng(int i, int j, double b);

	double (Mutagenesis::Impl:: *mKernel)(int, int, double);

	bool mIncludeGraphNeighbours;
	int  mNumSequenceNeighbours;
	
	struct Descriptor{
		boost::shared_ptr<AdjMat::AdjMatT> mA_ptr;
		vector<string> mNames;
		vector<string> mTypes;
		vector<string> mParam1;
		vector<string> mParam2;
		string         mName;
		int            mClassID;
	};
	map<string, Descriptor> mDescriptors;
	Descriptor              mCurDesc;
};

double Mutagenesis::Impl::kernelNienhuysCheng(int i, int j, double b){
	if (b==0) return 0;
	if(i>=mCurDesc.mTypes.size()) return 0;
	double d = 0.0;
	int    n = 0;
	if(mCurDesc.mTypes[i]  == mCurDesc.mTypes[j])  d+=1; n++;
	if(mCurDesc.mParam1[i] == mCurDesc.mParam1[j]) d+=1; n++;
	double x = pow(
			boost::lexical_cast<double>(mCurDesc.mParam2[i])-
			boost::lexical_cast<double>(mCurDesc.mParam2[j]),2.0);
	return 1.0 - 0.5*d/(2*n);
}
double Mutagenesis::Impl::kernelTypeSum(int i, int j, double b){
	if(b==0) return 0;
	if(i>=mCurDesc.mTypes.size()) return 0;
	double d = 0.0;
	map<string, double> m;
	m["n"]  = 0.3;
	m["f"]  = 0.3;
	m["cl"] = 0.3;
	m["o"]  = 0.2;
	m["h"]  = 0.1;
	m["c"]  = 0.4;

	d += fabs(m[mCurDesc.mTypes[i]] - m[mCurDesc.mTypes[j]]);

	return d;
}
double Mutagenesis::Impl::kernelTypeAndWeightEq(int i, int j, double b){
	if(b==0) return 0;
	if(i>=mCurDesc.mTypes.size()) return 0;
	double d = 0.0;
	if(mCurDesc.mTypes[i] == mCurDesc.mTypes[j]) d+=0.25;
	if(mCurDesc.mParam1[i] == mCurDesc.mParam1[j]) d+=0.25;
	return d;
}
double Mutagenesis::Impl::kernelNull(int i, int j, double b){
	return (b>0)?1.0:0.0;
}
string Mutagenesis::Impl::getPlainDescription(int ser_idx,const Serialization&s, const string&ref){
	Descriptor& d = (ref=="")?mCurDesc : mDescriptors[ref];
	stringstream str;
	AdjMat::AdjMatT& A = *d.mA_ptr;
	unsigned int idx = s.getRanks()[ser_idx];
	if(idx >= d.mTypes.size()) { str<<"Dummy"; return str.str(); }
	str << d.mTypes[idx];
	if(mIncludeGraphNeighbours){
		str << "[";
		string s = "";
		ublas::vector<double> col = ublas::column(A,idx);
		for(unsigned int i=0;i<A.size1();i++){
			if(i==idx)              continue;
			if(col(i) < 0.0001)     continue;
				s += d.mTypes[i];
		}
		sort(s.begin(),s.end());
		str << s;
		str<< "]";
	}
	return str.str();
}
string Mutagenesis::Impl::getPrologDescription(int ser_idx,const Serialization&s, const string& ref){
	Descriptor& d = (ref=="")?mCurDesc : mDescriptors[ref];
	unsigned int idx = s.getRanks()[ser_idx];
	if(idx >= d.mTypes.size()) { return string("dummy"); }
	ostringstream o;
	if(d.mTypes[idx] == "b"){
		// bond
		o << "chemBond("<<d.mParam1[idx] << ")";
	}else{
		// atom
		//o << "chemAtom(" << mName <<","<< mNames[idx] << "," << mTypes[idx] << ")";
		o << "chemAtom(" << d.mTypes[idx] << ")";
		o << ";weight(" << d.mParam1[idx] << ")";
	}
	for(int i=-mNumSequenceNeighbours;i<=mNumSequenceNeighbours;i++){
		if(i==0)                continue;
		int nidx = idx+i;              // nidx: position of neighbour in serialization
		if(nidx<0)              continue;
		if(nidx>=(int)s.getRanks().size()) continue;
		int idx = s.getRanks()[nidx];             // idx: original position of neighbour
		if(d.mTypes[idx]!="b"){
			o << ";neighChemAtom("<<d.mTypes[idx]<<")";
		}
	}
	if(mIncludeGraphNeighbours){
		AdjMat::AdjMatT& A = *d.mA_ptr;
		ublas::vector<double> col = ublas::column(A,idx);
		for(unsigned int i=0;i<A.size1();i++){
			if(i==idx)              continue;
			if(col(i) < 0.0001)     continue;
			o << ";graphNeigh("<<d.mTypes[i]<<")";
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
	int n = -1, actualn = -1;

	getline(mInputStream, line);
	mCurDesc.mName = line;
	getline(mInputStream, line);
	mCurDesc.mClassID = (line == "pos") ? 1 : 0;

	while(line.length()>1 && getline(mInputStream, line)){
		list<string> strvec;
		boost::split( strvec, line, boost::is_any_of(",") );
		int addParams=4;
		if(cnt==0){
			for(int i=0;i<addParams;i++) strvec.pop_front(); // ignore 1st 4 "None"s
			list<string>::iterator it = strvec.begin();
			copy(it,strvec.end(),back_inserter(names));
			actualn = names.size();
			n = 40;
			mCurDesc.mA_ptr.reset(new ublas::matrix<double>(n,n));
			*mCurDesc.mA_ptr = ublas::scalar_matrix<double>(n,n,0.0);
			cnt++;
			continue;
		}
		AdjMat::AdjMatT& A = *mCurDesc.mA_ptr;
		if(actualn != strvec.size()-addParams){
			if(cnt-1 == actualn)
				break;
			L("Breaking at cnt=%d n=%d A.size1()=%d, strvecsize=%d\n",cnt,actualn,A.size1(),strvec.size());
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
			//A(i,cnt-1) = A(cnt-1,i) = (val>0 ? 1:0);
			A(i,cnt-1) = A(cnt-1,i) = val;
		}
		cnt++;
	}
	mCurDesc.mNames = names;
	mCurDesc.mTypes = types;
	mCurDesc.mParam1 = param1;
	mCurDesc.mParam2 = param2;

	if(n>0){
		AdjMat::AdjMatT& A = *mCurDesc.mA_ptr;
		for(int i=0;i<n;i++)
			for(int j=i+1;j<n;j++)
				A(i,j) = A(j,i) = (this->*mKernel)(i,j, A(i,j));
	}

	ProbAdjPerm prob;
	prob.setAdjMat(mCurDesc.mA_ptr);
	prob.setId(mCurDesc.mName);
	mDescriptors[mCurDesc.mName] = mCurDesc;

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
	mImpl->mIncludeGraphNeighbours = gCfg().getBool("mutagenesis.include-graph-neighbours");
	mImpl->mNumSequenceNeighbours  = gCfg().getInt ("mutagenesis.num-seriation-neighbours");
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
string Mutagenesis::getPlainDescription(int ser_idx, const Serialization& s, const string& ref)
{
	return mImpl->getPlainDescription(ser_idx,s,ref);
}
string Mutagenesis::getPrologDescription(int ser_idx,const Serialization& s, const string& ref)
{
	return mImpl->getPrologDescription(ser_idx,s,ref);
}
boost::shared_ptr<AdjMat::AdjMatT>
Mutagenesis::getAdjMat(const string&ref)
{
	if(ref!="")
		return mImpl->mDescriptors[ref].mA_ptr;
	return mImpl->mCurDesc.mA_ptr;
}
std::string Mutagenesis::getGraphID(const string&ref)
{
	if(ref!="")
		return mImpl->mDescriptors[ref].mName;
	return mImpl->mCurDesc.mName;
}
int Mutagenesis::getClassID(const string&ref)
{
	if(ref!="")
		return mImpl->mDescriptors[ref].mClassID;
	return mImpl->mCurDesc.mClassID;
}


namespace{ registerInFactory<AdjMatGen, Mutagenesis> registerBase("Mutagenesis"); }
