#include <sstream>
#include <configuration.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <factory/factory.h>
#include <matlab_io.hpp>
#include "full_conn_adjmat_gen.hpp"
#include "jumbled_adjmat_gen.hpp"
#include <nana.h>

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;

void FullConnAdjmatGen::configure()
{
	AdjMatGen::configure(); 
	mSize    = gCfg().getInt("fullconn_adjmat_gen.size");
	mPatSize = gCfg().getInt("fullconn_adjmat_gen.patsize");
	mJumble  = gCfg().getBool("fullconn_adjmat_gen.jumble");

	mPatClass0 = ublas::scalar_matrix<int>(mSize, mSize,0);
	mPatClass1 = ublas::scalar_matrix<int>(mSize, mSize,0);

	int j;
	//for(int i=0;i<mPatSize; i++){
		//do { j = drand48() * mPatSize; } while(j==i);
		//mPatClass0(j,i) = mPatClass0(i,j) = 1;

		//do { j = drand48() * mPatSize; } while(j==i);
		//mPatClass1(j,i) = mPatClass1(i,j) = 1;
	//}

	for(int i=0;i<mPatSize; i++){
		do { j = drand48() * mPatSize; } while(j==i);
		mPatClass0(j,i) = mPatClass0(i,j) = 2;
		mPatClass1(j,i) = mPatClass1(i,j) = 1;
	}
}

FullConnAdjmatGen::FullConnAdjmatGen()
	:mSize(25),mRunningID(0),mPatSize(10),mAdj(mSize,mSize)
{
	mPatClass0 = AdjMat::AdjMatT(mSize, mSize);
	mPatClass1 = AdjMat::AdjMatT(mSize, mSize);
}
bool FullConnAdjmatGen::hasNext()
{
	return true;
}
FullConnAdjmatGen::~FullConnAdjmatGen()
{
}
std::string FullConnAdjmatGen::getGraphID(const string& ref)
{
	I(ref == "");
	stringstream str;
	str<<"fullconngraph_"<<mSize<<"_"<<mRunningID;
	return str.str();
}
std::string FullConnAdjmatGen::getGraphVizNodeAttribs(int idx, const string& ref)
{
	I(ref == "");
	return string("");
}
std::string FullConnAdjmatGen::getPlainDescription(int ser_idx, const Serialization& ser, const string& ref)
{
	I(ref == "");
	stringstream str;
	return str.str();
}
int FullConnAdjmatGen::getClassID(const string& ref)
{
	if(ref=="")
		return mRunningID % 2;
	return mDescriptors[ref].mClassID;
}

ProbAdjPerm FullConnAdjmatGen::operator()()
{
	mRunningID++;
	shared_ptr<AdjMat::AdjMatT> adj(new AdjMat::AdjMatT(mSize,mSize));

	Descriptor desc;
	desc.mA_ptr = adj;
	
	// determine what to copy
	if(mRunningID%2){
		*adj = mPatClass1;
		desc.mClassID = 1;
	}
	else{
		*adj = mPatClass0;
		desc.mClassID = 0;
	}

	// save info
	mDescriptors[getGraphID("")] = desc;

	// modify
	int j;
	for(int i=mPatSize; i<mSize; i++){
		do { j = drand48() * (mSize-mPatSize); } while(j==i);
		I(mPatSize+j >=mPatSize);
		I(mPatSize+j < mSize);
		int val = drand48() > 0.5 ? 2 : 1;
		(*adj)(i,mPatSize+j) = val;
		(*adj)(mPatSize+j,i) = val;
	}


	ProbAdjPerm pap;
	pap.setAdjMat(adj);
	pap.setId(getGraphID(""));

	if(mJumble){
		JumbledAdjMatGen jumb(pap);
		pap = jumb();
	}

	return pap;
}

namespace{ registerInFactory<AdjMatGen, FullConnAdjmatGen> registerBase("FullConnAdjmatGen"); }
