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
	mSize = gCfg().getInt("fullconn_adjmat_gen.size");
	mDim  = gCfg().getInt("fullconn_adjmat_gen.dim");
}

FullConnAdjmatGen::FullConnAdjmatGen()
	:mSize(25),mDim(2),mRunningID(0),mAdj(mSize,mSize)
{
}
bool FullConnAdjmatGen::hasNext()
{
	return true;
}
FullConnAdjmatGen::~FullConnAdjmatGen()
{
}
std::string FullConnAdjmatGen::getGraphID()
{
	stringstream str;
	str<<"fullconngraph_"<<mSize<<"_"<<mRunningID;
	return str.str();
}
std::string FullConnAdjmatGen::getGraphVizNodeAttribs(int idx)
{
	return string("");
}
std::string FullConnAdjmatGen::getPlainDescription(int ser_idx, const Serialization& ser)
{
	stringstream str;
	return str.str();
}
int FullConnAdjmatGen::getClassID()
{
	return (mRunningID - (mRunningID % 30))/30;
}

ProbAdjPerm FullConnAdjmatGen::operator()()
{
	if((mRunningID%30) == 0)
	{
		mVertexPos.clear();
		for(int i=0;i<mSize;i++){
			ublas::vector<double> v(mDim);
			for(int d=0;d<mDim;d++){
				v(d) = drand48();
			}
			mVertexPos.push_back(v); 
		}
	}

	vector<ublas::vector<double> > vp = mVertexPos;
	for(int i=0;i<mSize;i++){
		ublas::vector<double>& v = vp[i];
		for(int d=0;d<mDim;d++){
			const float f = 0.1;
			v(d) += f*drand48()-0.5*f; // move vertices a bit
		}
	}
	shared_ptr<AdjMat::AdjMatT> adj(new AdjMat::AdjMatT(mSize,mSize));
	for(int i=0;i<mSize;i++)
		for(int j=0;j<mSize;j++)
		{
			if((i+j)%2 != 0)
				(*adj)(i,j) = ublas::norm_2(vp[i]-vp[j]);
			else
				(*adj)(i,j) = 0;
		}

	ProbAdjPerm pap;
	pap.setAdjMat(adj);

	//JumbledAdjMatGen jumb(pap);

	//pap = jumb();

	for(int i=0;i<adj->size1();i++)
		for(int j=0;j<adj->size2();j++)
		{
			if(fabs((*adj)(i,j)) > 10)
				cout << "arrgh!"<<endl;
			if(((*adj)(i,j)) != ((*adj)(i,j)))
				cout << "arrgh!"<<endl;
		}

	mRunningID++;
	return pap;
}

namespace{ registerInFactory<AdjMatGen, FullConnAdjmatGen> registerBase("FullConnAdjmatGen"); }
