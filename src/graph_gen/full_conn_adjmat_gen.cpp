#include <sstream>
#include <configuration.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <factory/factory.h>
#include <matlab_io.hpp>
#include "full_conn_adjmat_gen.hpp"
#include "jumbled_adjmat_gen.hpp"

using namespace std;
using namespace boost;

void FullConnAdjmatGen::configure()
{
	AdjMatGen::configure(); 
}

FullConnAdjmatGen::FullConnAdjmatGen()
	:mSize(25),mRunningID(0)
{
}
bool FullConnAdjmatGen::hasNext()
{
	return mRunningID<10;
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
	switch(ser[ser_idx]){
		case 3: str << ser[ser_idx]; break;
		case 5: str << ser[ser_idx]; break;
		default:str << "_";
	}
	return str.str();
}
ProbAdjPerm FullConnAdjmatGen::operator()()
{
	ProbAdjPerm pap;
	shared_ptr<AdjMat::AdjMatT> adj(new AdjMat::AdjMatT(mSize, mSize));
	for(int i=0;i<mSize;i++)
		for(int j=i;j<mSize;j++)
		{
			if(j==i){
				(*adj)(j,i) = (*adj)(i,j) = 0;
				continue;
			}
			if(0);
			else if(i==3 || j==3)
				(*adj)(j,i) = (*adj)(i,j) = drand48()*0.3; 
			else if(i==5 || j==5)
				(*adj)(j,i) = (*adj)(i,j) = drand48()*0.1; 
			else
				(*adj)(j,i) = (*adj)(i,j) = drand48()*1;
		}


	pap.setAdjMat(adj);

	JumbledAdjMatGen jumb(pap);

	pap = jumb();

	return pap;
}

namespace{ registerInFactory<AdjMatGen, FullConnAdjmatGen> registerBase("FullConnAdjmatGen"); }
