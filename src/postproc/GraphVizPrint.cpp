#include <ostream>
#include <configuration.hpp>
#include <adjmat_gen.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/filesystem/path.hpp>
#include <factory/factory.h>
#include "GraphVizPrint.hpp"

using namespace std;
namespace fs = boost::filesystem;

// 2008-11-21 Hannes Schulz <mail at hannes-schulz dot de>

GraphVizPrint::GraphVizPrint()
{
}

GraphVizPrint::~GraphVizPrint()
{
  // cleanup
}


void GraphVizPrint::atStart()
{
}

void GraphVizPrint::atEnd()
{
}

void GraphVizPrint::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
	fs::path dir(gCfg().getString("output-dir").c_str());
	fs::path base_fn(gen.getGraphID() + ".dot");
	mOS.open((dir/base_fn).string().c_str());
	unsigned int n=prob.getAdjMat()->size1();
	seqhead(gen, n);
	// print node labels
	for(unsigned int i=0; i<n; i++)
	{
		int idx = prob.getOriginalIndex(i);
		int sernr = -1;
		for(unsigned int s=0;s<ser.size();s++)
			if(ser[s]==i)
				sernr=s;

		mOS	<< "  n" << i << "["
			<<"label=\"("<<sernr<<") "<<gen.getPlainDescription(idx)<<"\""
			<<gen.getGraphVizNodeAttribs(idx)<<"];"<<endl;
	}
	AdjMat::AdjMatT& A = *prob.getAdjMat();
	for(unsigned int i=0;i<n;i++)
		for(unsigned int j=i;j<n;j++)
		{
			if(A(i,j)>0.0001){
				mOS << "  n"<<i<<" -- n"<<j<<" ;"<<endl;
			}
		}
	seqfoot();
	mOS.close();
}

void GraphVizPrint::seqhead(AdjMatGen& gen, int size){
	mOS << "graph "<< gen.getGraphID()<<" {"<<endl
		<< "  size = \""<<size<<","<<size<<"\";"<<endl
		<< "  label = \""<<gen.getGraphID()<<"\\n\\n\";"<<endl
		;
}
void GraphVizPrint::seqfoot(){
	mOS  <<"}"<<endl; 
}

namespace{ registerInFactory<PostProc, GraphVizPrint> registerBase("GraphVizPrint"); }

