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
	Serialization cpy = ser;
	for(unsigned int i=0;i<ser.size();i++)
		cpy[i] = prob.getOriginalIndex(ser[i]);

	fs::path dir(gCfg().getString("output-dir").c_str());
	fs::path base_fn(gen.getGraphID() + ".dot");
	mOS.open((dir/base_fn).string().c_str());
	unsigned int n=prob.getAdjMat()->size1();
	seqhead(gen, n);
	// print node labels
	for(unsigned int i=0; i<n; i++)
	{
		mOS	<< "  n" << ser[i] << "[" // use ID in adj-matrix
			<<"label=\"("<<i<<") "<<gen.getPlainDescription(i, cpy)<<"\""
			<<gen.getGraphVizNodeAttribs(cpy[i])<<"];"<<endl;
	}
	AdjMat::AdjMatT& A = *prob.getAdjMat();
	for(unsigned int i=0;i<n;i++)
		for(unsigned int j=i;j<n;j++)
		{
			if(A(i,j)>0.0001){
				mOS << "  n"<<i<<" -- n"<<j<<" [weight="<<A(i,j)<<"];"<<endl;
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

