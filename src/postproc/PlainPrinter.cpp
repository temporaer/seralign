#include <ostream>
#include <configuration.hpp>
#include <adjmat_gen.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <factory/factory.h>
#include "PlainPrinter.hpp"

using namespace std;

// 2008-11-21 Hannes Schulz <mail at hannes-schulz dot de>

PlainPrint::PlainPrint()
{
}

PlainPrint::~PlainPrint()
{
  // cleanup
}


void PlainPrint::atStart()
{
	mOS.open(gCfg().getOutputFile("output").c_str());
	if(mOS.bad())
		throw runtime_error("Cannot open output file!");
}

void PlainPrint::atEnd()
{
	mOS.close();
}

void PlainPrint::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
	Serialization::RankT ranks(ser.getRanks().size());
	for(unsigned int i=0;i<ranks.size();i++)
		ranks[i] = prob.getOriginalIndex(ser.getRanks()[i]);
	seqhead(gen);
	for(unsigned int i=0; i<ranks.size(); i++)
	{
		mOS	<< gen.getPlainDescription(i, ranks);
	}
	seqfoot();
}

void PlainPrint::seqhead(AdjMatGen& gen){
	mOS  << gen.getGraphID() <<":";
}
void PlainPrint::seqfoot(){
	mOS  <<endl; 
}

namespace{ registerInFactory<PostProc, PlainPrint> registerBase("PlainPrint"); }

