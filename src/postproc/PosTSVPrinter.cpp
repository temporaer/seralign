#include <ostream>
#include <configuration.hpp>
#include <adjmat_gen.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <factory/factory.h>
#include "PosTSVPrinter.hpp"

using namespace std;

// 2008-11-21 Hannes Schulz <mail at hannes-schulz dot de>

PosTSVPrint::PosTSVPrint()
{
}

PosTSVPrint::~PosTSVPrint()
{
  // cleanup
}


void PosTSVPrint::atStart()
{
	mOS.open(gCfg().getOutputFile("output").c_str());
	if(mOS.bad())
		throw runtime_error("Cannot open output file!");
	for(int i=0;i<gCfg().getInt("fullconn_adjmat_gen.size")-1;i++)
		mOS << "s"<<i<< ",";
	mOS << "class";
	mOS<<endl;
}

void PosTSVPrint::atEnd()
{
	mOS.close();
}

void PosTSVPrint::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
	//Serialization::RankT ranks(ser.getRanks().size());
	//for(unsigned int i=0;i<ranks.size();i++)
		//ranks[i] = prob.getOriginalIndex(ser.getRanks()[i]);
	seqhead(gen);
	//for(unsigned int i=0; i<ranks.size(); i++)
		//mOS	<< gen.getPlainDescription(i, ranks);
	for(unsigned int i=0; i<ser.getPositions().size()-1; i++)
		mOS	<< ser.getPositions()(i+1) - ser.getPositions()(i) << ",";
	boost::any ref = prob.getBackground();
	mOS  << ((gen.getClassID(ref)==0)?"A":"B") <<endl;
	seqfoot();
}

void PosTSVPrint::seqhead(AdjMatGen& gen){
}
void PosTSVPrint::seqfoot(){
}

namespace{ registerInFactory<PostProc, PosTSVPrint> registerBase("PosTSVPrint"); }

