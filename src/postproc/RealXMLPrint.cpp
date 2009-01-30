#include <ostream>
#include <configuration.hpp>
#include <adjmat_gen.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <factory/factory.h>
#include "RealXMLPrint.hpp"

using namespace std;

// 2008-11-21 Hannes Schulz <mail at hannes-schulz dot de>

RealXMLPrint::RealXMLPrint()
{
}

RealXMLPrint::~RealXMLPrint()
{
  // cleanup
}


void RealXMLPrint::atStart()
{
	mOS.open(gCfg().getOutputFile("output").c_str());
	cout << "Writing to Output: "<< gCfg().getOutputFile("output").c_str()<<endl;
	if(mOS.bad())
		throw runtime_error("Could not open output file!");
	mOS<< "<sequences>"<<endl; 
}

void RealXMLPrint::atEnd()
{
	mOS<< "</sequences>"<<endl; 
	mOS.close();
}

void RealXMLPrint::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
	Serialization::RankT ranks = ser.getRanks();
	for(unsigned int i=0;i<ranks.size();i++)
		ranks[i] = prob.getOriginalIndex(ser.getRanks()[i]);
	mOS  << "  <sequence id=\"" << gen.getGraphID(prob.getBackground()) <<"\">"<<endl; 
	for(unsigned int i=0; i<ranks.size(); i++)
	{
		mOS	<<"    <atom>" <<endl
			<<"      <symbol>" << gen.getPrologDescription(i,ranks, prob.getBackground()) << "</symbol>" << endl
			<<"      <label>none</label>"   <<endl
			<<"    </atom>"<<endl;
	}
	seqfoot();
}

void RealXMLPrint::seqhead(AdjMatGen& gen){
}
void RealXMLPrint::seqfoot(){
	mOS  << "  </sequence>"<<endl; 
}

namespace{ registerInFactory<PostProc, RealXMLPrint> registerBase("RealXMLPrint"); }
