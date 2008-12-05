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
	Serialization cpy = ser;
	for(unsigned int i=0;i<ser.size();i++)
		cpy[i] = prob.getOriginalIndex(ser[i]);
	seqhead(gen);
	for(unsigned int i=0; i<cpy.size(); i++)
	{
		mOS	<<"    <atom>" <<endl
			<<"      <symbol>" << gen.getPrologDescription(i,cpy) << "</symbol>" << endl
			<<"      <label>none</label>"   <<endl
			<<"    </atom>"<<endl;
	}
	seqfoot();
}

void RealXMLPrint::seqhead(AdjMatGen& gen){
	mOS  << "  <sequence id=\"" << gen.getGraphID() <<"\">"<<endl; 
}
void RealXMLPrint::seqfoot(){
	mOS  << "  </sequence>"<<endl; 
}

namespace{ registerInFactory<PostProc, RealXMLPrint> registerBase("RealXMLPrint"); }
