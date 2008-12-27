#include <boost/numeric/ublas/matrix.hpp>
#include <seqan/align.h>
#include <Serialization.hpp>
#include <factory/factory.h>
#include <adjmat_gen.hpp>
#include <ProbData.hpp>
#include <configuration.hpp>
#include <matlab_io.hpp>
#include <progressbar.hpp>
#include <stats.hpp>
#include "localign_distmat.hpp"
LocAlignDistmat::Node::Node(const std::string& t)
	:s(t)
{
}
LocAlignDistmat::Node::Node()
{
}

LocAlignDistmat::LocAlignDistmat()
{
}
void LocAlignDistmat::atSeriation(AdjMatGen& gen, Serialization& ser, ProbAdjPerm& prob)
{
	Serialization::RankT ranks(ser.getRanks().size());
	for(unsigned int i=0;i<ranks.size();i++)
		ranks[i] = prob.getOriginalIndex(ser.getRanks()[i]);
	seqan::String<Node> str;
	for(unsigned int i=0; i<ranks.size(); i++)
	{
		seqan::appendValue(str, Node(gen.getPlainDescription(i, ranks)));
	}
	mGraphIDs.push_back(gen.getGraphID());
	mStrings.push_back(str);
}
void LocAlignDistmat::atStart(){
}
void LocAlignDistmat::atEnd()
{
	unsigned int n = mStrings.size(); 
	DistMatT distmat(n,n);
	ExactDescriptiveStatistics stats("score");

	ProgressBar pb(n, "aligning");
	for(unsigned int i=0;i<mStrings.size();i++){
		pb.inc();
		for(unsigned int j=i;j<mStrings.size();j++)
		{
			seqan::Align<seqan::String<Node> > align;
			if(seqan::length(mStrings[i]) < seqan::length(mStrings[j])){
				seqan::appendValue(seqan::rows(align), mStrings[j]);
				seqan::appendValue(seqan::rows(align), mStrings[i]);
			}else{
				seqan::appendValue(seqan::rows(align), mStrings[i]);
				seqan::appendValue(seqan::rows(align), mStrings[j]);
			}
			int score = seqan::localAlignment(align, seqan::Score<int>(7,-2,-2), seqan::SmithWaterman());
			stats += score;
			distmat(j,i) = distmat(i,j) = score;
		}
	}
	pb.finish();
	for(unsigned int i=0;i<n;i++)
		for(unsigned int j=i;j<n;j++)
			distmat(j,i) = distmat(i,j) = exp(-(distmat(i,j)-stats.getMin())/stats.getRange());
	std::string fn = gCfg().getOutputFile("output");

	// output class header
	std::ofstream os(fn.c_str());
	os << " ,"<<mGraphIDs[0] <<"";
	for(unsigned int i=1;i<n;i++)
		os << ","<<mGraphIDs[i]<<"";
	os << std::endl;
	csv_matrix_out(os,distmat);
	os.close();
}
LocAlignDistmat::~LocAlignDistmat()
{
}


namespace{ registerInFactory<PostProc, LocAlignDistmat> registerBase("LocAlignDistmat"); }
