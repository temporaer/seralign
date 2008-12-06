#include <memory>
#include <numeric>
#include "repetative_seriation_gen.hpp"
#include <configuration.hpp>
#include <jumbled_adjmat_gen.hpp>
#include <DegreeSort.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <factory/factory.h>
#include <matlab_io.hpp>
#include <nana.h>

using namespace std;
namespace ublas = boost::numeric::ublas;

// 2008-12-02 Hannes Schulz <mail at hannes-schulz dot de>

RepetativeSerGen::RepetativeSerGen()
{
}

void RepetativeSerGen::configure()
{
	SerGenAdj::configure();
	this->setSerializer(gCfg().getString("repetative-ser-gen.serializer"));
	this->setRepetitions(gCfg().getInt  ("repetative-ser-gen.repetition_num"));
}

struct VectorSorter{
	const ublas::vector<double>& v;
	VectorSorter(const ublas::vector<double>& x)
		:v(x){}
	bool operator()(int i, int j){
		I(i>=0 && j>=0);
		I(i<(int)v.size() && j<(int)v.size());
		return  v(i) < v(j) ;
	}
};

Serialization RepetativeSerGen::operator()(const ProbAdjPerm& pap)
{
	ProbAdjPerm mypap = pap;
	JumbledAdjMatGen jam(mypap);
	jam.configure();
	
	auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(mSerializer);
	if(!seriation_gen.get())
		throw logic_error(string("Supplied SerGenAdj `") + mSerializer + "' does not exist");
	seriation_gen->configure();

	int n = mypap.getAdjMat()->size1();
	
	ublas::matrix<double> positions(n, mRepetitions);
	
	Serialization ser(n);
	for(int rep=0;rep<mRepetitions;rep++){
		LG(isVerbose(),"RepetativeSerGen: Serializing %i of %i\n", rep+1, mRepetitions);
		mypap = jam();
		DegreeSort ds;
		ds.sort(mypap);
		ser = (*seriation_gen)(mypap);
		for(int i=0;i<n;i++)
			positions(mypap.getOriginalIndex(ser(i)), rep) = i;
	}
	ublas::vector<double> avgpos(n);
	for(int i=0;i<n;i++){
		ublas::matrix_row<ublas::matrix<double> > row(ublas::row(positions,i));
		avgpos(i) = accumulate(row.begin(),row.end(), 0.0);
		if(isVerbose())
			cout << "avgpos("<<i<<") = "<<( avgpos(i)/mRepetitions )<<endl;
	}
	//matlab_matrix_out(cout, "positions", positions);

	Serialization finalser(n);
	for(int i=0;i<n;i++)
		finalser(i)=i;

	sort(finalser.begin(), finalser.end(), VectorSorter(avgpos));

	return finalser;
}


RepetativeSerGen::~RepetativeSerGen()
{
  // cleanup
}

namespace{ registerInFactory<SerGenAdj, RepetativeSerGen> registerBase("RepetativeSerGen"); }
