#include <random_adjmat_gen.hpp>
#include <iostream>
#include <fstream>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <configuration.hpp>
#include <nana.h>

using namespace boost;
using namespace std;
void RandomAdjMatGen::configure()
{
	this->setMatrixSize(gCfg().getInt("rand_adj_mat_gen.size"));
	this->setConnectionProb(
			1.0/(mMatrixSize-1)
			*gCfg().getFloat("rand_adj_mat_gen.out-degree")
			);
	this->setWeightedEdges(gCfg().getBool("rand_adj_mat_gen.weighted"));
	this->setSeed(gCfg().getFloat("rand_adj_mat_gen.seed"));
}

shared_ptr<RandomAdjMatGen::AdjMatT> RandomAdjMatGen::operator()()
{
	L("Generating Random AdjMat with n=%d p=%2.2f w=%d...",mMatrixSize,mConnectionProb,mWeightedEdges);
	shared_ptr<AdjMatT> adjmat_ptr(new AdjMatT(mMatrixSize,mMatrixSize));
	
	AdjMatT& adjmat = *adjmat_ptr;
	int n = adjmat.size1();
	srand48(mSeed);
	for(int i=0;i<n;i++)
        for(int j=i;j<n;j++){
			if(mWeightedEdges){
				if(drand48()>mConnectionProb)
					adjmat(i,j) = 0;
				else
					adjmat(i,j) = drand48();
			}else{
				adjmat(i,j) = drand48()<mConnectionProb;
			}
            adjmat(j,i) = adjmat(i,j);
		}
	L("done.\n");
	
	ofstream o("adjmat.dat");
	o<<"A = [ ";
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			o << adjmat(i,j) <<" ";
		}
		o << "; ";
	}
	o << " ]; \n";
	return adjmat_ptr;
}
RandomAdjMatGen::~RandomAdjMatGen()
{
}

namespace{
	registerInFactory<AdjMatGen, RandomAdjMatGen> registerBase("RandomAdjMatGen");
}
