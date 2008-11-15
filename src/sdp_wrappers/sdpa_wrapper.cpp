/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Fri Oct 24 10:00 AM 2008 CEST
 */
#include <exception>
#include <fstream>
#include <signal.h>
#include <cstdlib>
#include <sdp_prob.hpp>
#include <sdpa_wrapper.hpp>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <configuration.hpp>
#include <nana.h>

#ifndef SDPA_BINARY
#  error "You need to define SDPA_BINARY, pointing me to where sdpa resides"
#endif

using namespace std;

/*!
* helper function to print a matrix in the format required by SDPA:
*  {{1,2,3},{4,5,6},{7,8,9}}
* for
*  /1 2 3\
*  |4 5 6|
*  \7 8 9/
*/
template<class T>
void sdpaPrintMat(ofstream& o,const T& m)
{
	SDPProb::MatT::const_iterator1 row=m.begin1(),row_begin=m.begin1(),row_end=m.end1();
	SDPProb::MatT::const_iterator2 col;
	o << "{";
	for (row = m.begin1(); row != row_end; row++) {
		if (row!=row_begin) o << " ";
		o << "{";
		copy(row.begin(),row.end(),ostream_iterator<double>(o," "));
		o << "}";
	}
	o << "}"<<endl;
}

void SDPAWrapper::writeSDPAInputFile(const SDPProb& p, const char* fn)
{
	ofstream o(fn);
	o   << "\"SDPA-Wrapper Problem\""     << endl
	<< " "<<p.F.size()<<" = mDIM"         << endl     // no. of constraint matrices
	<< " 1 = nBLOCK"                      << endl     // the number of blocks in the block diagonal structure of the matrices
	<< " "<<p.C.size1()<<" = bLOCKsTRUCT" << endl     // The third line after the comments contains a vector of numbers that give the sizes of the individual blocks
	                                                  // Negative numbers may be used to indicate that a block is actually a diagonal submatrix.
	<< "{";
	if (p.b.size() > 0)                               // The fourth line after the comments contains the objective function vector c
		o << p.b(0);
	for (unsigned int i=1;i<p.b.size();i++)
		o << ", " << p.b(i);
	o   << "}"                                << endl;

	/* sparse format:
	 * The remaining lines of the file contain entries in the constraint
	 * matrices, with one entry per line.  The format for each line is 
	 *  
	 *    <matno> <blkno> <i> <j> <entry>
	 *     
	 *     Here <matno> is the number of the matrix to which this entry belongs, 
	 *     <blkno> specifies the block within this matrix, <i> and <j> specify a
	 *     location within the block, and <entry> gives the value of the entry in
	 *     the matrix.  Note that since all matrices are assumed to be symmetric, 
	 *     only entries in the upper triangle of a matrix are given.  
	 * non-sparse: 
	 * all matrices in order.
	 */
	sdpaPrintMat(o,p.C);
	// TODO: Identity-Matrix could be handled separatly (just specify diag)
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdpaPrintMat(o,m);
	}
}
bool SDPAWrapper::readSDPAOutputFile(const SDPProb& p,const char*out, AnswerT& X)
{
	ifstream is(out);
	string line;
	bool found = false;
	while( !found && getline(is,line) ){
		string start ( line.substr(0,6) );
		if( start == "xVec =" )
			found = true;
	}
	if(found == false){
		IP(found,"Could not parse SDPA output, xVec not found");
		return false;
	}
	getline(is,line); // get xVec
	line = line.substr(1,line.size()-2);


	// extract numbers as strings
	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of(",") );

	// convert strings to doubles
#if 0
	int i=0;
	boost::numeric::ublas::vector<double> y = boost::numeric::ublas::vector<double>(strvec.size());
	for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,i++){
		y(i) = boost::lexical_cast<double>(*it);
	}
#endif

	found = false;
	while( !found && getline(is,line) ){
		string start ( line.substr(0,6) );
		if( start == "yMat =" )
			found = true;
	}
	if(found == false){
		IP(found,"Could not parse SDPA output, yMat not found");
		return false;
	}
	getline(is,line); // "{"

	X = SDPWrapper::AnswerT(p.C.size1(),p.C.size2()); 
	for(unsigned int i=0;i<X.size1();i++){
		getline(is,line); 
		boost::erase_all(line,"{");
		boost::erase_all(line,"}");
		vector<string> strvec;
		boost::split( strvec, line, boost::is_any_of(",") );
		int j=0;
		for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,j++){
			boost::trim(*it);
			if(it->length() < 1) continue;
			X(i,j) = boost::lexical_cast<double>(*it);
		}
	}


	return true;
}

void SDPAWrapper::configure()
{
	mParamFile = gCfg().get<string>("sdpa-param-file");
}


void SDPAWrapper::runSDPA(const char* in, const char* out, const char* param)
{
	char cmd[255];
	sprintf(cmd,"%s -p %s -dd %s -o %s 2>&1 > /dev/null",SDPA_BINARY,param,in,out);
	L("Exec Cmd: %s...",cmd);
	int res = system(cmd);
	L("done.");
	if(res == -1)
		throw runtime_error(std::string("SDPAWrapper could not execute sdpa."));
	if(WIFSIGNALED(res) &&
			(WTERMSIG(res) == SIGINT || WTERMSIG(res) == SIGQUIT)){
		cerr << "Got interrupt, calling exit." << endl;
		exit(0);
	}
}

SDPAWrapper::AnswerT SDPAWrapper::operator()(const SDPProb&p)
{
	if(!mParamFile.length()>0)
		throw invalid_argument("SDPA-Wrapper: you need to set ParamFile!");
	const char* fn_in  = "/tmp/x.dat";
	const char* fn_out = "/tmp/x.out";
	//const char* fn_par = gCfg().getString("sdpa-param-file").c_str();
	const char* fn_par = mParamFile.c_str();
	writeSDPAInputFile(p,fn_in);
	runSDPA(fn_in,fn_out,fn_par);
	AnswerT ret;
	readSDPAOutputFile(p,fn_out,ret);
	return ret;
}
SDPAWrapper::~SDPAWrapper()
{
}

namespace
{
registerInFactory<SDPWrapper, SDPAWrapper>  registerBase("SDPAWrapper");
}
