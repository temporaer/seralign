#include <fstream>
#include <signal.h>
#include <cstdlib>
#include <sdp_prob.hpp>
#include <sdplr_wrapper.hpp>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <configuration.hpp>
#include <nana.h>

#ifndef SDPLR_BINARY
#  error "You need to define SDPLR_BINARY, pointing me to where sdplr resides"
#endif

using namespace std;

/*!
* helper function to print a matrix in the format required by SDPLR
*/
template<class T>
void sdplrPrintMat(ofstream& o,int n, const T& m)
{
	// 1st line:
	// - number of matrix (0 for C)
	// - number of block (>=1)
	// - "s" for sparse
	// - number of entries
	o << n << " 1 s " << (m.size1()*(m.size1()+1)/2) << endl;
	for(unsigned int i=0;i<m.size1();i++)
		for(unsigned int j=i;j<m.size2();j++)
			o << i+1 << " " << j+1 << " " << m(i,j)<< endl;
}
/*!
* helper function to print a matrix in SDPA Sparse format
*/
template<class T>
void sdpaSparsePrintMat(ofstream& o,int n, const T& m)
{
	for(unsigned int i=0;i<m.size1();i++)
		for(unsigned int j=i;j<m.size2();j++)
			o << n <<" 1 " << i+1 << " " << j+1 << " " << m(i,j)<< endl;
}

void writeSDPASparseInputFile(const SDPProb& p, const char* fn)
{
	ofstream o(fn);
	o << "\"SDPASparse Sample Problem\"" << endl
	<< p.F.size()                        << endl  // num of constraint matr
	<< "1"                               << endl  // number of blocks in SDP
	<< p.C.size1()                       << endl; // size of the blocks from prev line (one line per block)
	if (p.b.size() > 0)                            // b, the righthandside vec, all on one line
		o << p.b(0);
	for (unsigned int i=1;i<p.b.size();i++)
		o << " " << p.b(i);
	o << endl;

	sdpaSparsePrintMat(o,0, p.C);
	// TODO: Identity-Matrix could be handled separatly (just specify diag)
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdpaSparsePrintMat(o,f+1,m);
	}
}
void SDPLRWrapper::writeSDPLRInputFile(const SDPProb& p, const char* fn)
{
	ofstream o(fn);
	o << p.F.size()                     << endl  // num of constraint matr
	<< "1"                               << endl  // number of blocks in SDP
	<< p.C.size1()                       << endl; // size of the blocks from prev line (one line per block)
	if (p.b.size() > 0)                            // b, the righthandside vec, all on one line
		o << p.b(0);
	for (unsigned int i=1;i<p.b.size();i++)
		o << " " << p.b(i);
	o << endl;
	o << "1" <<endl;   // ignored

	sdplrPrintMat(o,0, p.C);
	// TODO: Identity-Matrix could be handled separatly (just specify diag)
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdplrPrintMat(o,f+1,m);
	}
}
bool SDPLRWrapper::readSDPLROutputFile(const char*out, AnswerT& ret)
{
	ifstream is(out);

	vector<string> strvec;
	string s;
	int n;
	is >> s >> s;
	is >> n;
	for(int i=0;i<n;i++){
		is>>s;
		strvec.push_back(s);
	}

	// convert strings to doubles
	int i=0;
	boost::numeric::ublas::vector<double> y = boost::numeric::ublas::vector<double>(strvec.size());
	for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,i++){
		y(i) = boost::lexical_cast<double>(*it);
	}

	return true;
}

void SDPLRWrapper::runSDPLR(const char* in, const char* out, const char* param)
{
	char cmd[255];
	// use dummy for soln_in as described in manual
	sprintf(cmd,"%s %s %s soln_in %s 2>&1 > /dev/null",SDPLR_BINARY,in,param,out);
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

SDPLRWrapper::AnswerT SDPLRWrapper::operator()(const SDPProb&p)
{
	L("SDPLRWrapper::operator()\n");
	const char* fn_in  = "/tmp/x.dat";
	const char* fn_out = "/tmp/x.out";
	//const char* fn_par = gCfg().getString("sdplr-param-file").c_str();
	const char* fn_par = mParamFile.c_str();
	writeSDPASparseInputFile(p,fn_in);
	runSDPLR(fn_in,fn_out,fn_par);
	AnswerT ret;
	readSDPLROutputFile(fn_out,ret);
	return ret;
}
SDPLRWrapper::~SDPLRWrapper()
{
	L("Destroying SDPLRWrapper\n");
}

namespace
{
registerInFactory<SDPWrapper, SDPLRWrapper>  registerBase("SDPLRWrapper");
}
