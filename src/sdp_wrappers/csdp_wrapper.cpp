/*       Created   :  10/06/2008 12:52:01 AM
 *       Last Change: Tue Oct 28 01:00 AM 2008 CET
 */

#include <fstream>
#include <signal.h>
#include <csdp_wrapper.hpp>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <nana.h>

#ifndef CSDP_BINARY
#  error "You need to define CSDP_BINARY, pointing me to where csdp resides"
#endif

using namespace std;

struct CSDPWrapper::Impl{
	typedef CSDPWrapper::AnswerT AnswerT;
	AnswerT operator()(const SDPProb&);
	void writeSDPASparseInputFile(const SDPProb& p, const char* fn);
};

/*!
* helper function to print a matrix in SDPA Sparse format
*/
template<class T>
void sdpaSparsePrintMat(ofstream& o,int n, const T& m)
{
	for(unsigned int i=0;i<m.size1();i++)
		for(unsigned int j=i;j<m.size2();j++){
			I(m(i,j) == m(i,j)); // check for NaN
			o << n << " 1 " << (i+1) << " " << (j+1) << " " << m(i,j)<< endl;
		}
}

void CSDPWrapper::Impl::writeSDPASparseInputFile(const SDPProb& p, const char* fn)
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
	//       for example, specialize sdpaSparsePrintMat(...) for identity-matrix
	for (unsigned int f=0;f<p.F.size();f++) {
	    const SDPProb::MatT& m = p.F[f];
	    sdpaSparsePrintMat(o,f+1,m);
	}
}

bool readCSDPOutputFile(const SDPProb& p, const char* out, CSDPWrapper::AnswerT& X){
	ifstream is(out);
	string line;
	getline(is,line);
	boost::trim(line);

	vector<string> strvec;
	boost::split( strvec, line, boost::is_any_of(" ") );
	
	// convert strings to doubles
	int i=0;
#if 0
	ublas::vector<double> y(strvec.size());
	for(vector<string>::iterator it = strvec.begin();it!=strvec.end();it++,i++){
		y(i) = boost::lexical_cast<double>(*it);
	}
#endif

	X = SDPWrapper::AnswerT(p.C.size1(),p.C.size2()); 
	boost::numeric::ublas::symmetric_adaptor<SDPWrapper::AnswerT,
	boost::numeric::ublas::upper > X_(X);
	int n,b,j;
	string num;
	while(!is.eof()){
		is >> n;
		
		if(n==1) {getline(is,line); continue;}
		IP(n==2, "Syntax Error in DSDP Output File!");
		is >> b; // block
		is >> i; // idx
		is >> j; // idx
		is >> num;
		X(i-1,j-1) = boost::lexical_cast<double>(num);
		X(j-1,i-1) = X(i-1,j-1);
	}
	return true;
}

void runCSDP(const char* in, const char* out)
{
	char cmd[255];
	sprintf(cmd,"%s %s %s 2>&1 > /dev/null",CSDP_BINARY,in,out);
	L("Exec Cmd: %s...",cmd);
	int res = system(cmd);
	L("done.\n");
	if(res == -1)
		throw runtime_error(std::string("CSDPWrapper could not execute csdp."));
	if(WIFSIGNALED(res) &&
			(WTERMSIG(res) == SIGINT || WTERMSIG(res) == SIGQUIT)){
		cerr << "Got interrupt, calling exit." << endl;
		exit(0);
	}
}


CSDPWrapper::AnswerT CSDPWrapper::Impl::operator()(const SDPProb& p){
	const char* fn_in  = "/tmp/x.dat";
	const char* fn_out = "/tmp/x.out";
	writeSDPASparseInputFile(p,fn_in);
	runCSDP(fn_in,fn_out);
	CSDPWrapper::AnswerT ret;
	readCSDPOutputFile(p, fn_out,ret);
	return ret;
}
CSDPWrapper::CSDPWrapper()
	:mImpl(new Impl)
{
}


CSDPWrapper::AnswerT CSDPWrapper::operator()(const SDPProb& prob)
{
	return (*mImpl)(prob);
}
CSDPWrapper::~CSDPWrapper(){
}

namespace{
	registerInFactory<SDPWrapper, CSDPWrapper>  registerBase("CSDPWrapper");
}
