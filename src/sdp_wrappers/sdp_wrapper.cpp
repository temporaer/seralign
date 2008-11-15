#include <exception>
#include <sdp_wrapper.hpp>
#include <factory/factory.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <nana.h>

using namespace std;

SDPWrapper::AnswerT SDPWrapper::operator()(const SDPProb&)
{
	throw logic_error("Called SDPWrapper() w/o subclassing");
}
SDPWrapper::~SDPWrapper(){
}
void SDPWrapper::configure()
{
}


namespace{
	registerInFactory<SDPWrapper, SDPWrapper> registerBase("SDPWrapper");
}
