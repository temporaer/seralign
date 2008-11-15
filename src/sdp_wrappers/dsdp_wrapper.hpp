#ifndef __DSDP_WRAPPER_HPP__
#define __DSDP_WRAPPER_HPP__
#include <boost/shared_ptr.hpp>
#include <sdp_wrapper.hpp>

class DSDPWrapper: public SDPWrapper
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;

	public:
		DSDPWrapper();
		virtual ~DSDPWrapper();
		virtual AnswerT operator()(const SDPProb&);
};
#endif

