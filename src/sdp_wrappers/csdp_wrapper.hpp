#ifndef __CSDP_WRAPPER_HPP__
#define __CSDP_WRAPPER_HPP__
#include <sdp_wrapper.hpp>
class CSDPWrapper: public SDPWrapper
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;
	public:
		CSDPWrapper();
		virtual AnswerT operator()(const SDPProb&);
		virtual ~CSDPWrapper();
};
#endif

