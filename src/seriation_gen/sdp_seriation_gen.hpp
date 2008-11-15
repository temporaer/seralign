#ifndef __SDP_SERIATION_GEN_HPP__
#define __SDP_SERIATION_GEN_HPP__
#include <memory>
#include <SerGenAdj.hpp>

class SDPWrapper;

class SDPSeriationGen : public SerGenAdj
{
	private:
		struct Impl;
		boost::shared_ptr<Impl> mImpl;

	public:
        SDPSeriationGen();
		void setSDPWrapper(std::auto_ptr<SDPWrapper> w);
		virtual Serialization operator()(const ProbAdjPerm& pap);
		virtual void configure();
		virtual ~SDPSeriationGen();
};
#endif /* __SDP_SERIATION_GEN_HPP__ */

