#ifndef _TIMER_H_
#define _TIMER_H_

#include <boost/date_time/posix_time/posix_time.hpp>

namespace AprilTags
{

typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;

	class Timer
	{
	public:

		Timer()
		{}
		
		void Start()
		{
			lastTime = boost::posix_time::microsec_clock::universal_time();
		}
		
		void Lap( const std::string& pt )
		{
			Time now = boost::posix_time::microsec_clock::universal_time();
			TimeDuration diff = now - lastTime;
			std::cout << "Took " << diff << " to reach " << pt << std::endl;
			lastTime = now;
		}
		
	private:
		
		Time lastTime;
		
	};

}

#endif
