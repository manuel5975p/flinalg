#ifndef UTIL
#define UTIL
#include <chrono>
namespace std{
	/**
	@brief Class for time measurement
	*/
	class stopwatch{
		private:
		unsigned long long nanos;
		public:
		/**
		@brief Saves current time in nanoseconds
		*/
		stopwatch(){
			using namespace std::chrono;
			nanos = duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
		}
		/**
		@brief Returns the elapsed time in nanoseconds
		*/
		unsigned long long elapsed(){
			using namespace std::chrono;
			return duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count() - nanos;
		}
		/**
		@brief Resets the stopwatch
		*/
		void reset(){
			using namespace std::chrono;
			nanos = duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
		}
	};
}
#endif