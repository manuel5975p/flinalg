#ifndef UTIL
#define UTIL
#include <chrono>
namespace std{
	class stopwatch{
		private:
		unsigned long long nanos;
		public:
		stopwatch(){
			using namespace std::chrono;
			nanos = duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
		}
		unsigned long long elapsed(){
			using namespace std::chrono;
			return duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count() - nanos;
		}
		void reset(){
			using namespace std::chrono;
			nanos = duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
		}
	};
}
#endif