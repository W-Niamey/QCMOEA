#ifndef BASE_H
#define BASE_H
#include <string>
#include <vector>
using namespace std;

namespace Base
{
	template <typename T>
	struct Pair
	{
		int dim;
		T value;
	};

	template <typename T>
	struct PairGreater
	{
		bool operator () (Pair <T> a, Pair <T> b)
		{
			return a.value > b.value;
		}
	};

	template <typename T>
	struct PairLess {
		bool operator () (Pair <T> a, Pair <T> b)
		{
			return a.value < b.value;
		}
	};

	long GetElapsedProcessTime();//Return Process times in millionSeconds;


}
#endif


