#include <cstddef>
#include <cstdio>
#include <cstring>
#include "arithmetic.h"
#include "logic.h"

using namespace std;
namespace emath {
	static const shuffle_pair* SHUFFLE_TABLES[bitSize(uint64_t) - 1];
	static const long double* SINE_TABLES[bitSize(uint64_t) - 1];


	const shuffle_pair* shuffleTable (uint8_t magnitude) {
		if (SHUFFLE_TABLES[magnitude &= 0x3f] == NULL) {
			const uint64_t size = shuffleSize(magnitude);
			shuffle_pair* table = new shuffle_pair[size];
			for (uint64_t offset = 0, stop = 1ULL << magnitude, left = 0; left < stop; ++left) {
				const uint64_t right = reverse(left) >> (bitSize(uint64_t) - magnitude);
				if (left < right) {
					shuffle_pair pair;
					pair.left = left;
					pair.right = right;
					table[offset++] = pair; //Christoph: cast didnt work, init struct pair added
				}
			}
			SHUFFLE_TABLES[magnitude] = table;
		}
		return SHUFFLE_TABLES[magnitude];
	}


	const long double* sineTable (uint8_t magnitude) {
		if (SINE_TABLES[magnitude &= 0x3f] == NULL) {
			const double angleDelta = ldexp(M_TAUl, -magnitude);
			const uint64_t tableSize = (1ULL << magnitude >> 2) + 1;
			long double* table = new long double[tableSize];
			for (uint64_t index = 0; index < tableSize; ++index) {
				table[index] = sinl(angleDelta * index);
			}
			SINE_TABLES[magnitude] = table;
		}
		return SINE_TABLES[magnitude];
	}
}
