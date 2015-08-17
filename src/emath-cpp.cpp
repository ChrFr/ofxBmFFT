/*
 ============================================================================
 Name        : emath-cpp.cpp
 Author      : Sascha Baumeister
 Version     :
 Copyright   : all rights reserved
 Description : Test program for emath cpp-library
 ============================================================================
 */

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "arithmetic.h"
#include "logic.h"

#define STRESS_MAGNITUDE 13
#define STRESS_LOOP_COUNT 5000

namespace emath {
	/**
	 * Raises a SIGTERM if the given value is zero.
	 */
	void assert (const bool value) {
		if (value == 0) std::raise(SIGTERM);
	}


	/**
	 * Performs the bit operation tests.
	 */
	void testEmathBitOperations () {
		uint64_t word = 0x01;
		assert( getBit(word,  0));
		// const int8_t left = -64, right = 64; // crashes in O1-O3!
		// assert( getBit(word, left));
		// assert( getBit(word, right));
		std::cout << "bit test 1a ok" << std::endl;

		setBitOn(word, 1);
		assert( getBit(word,  0));
		assert( getBit(word,  1));
		setBitOff(word, 1);
		assert( getBit(word,  0));
		assert(!getBit(word,  1));
		setBit(word, 1, 1);
		assert( getBit(word,  0));
		assert( getBit(word,  1));
		setBit(word, 1, 0);
		assert( getBit(word,  0));
		assert(!getBit(word,  1));
		std::cout << "bit test 1b ok" << std::endl;

		setBitOn(word, 63);
		assert( getBit(word,  0));
		assert( getBit(word, 63));
		setBitOff(word, 63);
		assert( getBit(word,  0));
		assert(!getBit(word, 63));
		setBit(word, 63, 1);
		assert( getBit(word,  0));
		assert( getBit(word, 63));
		std::cout << "bit test 1c ok" << std::endl;

		assert(cardinality((uint64_t) 0x0000000000000000ULL) ==  0);
		assert(cardinality((uint64_t) 0x0101010101010101ULL) ==  8);
		assert(cardinality((uint64_t) 0x0303030303030303ULL) == 16);
		assert(cardinality((uint64_t) 0x0707070707070707ULL) == 24);
		assert(cardinality((uint64_t) 0x0f0f0f0f0f0f0f0fULL) == 32);
		assert(cardinality((uint64_t) 0xf0f0f0f0f0f0f0f0ULL) == 32);
		assert(cardinality((uint64_t) 0xf1f1f1f1f1f1f1f1ULL) == 40);
		assert(cardinality((uint64_t) 0xf3f3f3f3f3f3f3f3ULL) == 48);
		assert(cardinality((uint64_t) 0xf7f7f7f7f7f7f7f7ULL) == 56);
		assert(cardinality((uint64_t) 0xffffffffffffffffULL) == 64);
		std::cout << "bit test 2a ok" << std::endl;

		assert(cardinality((uint32_t) 0x00000000UL) ==  0);
		assert(cardinality((uint32_t) 0x01010101UL) ==  4);
		assert(cardinality((uint32_t) 0x03030303UL) ==  8);
		assert(cardinality((uint32_t) 0x07070707UL) == 12);
		assert(cardinality((uint32_t) 0x0f0f0f0fUL) == 16);
		assert(cardinality((uint32_t) 0xf0f0f0f0UL) == 16);
		assert(cardinality((uint32_t) 0xf1f1f1f1UL) == 20);
		assert(cardinality((uint32_t) 0xf3f3f3f3UL) == 24);
		assert(cardinality((uint32_t) 0xf7f7f7f7UL) == 28);
		assert(cardinality((uint32_t) 0xffffffffUL) == 32);
		std::cout << "bit test 2b ok" << std::endl;

		assert(cardinality((uint16_t) 0x0000) ==  0);
		assert(cardinality((uint16_t) 0x0101) ==  2);
		assert(cardinality((uint16_t) 0x0303) ==  4);
		assert(cardinality((uint16_t) 0x0707) ==  6);
		assert(cardinality((uint16_t) 0x0f0f) ==  8);
		assert(cardinality((uint16_t) 0xf0f0) ==  8);
		assert(cardinality((uint16_t) 0xf1f1) == 10);
		assert(cardinality((uint16_t) 0xf3f3) == 12);
		assert(cardinality((uint16_t) 0xf7f7) == 14);
		assert(cardinality((uint16_t) 0xffff) == 16);
		std::cout << "bit test 2c ok" << std::endl;

		assert(cardinality((uint8_t) 0x00) ==  0);
		assert(cardinality((uint8_t) 0x01) ==  1);
		assert(cardinality((uint8_t) 0x03) ==  2);
		assert(cardinality((uint8_t) 0x07) ==  3);
		assert(cardinality((uint8_t) 0x0f) ==  4);
		assert(cardinality((uint8_t) 0xf0) ==  4);
		assert(cardinality((uint8_t) 0xf1) ==  5);
		assert(cardinality((uint8_t) 0xf3) ==  6);
		assert(cardinality((uint8_t) 0xf7) ==  7);
		assert(cardinality((uint8_t) 0xff) ==  8);
		std::cout << "bit test 2d ok" << std::endl;

		assert(msb_index((uint64_t) 0x0000000000000000ULL) == -1);
		assert(msb_index((uint64_t) 0x0000000000000001ULL) ==  0);
		assert(msb_index((uint64_t) 0x8000000000000000ULL) == 63);
		assert(msb_index((uint32_t) 0x00000000UL) == -1);
		assert(msb_index((uint32_t) 0x00000001UL) ==  0);
		assert(msb_index((uint32_t) 0x80000000UL) == 31);
		assert(msb_index((uint16_t) 0x0000) == -1);
		assert(msb_index((uint16_t) 0x0001) ==  0);
		assert(msb_index((uint16_t) 0x8000) == 15);
		assert(msb_index((uint8_t) 0x00) == -1);
		assert(msb_index((uint8_t) 0x01) ==  0);
		assert(msb_index((uint8_t) 0x80) ==  7);
		std::cout << "bit test 3a ok" << std::endl;

		assert(lsb_index((uint64_t) 0x0000000000000000ULL) == -1);
		assert(lsb_index((uint64_t) 0x0000000000000001ULL) ==  0);
		assert(lsb_index((uint64_t) 0x8000000000000000ULL) == 63);
		assert(lsb_index((uint32_t) 0x00000000UL) == -1);
		assert(lsb_index((uint32_t) 0x00000001UL) ==  0);
		assert(lsb_index((uint32_t) 0x80000000UL) == 31);
		assert(lsb_index((uint16_t) 0x0000) == -1);
		assert(lsb_index((uint16_t) 0x0001) ==  0);
		assert(lsb_index((uint16_t) 0x8000) == 15);
		assert(lsb_index((uint8_t) 0x00) == -1);
		assert(lsb_index((uint8_t) 0x01) ==  0);
		assert(lsb_index((uint8_t) 0x80) ==  7);
		std::cout << "bit test 3b ok" << std::endl;

		assert(reverse((uint64_t) 0x0000000000000000ULL) == 0x0000000000000000ULL);
		assert(reverse((uint64_t) 0x0000000000000001ULL) == 0x8000000000000000ULL);
		assert(reverse((uint64_t) 0x8000000000000000ULL) == 0x0000000000000001ULL);
		assert(reverse((uint64_t) 0x0800000000000010ULL) == 0x0800000000000010ULL);
		assert(reverse((uint32_t) 0x00000000U) == 0x00000000U);
		assert(reverse((uint32_t) 0x00000001U) == 0x80000000U);
		assert(reverse((uint32_t) 0x80000000U) == 0x00000001U);
		assert(reverse((uint32_t) 0x08000010U) == 0x08000010U);
		assert(reverse((uint16_t) 0x0000U) == 0x0000U);
		assert(reverse((uint16_t) 0x0001U) == 0x8000U);
		assert(reverse((uint16_t) 0x8000U) == 0x0001U);
		assert(reverse((uint16_t) 0x0810U) == 0x0810U);
		assert(reverse((uint8_t)  0x00U) == 0x00U);
		assert(reverse((uint8_t)  0x01U) == 0x80U);
		assert(reverse((uint8_t)  0x80U) == 0x01U);
		assert(reverse((uint8_t)  0x81U) == 0x81U);
		std::cout << "bit test 4 ok" << std::endl;
	}


	/**
	 * Performs the algebraic operation tests.
	 */
	void testEmathAlgebraicOperations () {
		assert(abs(0LL) == 0);
		assert(abs(+1LL) == +1);
		assert(abs(-1LL) == +1);
		assert(abs(+10LL) == +10);
		assert(abs(-10LL) == +10);
		assert(abs(+1000LL) == +1000);
		assert(abs(-1000LL) == +1000);
		assert(abs(+1000000LL) == +1000000);
		assert(abs(-1000000LL) == +1000000);
		assert(abs(+1000000000000LL) == +1000000000000LL);
		assert(abs(-1000000000000LL) == +1000000000000LL);
		assert(abs(0) == 0);
		assert(abs(+1) == +1);
		assert(abs(-1) == +1);
		assert(abs(+10) == +10);
		assert(abs(-10) == +10);
		assert(abs(+1000) == +1000);
		assert(abs(-1000) == +1000);
		assert(abs(+1000000) == +1000000);
		assert(abs(-1000000) == +1000000);
		assert(abs((int16_t) 0) == 0);
		assert(abs((int16_t) +1) == +1);
		assert(abs((int16_t) -1) == +1);
		assert(abs((int16_t) +10) == +10);
		assert(abs((int16_t) -10) == +10);
		assert(abs((int16_t) +1000) == +1000);
		assert(abs((int16_t) -1000) == +1000);
		assert(abs((int8_t) 0) == 0);
		assert(abs((int8_t) +1) == +1);
		assert(abs((int8_t) -1) == +1);
		assert(abs((int8_t) +10) == +10);
		assert(abs((int8_t) -10) == +10);
		assert(abs((int8_t) +100) == +100);
		assert(abs((int8_t) -100) == +100);
		std::cout << "algebra test 1 ok" << std::endl;

		assert(signum(0LL) == 0);
		assert(signum(+1LL) == +1);
		assert(signum(-1LL) == -1);
		assert(signum(+10LL) == +1);
		assert(signum(-10LL) == -1);
		assert(signum(+1000LL) == +1);
		assert(signum(-1000LL) == -1);
		assert(signum(+1000000LL) == +1);
		assert(signum(-1000000LL) == -1);
		assert(signum(+1000000000000LL) == +1);
		assert(signum(-1000000000000LL) == -1);
		assert(signum(0) == 0);
		assert(signum(+1) == +1);
		assert(signum(-1) == -1);
		assert(signum(+10) == +1);
		assert(signum(-10) == -1);
		assert(signum(+1000) == +1);
		assert(signum(-1000) == -1);
		assert(signum(+1000000) == +1);
		assert(signum(-1000000) == -1);
		assert(signum((int16_t) 0) == 0);
		assert(signum((int16_t) +1) == +1);
		assert(signum((int16_t) -1) == -1);
		assert(signum((int16_t) +10) == +1);
		assert(signum((int16_t) -10) == -1);
		assert(signum((int16_t) +1000) == +1);
		assert(signum((int16_t) -1000) == -1);
		assert(signum((int8_t) 0) == 0);
		assert(signum((int8_t) +1) == +1);
		assert(signum((int8_t) -1) == -1);
		assert(signum((int8_t) +10) == +1);
		assert(signum((int8_t) -10) == -1);
		assert(signum((int8_t) +100) == +1);
		assert(signum((int8_t) -100) == -1);
		std::cout << "algebra test 2 ok" << std::endl;

		const std::complex<long double> ld0 = std::complex<long double>(+2.0l, +3.0l);
		const std::complex<long double> ld1 = std::complex<long double>(+2.0l, -3.0l);
		const std::complex<long double> ld2 = std::complex<long double>(-2.0l, -3.0l);
		const std::complex<long double> ld3 = std::complex<long double>(-2.0l, +3.0l);
		const std::complex<double> d0 = std::complex<double>(+2.0, +3.0);
		const std::complex<double> d1 = std::complex<double>(+2.0, -3.0);
		const std::complex<double> d2 = std::complex<double>(-2.0, -3.0);
		const std::complex<double> d3 = std::complex<double>(-2.0, +3.0);
		const std::complex<float> f0 = std::complex<float>(+2.0f, +3.0f);
		const std::complex<float> f1 = std::complex<float>(+2.0f, -3.0f);
		const std::complex<float> f2 = std::complex<float>(-2.0f, -3.0f);
		const std::complex<float> f3 = std::complex<float>(-2.0f, +3.0f);
		assert(ld0 == signum(ld0) * abs(ld0));
		assert(ld1 == signum(ld1) * abs(ld1));
		assert(ld2 == signum(ld2) * abs(ld2));
		assert(ld3 == signum(ld3) * abs(ld3));
		assert(d0 == signum(d0) * abs(d0));
		assert(d1 == signum(d1) * abs(d1));
		assert(d2 == signum(d2) * abs(d2));
		assert(d3 == signum(d3) * abs(d3));
		assert(f0 == signum(f0) * abs(f0));
		assert(f1 == signum(f1) * abs(f1));
		assert(f2 == signum(f2) * abs(f2));
		assert(f3 == signum(f3) * abs(f3));
		std::cout << "algebra test 3 ok" << std::endl;

		assert(div_e((int64_t) +7, (int64_t) +3) == +2);
		assert(div_e((int64_t) +7, (int64_t) -3) == -2);
		assert(div_e((int64_t) -7, (int64_t) +3) == -3);
		assert(div_e((int64_t) -7, (int64_t) -3) == +3);
		assert(mod_e((int64_t) +7, (int64_t) +3) == +1);
		assert(mod_e((int64_t) +7, (int64_t) -3) == +1);
		assert(mod_e((int64_t) -7, (int64_t) +3) == +2);
		assert(mod_e((int64_t) -7, (int64_t) -3) == +2);
		assert(div_e((int32_t) +7, (int32_t) +3) == +2);
		assert(div_e((int32_t) +7, (int32_t) -3) == -2);
		assert(div_e((int32_t) -7, (int32_t) +3) == -3);
		assert(div_e((int32_t) -7, (int32_t) -3) == +3);
		assert(mod_e((int32_t) +7, (int32_t) +3) == +1);
		assert(mod_e((int32_t) +7, (int32_t) -3) == +1);
		assert(mod_e((int32_t) -7, (int32_t) +3) == +2);
		assert(mod_e((int32_t) -7, (int32_t) -3) == +2);
		assert(div_e((int16_t) +7, (int16_t) +3) == +2);
		assert(div_e((int16_t) +7, (int16_t) -3) == -2);
		assert(div_e((int16_t) -7, (int16_t) +3) == -3);
		assert(div_e((int16_t) -7, (int16_t) -3) == +3);
		assert(mod_e((int16_t) +7, (int16_t) +3) == +1);
		assert(mod_e((int16_t) +7, (int16_t) -3) == +1);
		assert(mod_e((int16_t) -7, (int16_t) +3) == +2);
		assert(mod_e((int16_t) -7, (int16_t) -3) == +2);
		assert(div_e((int8_t) +7, (int8_t) +3) == +2);
		assert(div_e((int8_t) +7, (int8_t) -3) == -2);
		assert(div_e((int8_t) -7, (int8_t) +3) == -3);
		assert(div_e((int8_t) -7, (int8_t) -3) == +3);
		assert(mod_e((int8_t) +7, (int8_t) +3) == +1);
		assert(mod_e((int8_t) +7, (int8_t) -3) == +1);
		assert(mod_e((int8_t) -7, (int8_t) +3) == +2);
		assert(mod_e((int8_t) -7, (int8_t) -3) == +2);
		std::cout << "algebra test 4 ok" << std::endl;

		assert(floorLog2(1) == 0);
		assert(floorLog2(2) == 1);
		assert(floorLog2(3) == 1);
		assert(floorLog2(4) == 2);
		assert(floorLog2(4) == 2);
		assert(ceilLog2(1) == 0);
		assert(ceilLog2(2) == 1);
		assert(ceilLog2(3) == 2);
		assert(ceilLog2(4) == 2);
		assert(ceilLog2(5) == 3);
		std::cout << "algebra test 5 ok" << std::endl;
	}


	/**
	 * Performs the table access tests.
	 */
	void testEmathTableOperations () {
		uint64_t  size;
		const shuffle_pair* shuffles;
		size = shuffleSize(0);
		shuffles = shuffleTable(0);
		assert(size == 0);
		assert(shuffles != NULL);

		size = shuffleSize(1);
		shuffles = shuffleTable(1);
		assert(size == 0);
		assert(shuffles != NULL);

		size = shuffleSize(2);
		shuffles = shuffleTable(2);
		assert(size == 1);
		assert(shuffles != NULL);
		assert(shuffles[0].left  == 1);
		assert(shuffles[0].right == 2);

		size = shuffleSize(3);
		shuffles = shuffleTable(3);
		assert(size == 2);
		assert(shuffles != NULL);
		assert(shuffles[0].left  == 1);
		assert(shuffles[0].right == 4);
		assert(shuffles[1].left  == 3);
		assert(shuffles[1].right == 6);
		shuffleTable(15);
		std::cout << "table test 1 ok" << std::endl;

		const long double* sines;
		sines = sineTable(0);
		assert(sines != NULL);
		assert(sines[0] == 0.0l);
		sines = sineTable(1);
		assert(sines != NULL);
		assert(sines[0] == 0.0l);
		sines = sineTable(2);
		assert(sines != NULL);
		assert(sines[0] == 0.0l);
		assert(sines[1] == 1.0l);
		sines = sineTable(3);
		assert(sines != NULL);
		assert(sines[0] == 0.0l);
		assert(sines[2] == 1.0l);
		sineTable(15);
		std::cout << "table test 2 ok" << std::endl;
	}

	/**
	 * Performs the FFT tests.
	 */
	void testEmathFFTOperations () {
		const uint8_t magnitude = 15;
		const uint64_t size = 1 << magnitude;

		// note that C++ causes problems with stack allocation of arrays with a few
		// thousand elements; an amazig "feature" of a language that lacks active
		// heap management. C does not share this trait ...
		std::complex<float>* fvalues = new std::complex<float>[size];
		array_fill(fvalues, size, 1.0f);
		std::complex<double>* dvalues = new std::complex<double>[size];
		array_fill(dvalues, size, 1.0);
		std::complex<long double>* lvalues = new std::complex<long double>[size];
		array_fill(lvalues, size, 1.0l);

		fft(false, magnitude, fvalues);
		fft(true, magnitude, fvalues);
		fft(false, magnitude, dvalues);
		fft(true, magnitude, dvalues);
		fft(false, magnitude, lvalues);
		fft(true, magnitude, lvalues);
		std::cout << "fft function test ok" << std::endl;
		delete[] fvalues;
		delete[] dvalues;
		delete[] lvalues;
	}



	/**
	 * Performs a complex double FFT stress test.
	 */
	void stressEmathFloatFFT () {
		const uint64_t values_length = (uint64_t) 1ULL << STRESS_MAGNITUDE;
		std::complex<float>* values = new std::complex<float>[values_length];
		for (uint64_t index = 0; index < values_length; ++index) {
			values[index] = std::complex<float>((float) frand() - 0.5f, (float) frand() - 0.5f);
		}

		// hot test
		const clock_t clock0 = clock();
		for (uint64_t loop = 0; loop < STRESS_LOOP_COUNT; ++loop) {
			fft(false, STRESS_MAGNITUDE, values);
			fft(true, STRESS_MAGNITUDE, values);
		}
		const clock_t clock1 = clock();

		const double hash = abs(values[0]);
		printf("complex<float> (2*32bit): %lldµs per FFT (checksum %f).\n", (long long int) interval(clock0, clock1) / (1000 * STRESS_LOOP_COUNT * 2), hash);
		delete[] values;
	}


	/**
	 * Performs a complex double FFT stress test.
	 */
	void stressEmathDoubleFFT () {
		const uint64_t values_length = (uint64_t) 1ULL << STRESS_MAGNITUDE;
		std::complex<double>* values = new std::complex<double>[values_length];
		for (uint64_t index = 0; index < values_length; ++index) {
			values[index] = std::complex<double>(frand() - 0.5, frand() - 0.5);
		}

		// hot test
		const clock_t clock0 = clock();
		for (int64_t loop = STRESS_LOOP_COUNT; loop > 0; --loop) {
			fft(false, STRESS_MAGNITUDE, values);
			fft(true, STRESS_MAGNITUDE, values);
		}
		const clock_t clock1 = clock();

		const double hash = abs(values[0]);
		printf("complex<double> (2*64bit): %lldµs per FFT (checksum %f).\n", (long long int) interval(clock0, clock1) / (1000 * STRESS_LOOP_COUNT * 2), hash);
		delete[] values;
	}


	/**
	 * Performs a complex long double FFT stress test.
	 */
	void stressEmathLongDoubleFFT () {
		const uint64_t values_length = (uint64_t) 1ULL << STRESS_MAGNITUDE;
		std::complex<long double>* values = new std::complex<long double>[values_length];
		for (uint64_t index = 0; index < values_length; ++index) {
			values[index] = std::complex<long double>(frand() - 0.5l, frand() - 0.5l);
		}

		// hot test
		const clock_t clock0 = clock();
		for (int64_t loop = STRESS_LOOP_COUNT; loop > 0; --loop) {
			fft(false, STRESS_MAGNITUDE, values);
			fft(true, STRESS_MAGNITUDE, values);
		}
		const clock_t clock1 = clock();

		const double hash = abs(values[0]);
		printf("complex<long double> (2*80bit): %lldµs per FFT (checksum %f).\n", (long long int) interval(clock0, clock1) / (1000 * STRESS_LOOP_COUNT * 2), hash);
		delete[] values;
	}
}

/**
 * Test application.
 */
int main() {
	// initialize values
	srand(time(NULL));

	emath::testEmathBitOperations();
	emath::testEmathAlgebraicOperations();
	emath::testEmathTableOperations();
	emath::testEmathFFTOperations();
	emath::stressEmathFloatFFT();
	emath::stressEmathDoubleFFT();
	emath::stressEmathLongDoubleFFT();

	std::cout << "tests ok!" << std::endl;
	return EXIT_SUCCESS;
}
