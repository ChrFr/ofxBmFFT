//  ofxBmFFT
//  logic.h
//  Purpose: Inline functions for bit-set logic. Note that all indices are interpreted as Euclidean modulo 64 values.
//
//  slight changes to adapt functions to work with Visual Studio 2012
//  @author Sascha Baumeister, Christoph Franke
//  @version 0.8 01/09/15

#ifndef LOGIC_H_
#define LOGIC_H_ 1
#include <type_traits>
#include "meta.h"

namespace emath {
	/**
	 * Returns the given word's bit content at the given index.
	 */
	template<typename T>
	inline T getBit (const T word, const uint8_t index) {
		return (word >> index) & 1;
	}


	/**
	 * Flips the given word's bit content at the given index.
	 */
	template<typename T>
	inline void flipBit (T& word, const uint8_t index) {
		word ^= 1ULL << index;
	}


	/**
	 * Alters the given word by setting the bit at the given index to zero.
	 */
	template<typename T>
	inline void setBitOff (T& word, const uint8_t index) {
		word &= ~(1ULL << index);
	}


	/**
	 * Alters the given word by setting the bit at the given index to one.
	 */
	template<typename T>
	inline void setBitOn (T& word, const uint8_t index) {
		word |= 1ULL << index;
	}


	/**
	 * Alters the given word by setting the bit at the given index to the given value.
	 * While this operation is both more general than setBitOff() and setBitOn() and
	 * also branch-free, it requires around double as many operations.
	 */
	template<typename T>
	inline void setBit (T& word, const uint8_t index, const bool value) {
		word ^= (word ^ -(value & 1ULL)) & (1ULL << index);
	}


	/**
	 * Returns the number of one-bits in the given operand. Note: __cardinality_is and
	 * __cardinality_ib are solely required to ensure signed-to-unsigned casting without
	 * sign expansion.
	 */
	template<typename T>
	inline unsigned char cardinality (const T operand) {
		return __popcnt64((unsigned long long) operand);//ADAPTED TO VISUAL C++
	}
	template<>
	inline unsigned char cardinality (const unsigned long operand) {
		return __popcnt(operand);//ADAPTED TO VISUAL C++
	}
	template<>
	inline unsigned char cardinality (const unsigned int operand) {
		return __popcnt16(operand);//ADAPTED TO VISUAL C++
	}


	/**
	 * Returns the zero based index of the most significant bit set, or -1 if the
	 * given operand is zero. Note: msb_indexs and msb_indexb are solely required
	 * to ensure signed-to-unsigned casting without sign expansion.
	 */
	template<typename T>
	inline signed char msb_index (const T operand) {
		return (operand == 0) ? -1 : bitSize(unsigned long long) - 1 - __builtin_clzll((unsigned long long) operand);
	}
	template<>
	inline signed char msb_index (const unsigned long operand) {
		return (operand == 0) ? -1 : bitSize(unsigned long) - 1 - __lzcnt(operand); //ADAPTED TO VISUAL C++

	template<>
	inline signed char msb_index (const unsigned int operand) {
		return (operand == 0) ? -1 : bitSize(unsigned int) - 1 - __lzcnt(operand); //ADAPTED TO VISUAL C++
	}


	/**
	 * Returns the zero based index of the least significant bit set, or -1 if the
	 * given operand is zero. Note: lsb_indexs and lsb_indexb are solely required
	 * to ensure signed-to-unsigned casting without sign expansion.
	 */
	template<typename T>
	inline signed char lsb_index (const T operand) {
		return (operand == 0) ? -1 : __builtin_ctzll((unsigned long long) operand);  
	}
	template<>
	inline signed char lsb_index (const unsigned long operand) {
		return (operand == 0) ? -1 : __lzcnt(operand); //ADAPTED TO VISUAL C++
	}
	template<>
	inline signed char lsb_index (const unsigned int operand) {
		return (operand == 0) ? -1 : __lzcnt(operand); //ADAPTED TO VISUAL C++
	}


	/**
	 * Reverses all bits in the given value, combining Knuth's superscalar-friendly bit-parallel
	 * approach with processor-supported byte swap. 8bit algorithm devised by Sean Anderson,
	 * July 13, 2001.
	 */
	inline uint64_t reverse (uint64_t value) {
		value = ((value & 0x5555555555555555ULL) << 1) | ((value & 0xaaaaaaaaaaaaaaaaULL) >> 1);
		value = ((value & 0x3333333333333333ULL) << 2) | ((value & 0xccccccccccccccccULL) >> 2);
		value = ((value & 0x0f0f0f0f0f0f0f0fULL) << 4) | ((value & 0xf0f0f0f0f0f0f0f0ULL) >> 4);
		return _byteswap_uint64(value); //ADAPTED TO VISUAL C++,
	}
	inline uint32_t reverse (uint32_t value) {
		value = ((value & 0x55555555U) << 1) | ((value & 0xaaaaaaaaU) >> 1);
		value = ((value & 0x33333333U) << 2) | ((value & 0xccccccccU) >> 2);
		value = ((value & 0x0f0f0f0fU) << 4) | ((value & 0xf0f0f0f0U) >> 4);
		return _byteswap_ulong(value); //ADAPTED TO VISUAL C++,
	}
	inline uint16_t reverse (uint16_t value) {
		value = ((value & (uint16_t) 0x5555U) << 1) | ((value & (uint16_t) 0xaaaaU) >> 1);
		value = ((value & (uint16_t) 0x3333U) << 2) | ((value & (uint16_t) 0xccccU) >> 2);
		value = ((value & (uint16_t) 0x0f0fU) << 4) | ((value & (uint16_t) 0xf0f0U) >> 4);
		return _byteswap_ushort(value);  //ADAPTED TO VISUAL C++,
	}
	inline uint8_t reverse (const uint8_t value) {
		return ((value * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> bitSize(uint32_t);
	}
}

#endif
