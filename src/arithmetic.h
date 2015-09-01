//  ofxBmFFT
//  logic.h
//  Purpose: Macros and inline functions for integral arithmetics. Note that all indices are interpreted as Euclidean modulo 64 values. Also note that most of the operations are defined as generic macros in order to permit function overloading.spectrum analysis 
//  slight changes to adapt functions to work with Visual Studio 2012
//  @author Sascha Baumeister, Christoph Franke
//  @version 0.8 01/09/15

#ifndef ARITHMETIC_H_
#define ARITHMETIC_H_ 1
#include <csignal>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "meta.h"

#ifndef M_TAUl
# define M_TAUl		6.283185307179586476925286766559005768L /* tau = 2*pi */
#endif
 
#ifndef M_SQRT1_2l     //WARNING: added by me (Christoph), value is an educated guess (i think constant means 1/sqrt(2))
# define M_SQRT1_2l		0.707106781186547524401
#endif

namespace emath {
	/**
	 * Returns the next random floating point value within range [0.0, 1.0].
	 */
	inline double frand() {
		return std::rand() / (double) RAND_MAX;
	}


	/**
	 * Returns the number of nanoseconds that elapsed in between the left and right
	 * clock count.
	 */
	inline int64_t interval (const clock_t left, const clock_t right) {
		return ((int64_t) right - (int64_t) left) * 1000000000ULL / CLOCKS_PER_SEC;
	}


	/**
	 * Returns the largest (closest to positive infinity) integral value that is less than or
	 * equal to the binary logarithm of the given operand. Note that this operation will raise
	 * SIGFPE if the given operand is zero.
	 */
	inline uint64_t floorLog2 (const uint64_t operand) {
		if (operand == 0) std::raise(SIGFPE);
		return (bitSize(uint64_t) - 1) - __lzcnt(operand); //WARNING: changed by me (Christoph) without deeper research 
	}


	/**
	 * Returns the smallest (closest to negative infinity) integral value that is greater than
	 * or equal to the binary logarithm of the given operand. Note that this operation will raise
	 * SIGFPE if the given operand is zero.
	 */
	inline uint64_t ceilLog2 (const uint64_t operand) {
		if (operand == 0) std::raise(SIGFPE);
		if (operand == 1) return 0;
		return (uint64_t)bitSize(uint64_t) - __lzcnt(operand - 1); //WARNING: changed by me (Christoph) without deeper research 
	}


	/**
	 * Returns the absolute value of the given operand.
	 */
	inline signed char abs (const signed char operand) {
		return std::abs((signed int) operand);
	}
	inline signed short abs (const signed short operand) {
		return std::abs((signed int) operand);
	}
	inline signed int abs (const signed int operand) {
		return std::abs(operand);
	}
	inline signed long abs (const signed long operand) {
		return std::abs(operand);
	}
	inline signed long long abs (const signed long long operand) {
		return std::abs(operand);
	}
	inline float abs (const float operand) {
		return std::abs(operand);
	}
	inline double abs (const double operand) {
		return std::abs(operand);
	}
	inline long double abs (const long double operand) {
		return std::abs(operand);
	}
	inline float abs (const std::complex<float> operand) {
		return std::abs(operand);
	}
	inline double abs (const std::complex<double> operand) {
		return std::abs(operand);
	}
	inline long double abs (const std::complex<long double> operand) {
		return std::abs(operand);
	}


	/**
	 * Returns the signum of the given operand, which is -1 for negative operands, +1 for
	 * positive operands, and 0 for zero. Note that floating point versions may not return
	 * zero for zero arguments, because floating-point values lack a neutral zero
	 * representation, featuring +0.0 and -0.0 instead.
	 */
	inline char signum (const signed char operand) {
		return (operand > 0) - (operand < 0);
	}
	inline short signum (const signed short operand) {
		return (operand > 0) - (operand < 0);
	}
	inline int signum (const signed int operand) {
		return (operand > 0) - (operand < 0);
	}
	inline long signum (const signed long operand) {
		return (operand > 0) - (operand < 0);
	}
	inline long long signum (const signed long long operand) {
		return (operand > 0) - (operand < 0);
	}


	/**
	 * Returns the Euclidean quotient of the given dividend and divisor. Note that
	 * the result is guaranteed to be positive. See "Division and Modulus for Computer
	 * Scientists", Daan Leijen, 2001. Note that this operation may raise SIGFPE if
	 * the given divisor is zero.
	 */
	inline signed char  div_e (const signed char dividend, const signed char divisor) {
		const signed char div_t = dividend / divisor, mod_t = dividend % divisor;
		return mod_t >= 0 ? div_t : div_t - signum(divisor);
	}
	inline signed short div_e (const signed short dividend, const signed short divisor) {
		const signed short div_t = dividend / divisor, mod_t = dividend % divisor;
		return mod_t >= 0 ? div_t : div_t - signum(divisor);
	}
	inline signed int div_e (const signed int dividend, const signed int divisor) {
		const signed int div_t = dividend / divisor, mod_t = dividend % divisor;
		return mod_t >= 0 ? div_t : div_t - signum(divisor);
	}
	inline signed long div_e (const signed long dividend, const signed long divisor) {
		const signed long div_t = dividend / divisor, mod_t = dividend % divisor;
		return mod_t >= 0 ? div_t : div_t - signum(divisor);
	}
	inline signed long long div_e (const signed long long dividend, const signed long long divisor) {
		const signed long long div_t = dividend / divisor, mod_t = dividend % divisor;
		return mod_t >= 0 ? div_t : div_t - signum(divisor);
	}


	/**
	 * Returns the Euclidean modulo of the given dividend and divisor. Note that
	 * the result is guaranteed to be positive. See "Division and Modulus for Computer
	 * Scientists", Daan Leijen, 2001. Note that this operation may raise SIGFPE if
	 * the given divisor is zero.
	 */
	inline signed long long mod_e (const signed long long dividend, const signed long long divisor) {
		const signed long long mod_t = dividend % divisor;
		return mod_t >= 0 ? mod_t : mod_t + emath::abs(divisor);
	}
	inline int64_t imod_el (const signed long dividend, const signed long divisor) {
		const signed long mod_t = dividend % divisor;
		return mod_t >= 0 ? mod_t : mod_t + emath::abs(divisor);
	}
	inline signed int imod_e (const signed int dividend, const signed int divisor) {
		const signed int mod_t = dividend % divisor;
		return mod_t >= 0 ? mod_t : mod_t + emath::abs(divisor);
	}
	inline signed short imod_es (const signed short dividend, const signed short divisor) {
		const signed short mod_t = dividend % divisor;
		return mod_t >= 0 ? mod_t : mod_t + emath::abs(divisor);
	}
	inline signed char  imod_eb (const signed char dividend, const signed char divisor) {
		const signed char mod_t = dividend % divisor;
		return mod_t >= 0 ? mod_t : mod_t + emath::abs(divisor);
	}


	/**
	 * Returns the number of pairs in the shuffle table for the given magnitude.
	 */
	inline uint64_t shuffleSize (uint8_t magnitude) {
		return ((1ULL << magnitude) - (1ULL << magnitude >> (magnitude >> 1))) >> 1;
	}


	/**
	 * Binary-tree style recursive operations are highly expensive in terms of performance,
	 * primarily because they prevent the inlining of code, but also because they induce significant
	 * call costs due to their <tt>O(N·log(N))</tt> computational effort class. This can be avoided
	 * by using <tt>perfect shuffle</tt> index tables which mimic the combinatorial pattern caused by
	 * binary divide&conquer recursion, therefore allowing much faster calculations. Performance can
	 * be enhanced even further by <tt>caching</tt> the resulting shuffle indices within memory, as
	 * such recursive operations sharing the same binary magnitude are often called multiple times
	 * within a single process.<br />
	 * <br />
	 * This operation returns a lookup table containing <tt>2^m - 2^(m-m/2) perfect shuffle</tt>
	 * index pairs for the given binary magnitude <tt>m</tt>, interpreting the magnitude as a
	 * modulo 64 value. The results are lazily cached, and the algorithm ensures thread safe
	 * initialization of this cache for a given magnitude. Note that the resulting table only contains
	 * shuffle index pairs featuring a left index that is smaller than the corresponding right index,
	 * i.e. omitting shuffle index pairs where both are equal.
	 */
	typedef struct { uint64_t left; uint64_t right; } shuffle_pair;
	const shuffle_pair* shuffleTable (uint8_t magnitude);


	/**
	 * Returns the sine table for the given magnitude. If necessary, the
	 * table is lazily initialized, containing (2 ^ (magnitude - 2) + 1
	 * elements after initialization. The sine caches will contain sine
	 * values for an angle range of [0°, 90°], including both extremes.
	 * This is sufficient to quickly lookup the sine and cosine values
	 * over the whole angle range; tangent and co-tangent values can be
	 * rather quickly calculated as well.
	 */
	const long double* sineTable (uint8_t magnitude);


	/**
	 * Returns <tt>log<sub>b</sub>(z)</tt>, i.e. the principal branch solution for the infinitely
	 * multi-valued real based logarithm of a complex operand. The multi-valued nature implies a
	 * discontinuity, which is located at negative real values of the operand z. The exact behavior
	 * of this operation depends on the nature of the base:<ul>
	 * <li>if the base equals the operand z, then the result is always one.</li>
	 * <li>if the base is 0, then the result is zero.</li>
	 * <li>if the base is +1, then the result is infinite or NaN.</li>
	 * <li>if the base is Euler's unit e, then see log(complex) for details.</li>
	 * <li>if the base is any other positive number, then the principal branch solution is
	 * the one whose imaginary part is guaranteed to be within range <i>[-τ/2/ln(b),+τ/2/ln(b)[</i>.
	 * The solutions for the other branches may be obtained by adding <i>i·τ·k/ln(b)</i> (branch
	 * index k ∈ ℤ) to the principal branch solution. Regardless of the branch chosen, the
	 * real part of the result is guaranteed to be positive.</li>
	 * <li>if the base is -1, then the solutions for the other branches may be obtained by adding
	 * {@code 2·k} (with branch index k ∈ ℤ) to the principal branch solution.</li>
	 * <li>if the base is any other negative number, then the solutions for the other branches are
	 * also obtainable, but the branch correction becomes more complicated.</li>
	 * </ul>
	 * When combining this operation with it's inverse, as in <i>log(b, exp(b,z)) = exp(b, log(b,z))
	 * = z</i>, choose<ul>
	 * <li>if the base is positive: <i>k = floor(½ + Im(z)/τ · ln(b))</i></li>
	 * <li>if the base is -1:<ul>
	 * <li><i>k = floor(½ · (Re(z)+1))</i> if Re(z) < -1</li>
	 * <li><i>k = ceil(½ · (Re(z)-1))</i> if Re(z) >= -1</li>
	 * </ul>
	 * </ul>
	 */
	template<typename T>
	std::complex<T> log (const T base, const std::complex<T> z) {
		if (base == (T) M_El) return std::log(z);
		if (base == z.real() && z.imag() == (T) 0) return (T) 1;
		return base < (T) 0
			? std::log(z) / std::complex<T>(std::log(-base), (T) (-0.5l * M_TAUl))
			: std::log(z) / std::log(base);
	}


	/**
	 * Returns <tt>log<sub>z<sub>1</sub></sub>(z<sub>2</sub>)</tt>, i.e. the principal branch
	 * solution for the infinitely multi-valued complex based logarithm of a complex operand.
	 * The multi-valued nature implies a discontinuity, which is located at negative real
	 * values of the operand <tt>z<sub>2</sub></tt>. The exact behavior of this operation
	 * depends on the nature of the base <tt>z<sub>1</sub></tt>:
	 * <ul>
	 * <li>if the base z<sub>1</sub> equals the operand <tt>z<sub>2</sub></tt>, then the result is
	 * always one.</li>
	 * <li>if the base z<sub>1</sub> is Euler's unit e, then see log(complex) for details.</li>
	 * <li>if the base z<sub>1</sub> is a real number, then see log(double, complex).</li>
	 * <li>if the base z<sub>1</sub> is a complex number, the behavior of this operation is
	 * generally similar to log(double, complex) when passing negative bases. However, branch
	 * determination and correction become much more complicated.</li>
	 * </ul>
	 */
	template<typename T>
	std::complex<T> log (const std::complex<T> z1, const std::complex<T> z2) {
		if (z1 == z2) return (T) 1;
		if (z1.imag() == (T) 0) return emath::log(z1.real(), z2);
		return std::complex<T>(std::log(abs(z2)), arg(z2)) / std::complex<T>(std::log(abs(z1)), arg(z1));
	}


	/**
	 * Returns <tt>b<sup>z</sup></tt>, i.e. a real base raised to the power of a complex exponent.
	 * The result of <tt>1<sup>z</sup></tt> and <tt>b<sup>0</sup></tt> is always one, by definition
	 * including the special case <tt>0<sup>0</sup></tt> . The result of <tt>0<sup>z</sup></tt> is
	 * zero unless {@code z} is zero. Note that this operation is infinitely cyclic except for
	 * special bases like zero, one or infinity, and special exponents like 0 or infinity. This
	 * behavior causes the logarithm (as the inverse of this operation) to become multi-valued in ℂ,
	 * similarly to the sine function's cyclic behavior causing the arc sine function to become
	 * multi-valued in ℝ. The exact cycle behavior depends on the nature of the base:
	 * <ul>
	 * <li>if the given base is strictly positive except {@code +1}, then the cycling is simply
	 * vertical, with <i>exp(b, z + i·τ·k / ln(b)) = exp(b, z)</i> for any integer k.</li>
	 * <li>if the given base is +1, then the cycling is simply horizontal, with
	 * <i>exp(b, z + 2·k) = exp(b, z)</i> for any integer k.</li>
	 * <li>if the given base is strictly negative except -1, then the cycling is obliquely rotated
	 * from the vertical by an angle between ]0, τ/2[, depending on the absolute of the base.</li>
	 * </ul>
	 */
	template<typename T>
	std::complex<T> exp (const T base, const std::complex<T> z) {
		if (base < (T) 0) return emath::exp(base, z);
		if (base == (T) 0) return z == (T) 0 ? (T) 1 : (T) 0;
		if (base == (T) 1) return (T) 1;
		if (base == (T) M_El) return std::exp(z);

		const T logBase = std::log(base);
		const T reLog = logBase * z.real();
		if (z.imag() == (T) 0) return std::exp(reLog);
		const T imLog = logBase * z.imag();
		return std::complex<T>(cos(imLog), sin(imLog)) * std::exp(reLog);
	}


	/**
	 * Returns <tt>z<sub>1</sub><sup>z<sub>2</sub></sup></tt>, i.e. the principal branch solution
	 * for a complex base raised to the power of a complex exponent, which is infinitely
	 * multi-valued in ℂ. The exact behavior of this operation depends on the nature of the operands:
	 * <ul>
	 * <li>If the base <tt>z<sub>1</sub></tt> is a positive real number, then see
	 * exp(double, complex) for details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is an integer, then see pow(complex, int) for
	 * details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is or approximates a unit fraction, then see
	 * root(complex, int) for details.</li>
	 * <li>If the exponent <tt>z<sub>2</sub></tt> is a real number, but neither an integer nor a
	 * unit fraction, then it is technically nevertheless a dyadic fraction in lowest terms because
	 * any IEEE floating point number is. Therefore, this operation behaves as if root(complex, int)
	 * and pow(complex, int)} had been combined, but without the guarantee that the principal branch
	 * solution returned is the one whose argument (positive or negative) is closest to zero.</li>
	 * <li>Otherwise this operation displays a more delicate behavior: If the base is set constant
	 * and the exponent is varied, then this operation shows a cyclic pattern similarly to
	 * pow(complex, int), but rotated by varying degrees. If the base is varied and the exponent set
	 * constant, then this operation displays a pattern similarly to root(complex, int), but rotated
	 * by varying degrees resulting in cycle or spiral patterns. If z<sub>1</sub> = z<sub>2</sub>
	 * or z<sub>1</sub> = -z<sub>2</sub>, then this operation displays a cyclic pattern similarly to
	 * pow(complex, int)}, but bended around the origin. If z<sub>1</sub> = 1/z<sub>2</sub>, then
	 * this operation displays a cyclic shamrock pattern.</li>
	 */
	template<typename T>
	std::complex<T> exp (const std::complex<T> z1, const std::complex<T> z2) {
		if ((z1.imag() == (T) 0) & (z1.real() >= (T) 0)) return emath::exp(z1.real(), z2);

		if (z2.imag() == (T) 0) {
			const T re = z2.real();
			if (re == (int32_t) re) return power(z1, (int32_t) re);
			const T ire = (T) 1 / re;
			if (ire == (int32_t) ire) return root(z1, (int32_t) ire);
		}

		return std::exp(std::log(z1) * z2);
	}


	/**
	 * Returns <tt>z<sup>1/n</sup></tt>, i.e. the principal branch solution for the finitely
	 * multi-valued integer root of a complex operand z. It's multi-valued nature implies a
	 * discontinuity, which is located at negative real values of the operand z, except for
	 * n=±1 where this operation is single-valued, and n=0 where no valid branches exist.
	 * The principal branch solution is the one whose argument (positive or negative) is
	 * closest to zero; it's real part is guaranteed to be positive for |n| > 1. The
	 * solutions for the other branches may be obtained by multiplying the principal branch
	 * solution with <i>e<sup>i·τ·k/n</sup></i>, using branch index k ∈ ℤ within range
	 * <i>[0, |n|[</i>. When combining this operation with it's inverse, as in
	 * <i>(z<sup>1/n</sup>)<sup>n</sup> = (z<sup>n</sup>)<sup>1/n</sup> = z</i>, use
	 * k = (floor(abs(z) · n/τ + ½) % |n|, adding |n| if k < 0.<br />
	 * <br />
	 * Note that fractional exponentiation may be performed by first using this operation to
	 * calculate the fraction's divisor root of z (including choice of the proper branch),
	 * and subsequently raising the result to the power of the fraction's dividend using
	 * pow(complex,int)}.
	 */
	template<typename T>
	std::complex<T> root (const std::complex<T> z, const int32_t n) {
		switch (n) {
			case -1: {
				return (T) 1 / z;
			}
			case 0: {
				return std::complex<T>(NAN, NAN);
			}
			case 1: {
				return z;
			}
			case 2: {
				const T abs_value = abs(z);
				const std::complex<T> root = std::complex<T>(sqrt((T) 0.5l * (abs_value + z.real())), sqrt((T) 0.5l * (abs_value - z.real())));
				return z.imag() < (T) 0 ? conj(root) : root;
			}
			default: {
				const T inverse = (T) 1 / (T) n;
				return std::complex<T>(cos(pow(abs(z), inverse)), sin(arg(z) * inverse));
			}
		}
	}


	/**
	 * Returns <tt>z<sup>n</sup></tt>, i.e. a complex base raised to the power of an integer
	 * exponent. The special case <tt>z<sup>0</sup></tt> calculates to 1 in all cases, even for
	 * zero, infinite and NaN values of z, as defined in pow(double,double). Note that this
	 * operation is finitely cyclic, with (z · e<sup>i·τ·k/n</sup>)<sup>n</sup> =
	 * z<sup>n</sup> for every integer k within range [0, |n|[. This cyclic behavior causes the root
	 * operation (as the inverse of this operation) to become multi-valued in ℂ, similarly to the
	 * sine function's cyclic behavior causing the arc sine function to become multi-valued in ℝ.
	 * Also note that fractional exponentiation may be performed by first using root(complex, int)
	 * to calculate the fraction's divisor root of z (including choice of the proper branch), and
	 * subsequently using this operation to raise the result to the power of the fraction's dividend.
	 */
	template<typename T>
	std::complex<T> power (const std::complex<T> z, const int32_t n) {
		switch (n) {
			case -1: {
				return (T) 1 / z;
			}
			case 0: {
				return (T) 1;
			}
			case 1: {
				return z;
			}
			case 2: {
				return z*z;
			}
			default: {
				const T absolute = pow(abs(z), n), argument = arg(z) * n;
				return std::complex<T>(absolute * cos(argument), absolute * sin(argument));
			}
		}
	}


	/**
	 * Stores the Fast Fourier Transform of the given values in said values. The given
	 * array is expected to contain 2^magnitude complex numbers.
	 */
	template<typename T>
	void fft (const unsigned char magnitude, std::complex<T>* values) {
		const uint64_t size = 1ULL << magnitude;

		// note: live-calculation of shuffle indices is always slow compared to
		// table approach, regardless how fast the rotate function actually is ...
		const uint64_t  shufflePairCount = shuffleSize(magnitude);
		const shuffle_pair* shufflePairs = shuffleTable(magnitude);
		for (uint64_t index = 0; index < shufflePairCount; index += 2) {
			const shuffle_pair pair = shufflePairs[index];
			const std::complex<T> swap = values[pair.left];
			values[pair.left]  = values[pair.right];
			values[pair.right] = swap;
		}

		const long double* sine_cache = sineTable(magnitude);
		for (uint32_t mag = magnitude; mag > 0; --mag) {
			for (int64_t offset = (1ULL << mag >> 1) - 1; offset >= 0; --offset) {
				const uint64_t angle = (uint64_t) offset << (magnitude - mag);
				const std::complex<T> unit = angle <= (size >> 2)
					? std::complex<T>((T) +sine_cache[(size >> 2) - angle], (T) +sine_cache[angle])
					: std::complex<T>((T) -sine_cache[angle - (size >> 2)], (T) +sine_cache[(size >> 1) - angle]);

				for (uint64_t left = offset, right = offset + (1ULL << mag >> 1); right < size; left += 1ULL << mag, right += 1ULL << mag) {
					const std::complex<T> twiddle = values[right] * unit;
					values[right] = values[left] - twiddle;
					values[left] += twiddle;
				}
			}
		}

		const T norm = sqrt(ldexp((T) 1, -((int) magnitude)));
		for (uint64_t index = 0; index < size; ++index) {
			values[index] *= norm;
		}
	}


	/**
	 * Stores the Fast Fourier Transform of the given values in said values. If inverse is
	 * true, an FFT is performed; if it is false, an iFFT is performed instead. The given
	 * array is expected to contain 2^magnitude complex numbers.
	 */
	template<typename T>
	void fft (const bool inverse, unsigned char magnitude, std::complex<T>* values) {
		const uint64_t size = 1ULL << (magnitude &= 0x3f);
		if (inverse) {
			for (uint64_t index = 1; index < size; ++index) values[index] = -values[index];
		}

		fft(magnitude, values);

		if (inverse) {
			for (uint64_t index = 1; index < size; ++index) values[index] = -values[index];
		}
	}


	/**
	 * If inverse is true, the given natural spectrum is untangled/unfolded using the following
	 * operations for each corresponding pair of spectrum entries (r,l,t ∈ ℂ):
	 * forEach(l,r): t = r<sup>*</sup>; r = i·(t-l)/√2; l = (t+l)/√2;
	 * <br />
	 * If inverse is false, the given unfolded spectrum is folded again using the following
	 * operations (again r,l,t ∈ ℂ):
	 * forEach(l,r): t = i·r; r = (l-t)<sup>*</sup>/√2; l = (l+t)/√2</tt>
	 * <br />
	 * The untangling/unfolding of a spectrum created by an iFFT operation is required in order
	 * to cleanly separate it's left and right halves; in unfolded state, a spectrum's left half
	 * is solely influenced by the real parts of the values that went into the iFFT, and the
	 * right half is solely influenced by the corresponding imaginary parts. Note that the
	 * spectrum entries at index 0 and N/2 represent special cases, as their real parts belong
	 * to the left half of the spectrum, and their imaginary parts belong to the right half;
	 * this operation does not change them in any way.
	 */
	template<typename T>
	void fold (const bool inverse, unsigned char magnitude, std::complex<T>* spectrum) {
		const uint64_t size = 1ULL << (magnitude &= 0x3f);
		const T norm = (T) M_SQRT1_2l;
		const std::complex<T> imaginary_unit = std::complex<T>((T) 0.0l, (T) 1.0l);
		if (inverse) {
			for (uint64_t asc = 1, desc = size - 1; asc < desc; ++asc, --desc) {
				spectrum[desc] = conj(spectrum[desc]);
				const std::complex<T> left = spectrum[desc], right = spectrum[asc];
				spectrum[asc]  = (left + right) * norm;
				spectrum[desc] = (left - right) * norm * imaginary_unit;
			}
		} else {
			for (uint64_t asc = 1, desc = size - 1; asc < desc; ++asc, --desc) {
				const std::complex<T> left = spectrum[asc], right = spectrum[desc] * imaginary_unit;
				spectrum[asc]  = (left + right) * norm;
				spectrum[desc] = (left - right) * norm;
				spectrum[desc] = conj(spectrum[desc]);
			}
		}
	}
}

#endif
