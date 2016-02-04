//  ofxBmFFT
//  bmFFT.h
//  Purpose: spectrum analysis with custom FFT (provided by Sascha Baumeister)
//
//  @author Christoph Franke
//  @version 0.2 01/09/15

#include "arithmetic.h"

/// <summary>
/// contains functions to call the emath-library and build audio-spectrums
/// </summary>
class BmFFT {
	private:

	public:

		/// <summary>
		/// Performs an FFT of the given samples and creates a the spectrum with logarithmic averages (stored in bandVolumes)
		/// </summary>
		/// <param name='numSamples'>determines the number of samples</param>
		/// <param name='samples'> the stereo samples to be analyzed, stereo channels are interleaved (even numbers left channel, uneven right channel), expected to contain 2 * buffersize complex numbers</param>
		/// <param name='bandVolumes'>array of floats, where the volumes of the bands are stored in (between 0 an 1), length is 2 * log2(buffersize)</param>
		static void getSimpleSpectrum(const int numSamples, const float * samples, float * bandVolumes);

		/// <summary>
		/// Performs an FFT of the given samples and creates a the spectrum with logarithmic averages (stored in bandVolumes)
		/// </summary>
		/// <param name='numSamples'>determines the number of samples</param>
		/// <param name='samples'> the stereo samples to be analyzed, stereo channels are interleaved (even numbers left channel, uneven right channel), expected to contain 2 * buffersize complex numbers</param>
		/// <param name='bandVolumes'>bandVolumes  array of floats, where the volumes of the bands are stored in(between 0 an 1), stereo channels are interleaved (even numbers left channel, uneven right channel), length is 2 * (log2(numSamples/2) - 1) * 3 (e.g. 42 for a buffersize of 512, 21 for each channel)</param>
		static void getSpectrum(const int numSamples, const float * samples, float * bandVolumes);
};
