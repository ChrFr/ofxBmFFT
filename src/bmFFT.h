//  ofxBmFFT
//  bmFFT.h
//  Purpose: spectrum analysis with custom FFT (provided by Sascha Baumeister)
//
//  @author Christoph Franke
//  @version 0.8 01/09/15

#include "arithmetic.h"

class BmFFT
{
private:

public:	

	//  Performs a FFT of the given samples and analyses the spectrum
	//
	//  @param magnitude <no idea why baumeister calls it magnitude, in fact it determines the number of samples, see parameter samples>
	//  @param samples the stereo samples to be analyzed, stereo channels are interleaved (even numbers left channel, uneven right channel), expected to contain 2 * 2^magnitude complex numbers.
	//  @param bandVolumes  array of floats, where the volumes of the bands are stored in(between 0 an 1), 
	//		   stereo channels are interleaved (even numbers left channel, uneven right channel), 
	//		   length is       2 * (log2(numSamples / 2) - 1) * 3	        
	//                    nr channels       octaves        1/3 octave
	//		   (e.g. 42 for a buffersize of 512, 21 for each channel)
	static void getSpectrum(const int buffersize, const float * samples, float * bandVolumes);
};
