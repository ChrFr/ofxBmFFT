class BmFFT
{
private:

public:
	
	
	/**
	 * Performs a Fast Fourier Transformation with the given samples. Results are stored in given bandVolumes.
	 *
	 * @param magnitude	  <no idea why baumeister calls it magnitude, in fact it determines the number of samples, see parameter samples>
	 * @param samples	  the stereo samples to be analyzed,
	 *					  stereo channels are interleaved (even numbers left channel, uneven right channel),
	 *					  expected to contain 2 * 2^magnitude complex numbers.
	 * @param numBands    the number of bands, the spectrum will be divided into
	 * @param bandVolumes array of floats, where the volumes of the bands are stored in(between 0 an 1),
	 *					  stereo channels are interleaved (even numbers left channel, uneven right channel),
	 *					  length is numBands * 2
	 */
	static void getSpectrum(const int magnitude, const float * samples, const int numBands, float * bandVolumes);
};
