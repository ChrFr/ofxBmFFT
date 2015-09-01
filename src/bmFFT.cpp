#include "bmFFT.h"

void BmFFT::getSpectrum(const int magnitude, const float * samples, const int numBands, float * bandVolumes){
	int numSamples = pow(2, magnitude); 

	std::complex<float>* complexSamples = new std::complex<float>[numSamples];

	for (int i = 0; i < numSamples; i++){
		complexSamples[i] = std::complex<float>(samples[i * 2], samples[i * 2 + 1]);
	}
	
	// do the fft
	emath::fft(false, magnitude, complexSamples);
	// seperate left and right channel
	emath::fold(false, magnitude, complexSamples);
	
	float * leftMagnitude = new float[numSamples / 2];
	float * rightMagnitude = new float[numSamples / 2];
	float out_real, out_img, power;

	//first element of each half contains unnessecary information
	// ascending first half
	for (int i = 1; i < numSamples / 2; i++) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		leftMagnitude[i] = 2.0*sqrt(power);
		//phase[i] = atan2(out_img, out_real);
	}

	// ascending 2nd half
	for (int i = numSamples - 1; i > numSamples / 2; i--) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		rightMagnitude[numSamples - i] = 2.0*sqrt(power);
		//phase[i] = atan2(out_img, out_real);
	}


	int step = (10 * (float)(numSamples / 2) / (float)numBands + 5) / 10;
	float leftSum = 0;
	float rightSum = 0;
	int counter = 0;
	int band = 0;
	for (int i = 1; i < numSamples / 2; i++){
		counter ++;
		leftSum += leftMagnitude[i];
		rightSum += rightMagnitude[i];
		if(counter == step){
			bandVolumes[band * 2] = leftSum / counter;
			bandVolumes[band * 2 + 1] = rightSum / counter;
			counter = 0;
			leftSum = 0;
			rightSum = 0;
			band++;
		}
	}

	
	/*
	int half = numSamples / 2;
	int counter = 0;
	int bandCount = 0;
	float leftSum = 0;
	float rightSum = 0;
	for (int i = numSamples; i > 0; i--){
		counter++;
		leftSum += leftMagnitude[i];
		rightSum += rightMagnitude[i];
		if(i == half){
			bandVolumes[numBands  - bandCount * 2] = leftSum / counter;
			bandVolumes[numBands  - bandCount - 1] = rightSum / counter;
			half /= 2;
			bandCount++;
			counter = 0;
			leftSum = 0;
			rightSum = 0;
		}

		if(bandCount == numBands)
			break;
	}*/

	delete[]complexSamples;
	delete[]leftMagnitude;
	delete[]rightMagnitude;

}
