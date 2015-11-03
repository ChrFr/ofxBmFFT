#include "bmFFT.h"
#include <iostream>

double log2(double x){
	return log(x)/log(2);
}

void BmFFT::getSimpleSpectrum(const int numSamples, const float * samples, float * bandVolumes){
	int magnitude = log2(numSamples / 2);

	std::complex<float>* complexSamples = new std::complex<float>[numSamples / 2];

	for (int i = 0; i < numSamples / 2; i++){
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
	//ascending first half
	for (int i = 1; i < numSamples / 2; i++) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		leftMagnitude[i - 1] = 2.0*sqrt(power);
	}

	// ascending 2nd half
	for (int i = numSamples - 1; i > numSamples / 2; i--) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		rightMagnitude[numSamples - i - 1] = 2.0*sqrt(power);
	}

	int half = numSamples / 4;
	int counter = 0;
	int band = magnitude - 1;
	float leftSum = 0;
	float rightSum = 0;
	for (int i = numSamples / 2; i > 0; i--){
		leftSum += leftMagnitude[i];
		rightSum += rightMagnitude[i];
		counter++;
		if(i == half){
			bandVolumes[band * 2] = leftSum / counter;
			bandVolumes[band * 2 + 1] = rightSum / counter;
			half /= 2;
			band--;
			leftSum = 0;
			rightSum = 0;
			counter = 0;
		}
	}
	/*
	int numOct = log2(numSamples / 2) - 1;
	int half = (numSamples / 2 - 1) / 2;
	// octaves
	for (int i = numSamples / 2 - 1; i >= 0; i--){
		int half = high / 2;		
		// special case: next to last octave only has 2 values-> take last octave(single value) as last third (-> last band is duplicated)
		if(high - low == 2){
			low--;
		}
		float leftSum = 0;
		float rightSum = 0;
		int step = (high - low) / 3;
		int counter = 0;
		int oBand = 0;
		// divide octave into equal thirds
		for(int k = high; k > low; k--){
			counter ++;			
			leftSum += leftMagnitude[k];
			rightSum += rightMagnitude[k];
			if(counter >= step){
				int band = 3 * i - oBand;
				// interchanging
				bandVolumes[band * 2] = leftSum / counter;
				bandVolumes[band * 2 + 1] = rightSum / counter;
				counter = 0;
				leftSum = 0;
				rightSum = 0;
				oBand++;
			}
			if(oBand >= 3)
				break;
		}
		high /= 2;
	}*/

	delete[]complexSamples;
	delete[]leftMagnitude;
	delete[]rightMagnitude;

}

void BmFFT::getSpectrum(const int numSamples, const float * samples, float * bandVolumes){
	int magnitude = log2(numSamples / 2);

	std::complex<float>* complexSamples = new std::complex<float>[numSamples / 2];

	for (int i = 0; i < numSamples / 2; i++){
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
	//ascending first half
	for (int i = 1; i < numSamples / 2; i++) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		leftMagnitude[i - 1] = 2.0*sqrt(power);
	}

	// ascending 2nd half
	for (int i = numSamples - 1; i > numSamples / 2; i--) {
		out_real = complexSamples[i].real();
		out_img = complexSamples[i].imag();

		power = out_real * out_real + out_img * out_img;
		rightMagnitude[numSamples - i - 1] = 2.0*sqrt(power);
	}
	
	int numOct = log2(numSamples / 2) - 1;
	int high = numSamples / 2 - 1;
	// octaves
	for (int i = numOct; i >= 0; i--){
		int low = high / 2;		
		// special case: next to last octave only has 2 values-> take last octave(single value) as last third (-> last band is duplicated)
		if(high - low == 2){
			low--;
		}
		float leftSum = 0;
		float rightSum = 0;
		int step = (high - low) / 3;
		int counter = 0;
		int oBand = 0;
		// divide octave into equal thirds
		for(int k = high; k > low; k--){
			counter ++;			
			leftSum += leftMagnitude[k];
			rightSum += rightMagnitude[k];
			if(counter >= step){
				int band = 3 * i - oBand;
				// interchanging
				bandVolumes[band * 2] = leftSum / counter;
				bandVolumes[band * 2 + 1] = rightSum / counter;
				counter = 0;
				leftSum = 0;
				rightSum = 0;
				oBand++;
			}
			if(oBand >= 3)
				break;
		}
		high /= 2;
	}

	delete[]complexSamples;
	delete[]leftMagnitude;
	delete[]rightMagnitude;

}
