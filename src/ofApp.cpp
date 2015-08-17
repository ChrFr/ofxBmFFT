#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
	ofEnableSmoothing();
	ofSetFrameRate(30);
	ofBackground(0);

	
	soundStream.listDevices();

	//if you want to set a different device id, default is microphone-in i guess (at least it is on all of my PCs)?
	//soundStream.setDeviceID(1);  

	// 0 output channels, 
	// 2 input channels
	// 44100 samples per second
	// 4 num buffers (latency)
	soundStream.setup(this, 0, 2, 44100, BUFFER_SIZE, 4);

	samples = new float[BUFFER_SIZE * 2];
	bandVolumes = new float[N_BANDS * 2];

	ofSetColor(0x666666);
}


//--------------------------------------------------------------
void ofApp::update(){
	ofBackground(80, 80, 20);
}

//--------------------------------------------------------------
void ofApp::draw(){


	ofNoFill();

	ofPushStyle();
	ofPushMatrix();
	ofTranslate(200, 150, 0);

	ofSetColor(245, 58, 135);
	ofSetLineWidth(3);
	
	// draw the left channel:
	ofBeginShape();
	for (unsigned int i = 0; i < BUFFER_SIZE / 2; i++){
		ofVertex(50 + i * 2, 100 - samples[i * 2] * 180.0f);
	}	
	ofEndShape(false);
	
	// draw the right channel:
	ofBeginShape();
	for (unsigned int i = 0; i < BUFFER_SIZE / 2; i++){
		ofVertex(650 + i * 2, 100 - samples[i * 2 + 1] * 180.0f);
	}	
	ofEndShape(false);

	// do the FFT and compute spectrum
	BmFFT::getSpectrum(MAGNITUDE, samples, N_BANDS, bandVolumes);


	// draw the FFT
	for (int i = 0; i < N_BANDS; i++){
		ofLine(50 + (i * 3), 600, 50 + (i * 3), 600 - bandVolumes[i * 2] * 300); //left
		ofLine(650 + (i * 3), 600, 650 + (i * 3), 600 - bandVolumes[i * 2 + 1] * 300); //right
	}
	
	ofPopMatrix();
	ofPopStyle();
}

void ofApp::audioReceived(float * input, int bufferSize, int nChannels){
	for (int i = 0; i < bufferSize * 2; i++){
		samples[i] = input[i];
	}
	bufferCounter++;
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
