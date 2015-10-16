# ofxBmFFT

Adapted version of Sascha Baumeisters FFT. Wraps it for using in Openframeworks under Visual Studio 2012. 
Performs a spectrum analysis of the transformed signal.

## usage

download and copy into <your oF directory>/addons

### new project
can be included as an addon via the project-generator provided by openFrameworks

### add to existing project
* create filters \addons\ofxBmFFT\src in solution map
* add existing files from <your oF directory>\addons\ofxBmFFT\src\ (by rightklicking or via drag and drop into src-filter)
* add additional dependency ..\..\..\addons\ofxBmFFT\src\ (rightklick project, select properties->configuration->c/c++->general)

## example

the example performs an analysis of a line-in signal and renders it's stereo waveform and spectrum