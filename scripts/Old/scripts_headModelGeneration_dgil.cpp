// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/objdetect/objdetect.hpp>

// General
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <array>
#include <ctime>

// Impress...
#include "../include/humanmodel.hpp"

using namespace std;
using namespace cv;

// Function main
int main(int argc, char * argv[])
{
	// Capture the video
	if (argc != 2)
	{
		cerr << "Video filename unknown" << endl;
		exit(EXIT_FAILURE);
	}
	const string videoFilename = argv[1];
	VideoCapture videoObj(videoFilename);

	// Error in opening the video input
	if (!videoObj.isOpened())
	{
		cerr << "Unable to open video file: " << videoFilename << endl;
		exit(EXIT_FAILURE);
	}

	// Get frame range to detect head profiles
	
	array<int, 2> startEnd = getStartEnd( videoObj, "python/time.txt" );
	cout << startEnd[0] << " " << startEnd[1] << endl;

	clock_t t;
	t = clock();

	// Detect key frames and save head masks
	HeadSequence extHeads;
	detectHeadProfiles( videoObj, startEnd, extHeads );

	t = clock() - t;
	double timeSecs = ((double)t)/CLOCKS_PER_SEC;
	cout << "Time detecting key frames = " << timeSecs << " s "
		<< "--> " << timeSecs/60 << " minutes" << endl;

	cout << "Frame numbers of head profiles = " << extHeads.headFrameNumbers.at(LEFT) << " & "
		 << extHeads.headFrameNumbers.at(RIGHT) << endl
		 << "Frame numbers of frontal views = " << extHeads.headFrameNumbers.at(FRONT1) << " & "
		 << extHeads.headFrameNumbers.at(FRONT2) << endl;

	// Sort profiles
	// Left profile, center view, right profile, center view
	//sortProfiles( originalHeads, extHeads.headFrameNumbers );

	// Create a vector of head images from a video
	std::vector<cv::Mat> originalHeads(startEnd[1] - startEnd[0] + 1);
	video2HeadImages( videoObj, startEnd, originalHeads );

	// Get heads without background
	std::vector<cv::Mat> finalHeads;
	getSegmentedHeads( originalHeads, extHeads.headImages, finalHeads );
	

	return 0;
}

