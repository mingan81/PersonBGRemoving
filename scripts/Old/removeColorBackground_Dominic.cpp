///////////////////////////////////////////////
//File name: removeColorBackground.cpp
//Author: Dominic de Lanauze
//Date: 2014-10-23
//Description:	The function removeColorBackground takes a color portrait image of a
//				person and removes a uniform background by removing HSV pixel values within the range 
//				of the uniform background pixels.
///////////////////////////////////////////////

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

using namespace cv;
using namespace std;

//function prototypes
void removeColorBackground(Mat& colorImageInput, Mat& binaryFaceMask);


// Main function to test the function
int main(int argc, char** argv )
{
	Mat colorImageInput, binaryFaceMask;

	// ***PLEASE SPECIFY THE PATH OF THE IMAGE TO BE USED***
	colorImageInput = imread("C:\\Users\\Do\\Documents\\ImpressView\\programming\\impressView\\Dominic.jpg", CV_LOAD_IMAGE_COLOR);   // Read the file
	
    // Try to read the image
	if(! colorImageInput.data )  // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }
	else
	{
		// call the remove background function
		//removeColorBackground(colorImageInput, binaryFaceMask);
		removeUniformColorBackground_D3(colorImageInput, binaryFaceMask);
		//Display initial Image
		namedWindow( "Input Color Image", WINDOW_NORMAL );// Create a window for display.
		imshow( "Input Color Image", colorImageInput );  

		//Display result Mask
		namedWindow( "Result Mask", WINDOW_NORMAL );// Create a window for display.
		imshow( "Result Mask", binaryFaceMask );  

		waitKey(0); // Wait for a keystroke in the window
		return 0;
	}	
}
/*
void removeColorBackground(Mat& colorImageInput, Mat& binaryFaceMask)
{
	//Variable Declarations
	uchar meanRValue, meanGValue, meanBValue, RValueStdDev, GValueStdDev, BValueStdDev;
	uchar meanHValue, meanSValue, meanVValue;
	double HValueSum=0, SValueSum=0, VValueSum=0;
	uchar HValueStdDev, SValueStdDev, VValueStdDev;
	double HValueVariance, SValueVariance, VValueVariance;
	Vec3b refPixel1, refPixel2, refPixel3, refPixel4, BGRPixel;
	Vec3b HSVRefPixel1, HSVRefPixel2, HSVRefPixel3, HSVRefPixel4,HSVRefPixel5,HSVRefPixel6, HSVPixel;
	int cornerSideLength = 0, numberHSVValues=0;
	
	// Declare variable 2D array to store HSV background pixel reference values
	uchar** refPixelHSVValues;

	Mat blurredColorImageInput,HSVColorImageInput,HSVHistogram,HSVBlurredColorImageInput;
	int hBins=500,sBins=600;
	int channels[] = {0,1};
	int histSizes[] = {hBins,sBins};
	float hueRange[] = {0,179}, saturationRange[]={0,255};
	const float* ranges[] = {hueRange, saturationRange};
	float HFactor, SFactor, VFactor;
	int counter=0;
	Mat element;
	int kernelSize = 0;
	
	//Convert image to HSV format
	cvtColor(colorImageInput,HSVColorImageInput,COLOR_BGR2HSV);

	// Calculate Histogram for HSV image
	calcHist(&HSVColorImageInput, 1, channels, Mat(), HSVHistogram, 2, histSizes,ranges, true, false );
    normalize( HSVHistogram, HSVHistogram, 0, 255, NORM_MINMAX, -1, Mat() );

	// Extract pixels from two windows located in the top corners (2 corners) of the image (1/16 x 1/16) 
	if (HSVColorImageInput.rows>=HSVColorImageInput.cols)
		cornerSideLength = HSVColorImageInput.cols/16;
	else
		cornerSideLength = HSVColorImageInput.rows/16;

	numberHSVValues = cornerSideLength*cornerSideLength*2;

	refPixelHSVValues = new uchar*[numberHSVValues];
	for (int i = 0; i < numberHSVValues; i++)
	{
		refPixelHSVValues[i] = new uchar[3];
	}

	for (int x=0; x<cornerSideLength; x++)
	{
		for (int y=0; y<cornerSideLength; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y,x);
			refPixelHSVValues[counter][0] = HSVPixel.val[0];
			refPixelHSVValues[counter][1] = HSVPixel.val[1];
			refPixelHSVValues[counter][2] = HSVPixel.val[2];			
			counter++;
		}
	}
	
	for (int x=HSVColorImageInput.cols-1; x>HSVColorImageInput.cols-1-cornerSideLength; x--)
	{
		for (int y=0; y<cornerSideLength; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y,x);
			refPixelHSVValues[counter][0] = HSVPixel.val[0];
			refPixelHSVValues[counter][1] = HSVPixel.val[1];
			refPixelHSVValues[counter][2] = HSVPixel.val[2];			
			counter++;
		}
	}

	//Sum each values for mean calculation
	for (int i = 0; i < numberHSVValues; i++)
	{
		HValueSum = HValueSum + refPixelHSVValues[i][0];
		SValueSum = SValueSum + refPixelHSVValues[i][1];
		VValueSum = VValueSum + refPixelHSVValues[i][2];
	}

	// Calculate mean H,S and V values of background
	meanHValue = HValueSum/(numberHSVValues);
	meanSValue = SValueSum/(numberHSVValues);
	meanVValue = VValueSum/(numberHSVValues);

	//calculate variance of H,S,V values
	HValueVariance = 0;
	SValueVariance = 0;
	VValueVariance = 0;
	for (int i = 0; i < cornerSideLength*2; i++)
	{
		HValueVariance = HValueVariance + pow(double(refPixelHSVValues[i][0])-double(meanHValue),2.0);
		SValueVariance = SValueVariance + pow(double(refPixelHSVValues[i][1])-double(meanSValue),2.0);
		VValueVariance = VValueVariance + pow(double(refPixelHSVValues[i][2])-double(meanVValue),2.0);
	}
	HValueVariance = HValueVariance/(numberHSVValues);
	SValueVariance = SValueVariance/(numberHSVValues);
	VValueVariance = VValueVariance/(numberHSVValues);

	// Calculate standard deviations of H,S,V values of background
	HValueStdDev = sqrt(HValueVariance);
	if (HValueStdDev == 0)
		HValueStdDev = 1;
	SValueStdDev = sqrt(SValueVariance);
	if (SValueStdDev == 0)
		SValueStdDev = 1;
	VValueStdDev = sqrt(VValueVariance);
	if (VValueStdDev == 0)
		VValueStdDev = 1;

	//Convert source image to grayscale
	cvtColor(colorImageInput,binaryFaceMask,CV_RGB2GRAY);

	//Filtering factors (must be modified for each image!)
	if(HValueStdDev <2.0)
		HFactor = 18.0;
	else
		if(HValueStdDev <5.0)
			HFactor = 8.0;
		else
			HFactor = 4.0;
	if(SValueStdDev <2.0)
		SFactor = 22.0;
	else
		if(SValueStdDev <5.0)
			SFactor = 8.0;
		else
			SFactor = 7.0;
	
	VFactor = 1000.0;

	cout<<"Mean Hue Value: " << float(meanHValue)<<endl;
	cout<<"Mean Saturation Value: " << float(meanSValue)<<endl;
	cout<<"Mean V Value: " << float(meanVValue)<<endl;
	cout<<"Hue stddev: " << float(HValueStdDev)<<endl;
	cout<<"Saturation stddev: " << float(SValueStdDev)<<endl;
	cout<<"Values stddev: " << float(VValueStdDev)<<endl;

	//Generate binary mask of image with background removed
	for (int y=0;y<colorImageInput.rows;y++)
	{
		for (int x=0;x<colorImageInput.cols;x++)
		{
			BGRPixel = colorImageInput.at<Vec3b>(y,x);
			HSVPixel = HSVColorImageInput.at<Vec3b>(y,x);
			//if (abs(BGRPixel.val[0]-meanBValue)<(6*BValueStdDev) && abs(BGRPixel.val[1]-meanGValue)<(6*GValueStdDev) && abs(BGRPixel.val[2]-meanRValue)<(6*RValueStdDev))
			if (abs(HSVPixel.val[0]-meanHValue)<=(HFactor*HValueStdDev) && abs(HSVPixel.val[1]-meanSValue)<=(SFactor*SValueStdDev)) //&& abs(HSVPixel.val[2]-meanVValue)<=(VFactor*VValueStdDev))
			{
				binaryFaceMask.at<uchar>(y,x) = 0;
			}
			else
			{
				binaryFaceMask.at<uchar>(y,x) = 255;
			}
		}
	}
	//Display result Mask before morphological "closing"
	namedWindow( "Result Mask before morphological closing", WINDOW_NORMAL );// Create a window for display.
	imshow( "Result Mask before morphological closing", binaryFaceMask ); 

	//Perform morphological closing operation on mask
	kernelSize = double(binaryFaceMask.rows)*0.01;
	if (kernelSize%2 != 0 )
		kernelSize = kernelSize+1;
	
	//Display kernel size
	cout<<"Kernel Size: "<<kernelSize<<endl;
	
	//Create structuring element
	element = getStructuringElement(MORPH_ELLIPSE,Size(kernelSize,kernelSize),Point(-1,-1));
	
	//Perform morphological closing on mask
	morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_CLOSE, element ,Point(-1,-1),1,0,morphologyDefaultBorderValue());
}
*/