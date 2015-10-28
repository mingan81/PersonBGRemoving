/*
 * MotionDetector.cpp
 *
 *  Created on: 01/03/2015
 *      Author: Diana Gil (dianacgil@ieee.org)
 */

#include "MotionDetector.h"

/*!
* MotionDetector inherits from OpenCV BackgroundSubtractorMOG2
*/
MotionDetector::MotionDetector() : BS_HISTORY(200.0), BS_THRESHOLD(16.0), BS_SHADOW(false), LEARNING_RATE(0.0), MORPHO_SIZE(4), MORPHO_TYPE(cv::MORPH_ELLIPSE)
{

	// Parameters for background subtraction:
	// As they are constants, use initialization lists to initialize member variables
	// TODO: verify if these values make sense

	// Create Background Subtractor with MOG2 approach
	cv::BackgroundSubtractorMOG2(BS_HISTORY, BS_THRESHOLD, BS_SHADOW);
}

MotionDetector::~MotionDetector()
{
}

/*!
* @brief Get a foreground/motion mask based on background subtraction from OpenCV
* @param[in] frame : original image
* @param[out] fgMask : the resulting motion mask
*/
void MotionDetector::getMotionMask(const DipColorImage<uchar> & frame, DipColorImage<uchar> fgMask)
{
	// Remove background and get foreground mask
	// PROBLEM HERE!!!
	this->operator()(frame, fgMask, LEARNING_RATE);

	// Apply opening to reduce noise
	cv::Mat element = cv::getStructuringElement( MORPHO_TYPE, cv::Size( 2*MORPHO_SIZE + 1, 2*MORPHO_SIZE + 1 ), cv::Point( MORPHO_SIZE, MORPHO_SIZE ) );
	cv::morphologyEx( fgMask, fgMask, 2, element ); // 2 indicates opening operation
}
