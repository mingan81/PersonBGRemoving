/*
 * MotionDetector.h
 *
 *  Created on: 01/03/2015
 *      Author: Diana Gil (dianacgil@ieee.org)
 */

#ifndef MOTIONDETECTOR_H_
#define MOTIONDETECTOR_H_

#include "DipColorImage.h"
#include <opencv2/core/core.hpp>
#include <opencv2/video/background_segm.hpp>

/*!
* MotionDetector inherits from OpenCV BackgroundSubtractorMOG2
*/
class MotionDetector: public cv::BackgroundSubtractorMOG2
{
public:
	MotionDetector();
	~MotionDetector();

	void getMotionMask(const DipColorImage<uchar> &, DipColorImage<uchar>);

private:
	// Parameters for background subtraction
	// TODO: the first 3 are already there!!!
	const int BS_HISTORY;
	const float BS_THRESHOLD;
	const bool BS_SHADOW;
	const double LEARNING_RATE;
	const int MORPHO_SIZE;
	const int MORPHO_TYPE;
};

#endif /* MOTIONDETECTOR_H_ */
