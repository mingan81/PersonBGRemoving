

#include "../include/Core.hpp"

#define DBIDEF

// Declaration of extern dbi declared in Utils.hpp
DebugInfo dbi;

using namespace std;

const double HUGE_FPS = 500.0;
const double DEFAULT_FPS = 30.0;

// CONSTRUCTORS
// ---------------------
DipVideo::DipVideo() : cv::VideoCapture()
{
}

DipVideo::DipVideo(const string& videoFilename) : cv::VideoCapture(videoFilename)
{
	filename.assign(videoFilename);
}

DipVideo::DipVideo(int device) : cv::VideoCapture(device)
{
	camera = device;
}

// DESTRUCTOR
// ---------------------
DipVideo::~DipVideo()
{
	cv::VideoCapture::release();
}

// METHODS
// ---------------------
string DipVideo::getFilename()
{
	return filename;
}

int DipVideo::getCamera()
{
	return camera;
}

int DipVideo::getFrameWidth()
{
	return cv::VideoCapture::get(CV_CAP_PROP_FRAME_WIDTH);
}

int DipVideo::getFrameHeight()
{
	return cv::VideoCapture::get(CV_CAP_PROP_FRAME_HEIGHT);
}

int DipVideo::getFPS()
{
	double fps = cv::VideoCapture::get(CV_CAP_PROP_FPS);
	if (fps > HUGE_FPS)
	{
		cout << "**** Frame rate " << fps << " not valid ****" << endl
			 << "Default value = " << DEFAULT_FPS << endl;
			 
		fps = DEFAULT_FPS;
	}
	return fps;
}

int DipVideo::getNumFrames()
{
	return cv::VideoCapture::get(CV_CAP_PROP_FRAME_COUNT);
}

void DipVideo::checkInit()
{
	if (!cv::VideoCapture::isOpened())
	{
		cerr << "Unable to open video file: " << filename << endl;
		exit(EXIT_FAILURE);
	}
}

void DipVideo::restartVideo()
{
	// Restart the video, just in case
	cv::VideoCapture::set(CV_CAP_PROP_POS_AVI_RATIO, 0);
}

/*!
* @brief Store the video as a vector of frames
* @param[out] vecFrames : vector of original images
*/
void DipVideo::decodeFrames(std::vector<DipColorImage<uchar> > & vecFrames)
{
	this->restartVideo();
	double width  = this->getFrameWidth();
	double height  = this->getFrameHeight();

	DipColorImage<uchar> frame;

	for (;;)
	{
		this->read(frame);

		// There is a special behavior with Mat variables
		// It is necessary to have a separate image and copy it to the vector
		DipColorImage<uchar> tempFrame(width,height);
		tempFrame.init(0);

		// -- Save original frame in a vector
		if (!frame.empty())
		{
			//~ frame.copyTo(tempFrame);
			cvtColor(frame,tempFrame,CV_BGR2RGB);
			vecFrames.push_back(tempFrame);
		}
		else
		{
			// End of the video file
			break;
		}
	}
}
