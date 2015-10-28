#ifndef HUMANMODEL_HPP
#define HUMANMODEL_HPP

#include <opencv2/video/background_segm.hpp>
#include <array>
using namespace cv;

enum ZoneInfoType{ Width, Height, Val_min, Val_max, Dynamic, Average, Area, Volume, StdDev, NbPtNoZ };

enum headViews_t {LEFT, FRONT1, RIGHT, FRONT2};


template<typename Tin, typename Tout> int IsMinOrMax(const Tin valeur, Tout & min, Tout & max)
{
	if (valeur < min) min = Tout(valeur);
	if (valeur > max) max = Tout(valeur);
	return 0;
}


 



 // Result of detecting head profiles
 struct HeadSequence
 {
	 std::vector<cv::Mat> headImages;
	 std::vector<int> headFrameNumbers;
 };



cv::Rect detectFace(cv::Mat&, cv::Mat&);
cv::Rect detectFace_dgil(cv::Mat& frame,cv::Mat& crop);
cv::Rect detectNose(cv::Mat&,cv::Mat&);
void GetVerticalEdge(cv::Mat&,cv::Mat&);
void FindBlobs(const cv::Mat&, std::vector < std::vector<cv::Point2i> >&);
std::vector<int> CountBrightPixels(cv::Mat, int, int);
std::vector<int> CountBrightPixelsCol(cv::Mat&);
cv::Mat HeadForegroundExtraction(cv::Mat&);
int ImLabel_QG(cv::Mat& , int, cv::Mat& );
std::array<int, 2> getStartEnd(cv::VideoCapture&, const std::array<int, 2> );
void video2HeadImages( cv::VideoCapture&, const std::array<int, 2>&, std::vector<cv::Mat>& );
cv::Ptr<BackgroundSubtractor> createBgSubtractor();
cv::Mat removeBackground( const cv::Mat&, cv::Ptr<BackgroundSubtractor> );
cv::Rect enlargeBBox( const cv::Rect&, double, int );
std::vector<int> getPosProfiles( const vector<int>&, const vector<int>&, int );
cv::Mat getHeadImg( const cv::Mat&, const cv::Mat&, const int&, string );
bool detectEar( cv::Mat&, std::string );
int identifyProfiles( cv::Mat&, cv::Mat&, std::vector<int>& );
void detectHeadProfiles( cv::VideoCapture&, const std::array<int, 2>&, HeadSequence& );
void sortProfiles( const std::vector<cv::Mat>&, std::vector<int>& );
void getSegmentedHeads( const std::vector<cv::Mat>&, const std::vector<cv::Mat>&, std::vector<cv::Mat>& );
void removeColorBackground( cv::Mat&, cv::Mat& );
void removeUniformColorBackground_D2(Mat&, Mat&);
void removeUniformColorBackground_D3(Mat&, Mat&);
void removeUniformColorBackground_Q3(Mat&, Mat&);
void removeUniformColorBackground(Mat&, Mat&);
void ImLabelFlatZones_WithCriterion(const Mat &, const Mat &, const ZoneInfoType, Mat &);
// Constants
// TODO: verify if these values make sense
// Parameters for background subtraction
const int BS_HISTORY = 200;
const float BS_THRESHOLD = 16.0;
const bool BS_SHADOW = false;
const double LEARNING_RATE = 0.0;
const int MORPHO_SIZE = 4;

// For detection of head profiles
const double INCREASE_WIDTH = 0.9;
const double MINBOX_PERCENTAGE = 0.05;

#endif
