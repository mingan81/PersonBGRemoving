#include "opencv2/core/core.hpp"
#include "opencv2/contrib/contrib.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/objdetect/objdetect.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <array>

//#include "../include/core.hpp"
#include "../include/humanmodel.hpp"

using namespace cv;
using namespace std;

#define COL 1
#define ROW 2

// Please specify your own file path here.
// Please do not hardcode any self-defined variables in the code. 
string HAARCASCADE_PATH = "/Users/kevinwang264/QtProject/ProjectDemo/";
string IMG_PATH = HAARCASCADE_PATH + "images/";

std::pair<int, int> MeasMinMax(const Mat & imIn)
{	

	int mini = std::numeric_limits<int>::max();
	int maxi = 0;

	const int w = imIn.cols, h = imIn.rows;
	const int * p_in = &imIn.at<int>(0, 0);
	const int * p_inend = p_in + w*h;

	int counter = 0;

	for (; p_in != p_inend; p_in++)
	{
		if (*p_in < mini) mini = *p_in;
		if (*p_in > maxi) maxi = *p_in;
		counter++;
	}
	std::pair<int, int> minmax;
	minmax.first = mini;
	minmax.second = maxi;
	return minmax;
}




cv::Rect detectFace(cv::Mat& frame,cv::Mat& crop)
{
	//string face_cascade_name = "c:/opencv/sources/data/haarcascades/haarcascade_frontalface_default.xml";
	//std::string face_cascade_name = "/usr/share/opencv/haarcascades/haarcascade_frontalface_default.xml";
	std::string face_cascade_name = HAARCASCADE_PATH + "haarcascade_frontalface_default.xml";
	CascadeClassifier face_cascade;
	if (!face_cascade.load(face_cascade_name))
	{
		printf("--(!)Error loading\n");
	};

	
	std::vector<Rect> faces;
	cv::Mat frame_gray; //grayscale image 

	cv::Mat res;
	cv::Mat gray;
	string text;
	stringstream sstm;
	
	// RGB TO GRAYSCALE IMAGE
	cv::cvtColor(frame, frame_gray, COLOR_BGR2GRAY);
	// HISTGRAM EQUALIZATION 
	equalizeHist(frame_gray, frame_gray);
	imwrite(IMG_PATH + "gray.png",frame_gray);

	// DETECT FACE USING CASCADE CLASSIFIER
	face_cascade.detectMultiScale(frame_gray, faces, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(frame.cols/10, frame.rows/10));
	//face_cascade.detectMultiScale(frame_gray, faces, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));

	//SET ROI
	cv::Rect roi_b; //Rect for biggest ROI 
	cv::Rect roi_c; // temp Rect for current ROI
	cv::Rect roi_head; // Rect which contains head

	size_t ic = 0; // ic is index of current ROI
	int ac = 0; // ac is area of current ROI

	size_t ib = 0; // ib is index of biggest ROI
	int ab = 0; // ab is area of biggest ROI
	std::cout << faces.size() << endl;
	//~ cin.ignore(1); //pause for debugging

	for (ic = 0; ic < faces.size(); ic++) // Iterate through all current ROIs (detected faces)
	{
		roi_c.x = faces[ic].x;
		roi_c.y = faces[ic].y;
		roi_c.width = (faces[ic].width);
		roi_c.height = (faces[ic].height);

		ac = roi_c.width * roi_c.height; // Get the area of current ROI (detected face)

		roi_b.x = faces[ib].x;
		roi_b.y = faces[ib].y;
		roi_b.width = (faces[ib].width);
		roi_b.height = (faces[ib].height);

		ab = roi_b.width * roi_b.height; // Get the area of biggest ROI, at beginning it is same as "current" ROI

		if (ac > ab)
		{
			ib = ic;
			roi_b.x = faces[ib].x;
			roi_b.y = faces[ib].y;
			roi_b.width = (faces[ib].width);
			roi_b.height = (faces[ib].height);
		}
	}
	
	roi_head.x = roi_b.x; 
	roi_head.y = 0; // modify it when the suitable method to set this point
	roi_head.width = roi_b.width;
	roi_head.height = roi_b.y + roi_b.height;

	crop = frame(roi_head);
	return roi_head;
	
}

cv::Rect detectFace_dgil(cv::Mat& frame,cv::Mat& crop)
{
	// Diana modified this ------------------
	std::string face_cascade_name = HAARCASCADE_PATH + "haarcascade_frontalface_default.xml";
	CascadeClassifier face_cascade;
	if (!face_cascade.load(face_cascade_name))
	{
		cerr << "detectFace_dgil(): Error loading cascade file\n";
		exit(EXIT_FAILURE);
	};
	// ---------------------------------------

	std::vector<Rect> faces;
	cv::Mat frame_gray; //grayscale image

	cv::Mat res;
	cv::Mat gray;
	string text;
	stringstream sstm;

	// RGB TO GRAYSCALE IMAGE
	cv::cvtColor(frame, frame_gray, COLOR_BGR2GRAY);
	// HISTGRAM EQUALIZATION
	equalizeHist(frame_gray, frame_gray);
	imwrite(IMG_PATH + "gray.png",frame_gray);

	// DETECT FACE USING CASCADE CLASSIFIER
	face_cascade.detectMultiScale(frame_gray, faces, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(frame.cols/10, frame.rows/10));
	//face_cascade.detectMultiScale(frame_gray, faces, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));

	//SET ROI
	cv::Rect roi_b; //Rect for biggest ROI
	cv::Rect roi_c; // temp Rect for current ROI
	cv::Rect roi_head; // Rect which contains head

	size_t ic = 0; // ic is index of current ROI
	int ac = 0; // ac is area of current ROI

	size_t ib = 0; // ib is index of biggest ROI
	int ab = 0; // ab is area of biggest ROI
	std::cout << faces.size() << endl;
	//~ cin.ignore(1); //pause for debugging

	for (ic = 0; ic < faces.size(); ic++) // Iterate through all current ROIs (detected faces)
	{
		roi_c.x = faces[ic].x;
		roi_c.y = faces[ic].y;
		roi_c.width = (faces[ic].width);
		roi_c.height = (faces[ic].height);

		ac = roi_c.width * roi_c.height; // Get the area of current ROI (detected face)

		roi_b.x = faces[ib].x;
		roi_b.y = faces[ib].y;
		roi_b.width = (faces[ib].width);
		roi_b.height = (faces[ib].height);

		ab = roi_b.width * roi_b.height; // Get the area of biggest ROI, at beginning it is same as "current" ROI

		if (ac > ab)
		{
			ib = ic;
			roi_b.x = faces[ib].x;
			roi_b.y = faces[ib].y;
			roi_b.width = (faces[ib].width);
			roi_b.height = (faces[ib].height);
		}
	}

	roi_head.x = roi_b.x;
	roi_head.y = 0; // modify it when the suitable method to set this point
	roi_head.width = roi_b.width;
	roi_head.height = roi_b.y + roi_b.height;

	crop = frame(roi_head);
	return roi_head;

}

cv::Rect detectNose(cv::Mat& frame,cv::Mat& crop)
{
	//dbi.WriteEnter("Enter the Detectnose\n");
	//dbi.WriteImage(frame,"frame.png");
	string nose_cascade_name = HAARCASCADE_PATH + "haarcascade_mcs_nose.xml";
	//std::string nose_cascade_name = "/usr/share/opencv/haarcascades/haarcascade_mcs_nose.xml";
	CascadeClassifier nose_cascade;
	if (!nose_cascade.load(nose_cascade_name))
	{
		printf("--(!)Error loading\n");
	};

	
	std::vector<Rect> noses;
	cv::Mat frame_gray; //grayscale image
	//cv::Mat crop; // the image which has human head  

	cv::Mat res;
	cv::Mat gray;
	string text;
	stringstream sstm;
	
	// RGB TO GRAYSCALE IMAGE
	cv::cvtColor(frame, frame_gray, COLOR_BGR2GRAY);
	// HISTGRAM EQUALIZATION 
	equalizeHist(frame_gray, frame_gray);
	//imshow(IMG_PATH + "gray",frame_gray);

	// DETECT NOSE USING CASCADE CLASSIFIER
	nose_cascade.detectMultiScale(frame_gray, noses, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(30, 30));	

	//SET ROI
	cv::Rect roi_b; //Rect for biggest ROI 
	cv::Rect roi_c; // temp Rect for current ROI
	cv::Rect roi_nose; // Rect which contains head

	size_t ic = 0; // ic is index of current ROI
	int ac = 0; // ac is area of current ROI

	size_t ib = 0; // ib is index of biggest ROI
	int ab = 0; // ab is area of biggest ROI
	std::cout << noses.size() << endl;
	//~ cin.ignore(1); //pause for debugging

	for (ic = 0; ic < noses.size(); ic++) // Iterate through all current ROIs (detected faces)
	{
		roi_c.x = noses[ic].x;
		roi_c.y = noses[ic].y;
		roi_c.width = (noses[ic].width);
		roi_c.height = (noses[ic].height);

		ac = roi_c.width * roi_c.height; // Get the area of current ROI (detected face)

		roi_b.x = noses[ib].x;
		roi_b.y = noses[ib].y;
		roi_b.width = (noses[ib].width);
		roi_b.height = (noses[ib].height);
		

		ab = roi_b.width * roi_b.height; // Get the area of biggest ROI, at beginning it is same as "current" ROI

		Point pt1(noses[ic].x, noses[ic].y); // Display detected faces on main window
		Point pt2((noses[ic].x + noses[ic].height), (noses[ic].y + noses[ic].width));
		cv::rectangle(frame, pt1, pt2, Scalar(0, 255, 0), 2, 8, 0);

		if (ac < ab)
		{
			ib = ic;
			roi_b.x = noses[ib].x;
			roi_b.y = noses[ib].y;
			roi_b.width = (noses[ib].width);
			roi_b.height = (noses[ib].height);
		}
	}
	
	roi_nose.x = roi_b.x; 
	roi_nose.y = roi_b.y; // modify it when the suitable method to set this point
	roi_nose.width = roi_b.width;
	roi_nose.height = roi_b.height;

	crop = frame(roi_nose);
	
	
	imshow(IMG_PATH + "original", frame);
	return roi_nose;
}

void GetVerticalEdge(cv::Mat& binary,cv::Mat& binary_return){ 
	
	imwrite(IMG_PATH + "binary.png",binary);	
	binary_return = Mat::zeros(binary.rows,binary.cols,CV_8U);
	//imwrite(IMG_PATH + "binary_return.png",binary_return);	
	int i,j;
	cout<<binary.cols <<"    "<< binary.rows<<endl;
	cout<<binary.at<uchar>(510,310)<< endl;
	for(i=1;i<binary.rows-1;i++){
		for(j=1;j<binary.cols-1;j++){
			if(binary.at<uchar>(i,j) !=0 && binary.at<uchar>(i-1,j) ==255 && binary.at<uchar>(i+1,j) ==255 ){ 
				binary_return.at<uchar>(i,j) = 255;
			}
			
		}
		
	}
	imwrite(IMG_PATH + "binary_return.png",binary_return);	
}

void FindBlobs(const cv::Mat &binary, std::vector < std::vector<cv::Point2i> > &blobs)
{
    blobs.clear();

    // Fill the label_image with the blobs
    // 0  - background
    // 1  - unlabelled foreground
    // 2+ - labelled foreground

    cv::Mat label_image;
    binary.convertTo(label_image, CV_32SC1);

    int label_count = 2; // starts at 2 because 0,1 are used already

    for(int y=0; y < label_image.rows; y++) {
        int *row = (int*)label_image.ptr(y);
        for(int x=0; x < label_image.cols; x++) {
            if(row[x] != 1) {
                continue;
            }

            cv::Rect rect;
            cv::floodFill(label_image, cv::Point(x,y), label_count, &rect, 0, 0, 4);

            std::vector <cv::Point2i> blob;

            for(int i=rect.y; i < (rect.y+rect.height); i++) {
                int *row2 = (int*)label_image.ptr(i);
                for(int j=rect.x; j < (rect.x+rect.width); j++) {
                    if(row2[j] != label_count) {
                        continue;
                    }

                    blob.push_back(cv::Point2i(j,i));
                }
            }

            blobs.push_back(blob);

            label_count++;
        }
    }
}

std::vector<int> CountBrightPixels(cv::Mat image, int ColOrRow = 1, int numbers = 0) { //specify the col or row want to count

	Mat binary_image = image.clone();

    //Convert color image to binary image if the original image is not grayscale
	if (image.channels() != 1)
    {
    	//Grayscale matrix
	    Mat gray_image(image.size(), CV_8U);

	    //Convert RGB to Gray
	    cvtColor(image, gray_image, CV_BGR2GRAY);
	    
	    //Binary image
	    //Mat binary_image(gray_image.size(), gray_image.type());

	    //Apply thresholding
	    threshold(gray_image, binary_image, 100.0, 255.0, THRESH_BINARY);	//TODO
    }

    //cout << "M=" << endl << " " << binary_image << endl << endl;

    imshow(IMG_PATH +  "Binary Image", binary_image );			// Show our image inside it.
 
    vector<int> counts;
    counts.clear();

    switch(ColOrRow)
    {
        case COL:
        {
            int colCounts;  //Accumulated if > 0
            
            if (numbers == 0) numbers = binary_image.cols;
           	
           	//Set numbers to MAX column numbers if the parameter 'numbers' specified is out of range
           	try{
           		if (numbers > binary_image.cols) {
           			throw "Column numbers out of range, set 'numbers' to image max column numbers";
           		}
           	} catch (const char* msg){
           		cerr << msg << endl;
           		numbers = binary_image.cols;
           	}

            for (int i = 0; i < numbers; ++i){
                colCounts = 0;
                for (int j = 0; j < binary_image.rows; ++j){       
                    if (binary_image.at<uchar>(i,j) != 0){
                        colCounts += 1;
                    }
                }
                counts.push_back(colCounts);
            }
        }
        break;

        case ROW:
        {
            int rowCounts;  //Accumulated if > 0
           
            if (numbers == 0) numbers = binary_image.rows;

            //Set numbers to MAX Row numbers if the parameter 'numbers' specified is out of range
           	try{
           		if (numbers > binary_image.rows) {
           			throw "ow numbers out of range, set 'numbers' to image max Row numbers";
           		}
           	} catch (const char* msg){
           		cerr << msg << endl;
           		numbers = binary_image.rows;
           	}
            
            for (int i = 0; i < numbers; ++i){
                rowCounts = 0;
                for (int j = 0; j < binary_image.cols; ++j){       
                    if (binary_image.at<uchar>(i,j) != 0){
                        rowCounts += 1;
                    }
                }
                counts.push_back(rowCounts);
            }
        }
        break;

        default:
            break;
    }

    return counts;
}

std::vector<int> CountBrightPixelsCol(cv::Mat& binary){ //specify the col or row want to count
	
	imshow(IMG_PATH + "binary",binary);
	std::vector<int> count;count.clear();
	int countCols;
	int i,j;
	for(i=0;i<binary.cols;i++){
		countCols = 0;
		for(j=0;j<binary.rows;j++){
			if(binary.at<uchar>(j,i) !=0){ 
				countCols = countCols +1;
			}
			
		}
		cout<< countCols <<endl;
		count.push_back(countCols);
	}
	
	cout<<"size" << count.size()<<endl;
	return count;
	
	
	
}


cv::Mat HeadForegroundExtraction(cv::Mat& frame){
	
	//imshow(IMG_PATH + "cont", frame);

	Mat frame2 = frame.clone(); 
	Mat frame_gray;
	cv::cvtColor(frame, frame_gray, COLOR_BGR2GRAY);
	equalizeHist(frame_gray, frame_gray);
	
	Mat edge_detected;
	cv::Canny(frame_gray, edge_detected, 40, 120);
	
	//imwrite(IMG_PATH + "edge.png",edge_detected);
	//Mat vertical_edge;
	//GetVerticalEdge(edge_detected,vertical_edge);
	//imwrite(IMG_PATH + "vertical_edge.png",vertical_edge);


	
	//std::vector<int> count;
	//count = CountBrightPixelsCol(vertical_edge);
	
	
	Mat result;
	//Mat result(frame.rows,frame.cols,CV_8U);
	//for (int i = frame.rows)
	Mat backgroundModel, foregroundModel;
	Rect rectangle(1, 1, frame.cols-1, frame.rows-1);	
	
	//~ Mat YCrCb;
	//~ cv::cvtColor(frame,YCrCb,CV_BGR2YCrCb);
	//~ Mat color_Mask;
	//~ cv::inRange(YCrCb,Scalar(80,133,77),Scalar(255,173,127),color_Mask);
	//~ imshow(IMG_PATH + "skincolor", color_Mask);
	
	//~ ellipse(skinCrCbHist, Point(113, 155.6), Size(23.4, 15.2), 43.0, 0.0, 360.0, Scalar(255, 255, 255), -1);
//~ 
	//~ bool skincolor;
	//~ Mat YCrCb;
	//~ cv::cvtColor(frame,YCrCb,CV_RGB2YCrCb);
	//~ Mat color_Mask;
	//~ for(int row = 0; row < YCrCb.rows; ++row) {
		//~ for(int col = 0; col < YCrCb.cols; ++col) {
			//~ Vec3f colorY = YCrCb.at<Vec3f>(col,row);
			//~ skincolor = isSkin(Scalar(colorY.val[0],colorY.val[1],colorY.val[2]));
			//~ detected_gray.at<uchar>(col,row) = 254*(int)skincolor;
		//~ }
	//~ }
	//~ imshow(IMG_PATH + "skincolor", detected_gray);
	 
	//~ 
	//~ 
	//~ 
//~ 
//~ 
	cv::grabCut(frame, result, rectangle, backgroundModel, foregroundModel, 7, cv::GC_INIT_WITH_RECT);
	cv::compare(result, cv::GC_PR_FGD, result, cv::CMP_EQ);	
	//imshow(IMG_PATH + "result", result);
	return result;
//~ 
	//~ cv::Mat foreground(frame.size(), CV_8UC3, cv::Scalar(255, 255, 255));
	//~ //cv::Mat background(image.size(),CV_8UC3,cv::Scalar(255,255,255));
	//~ frame.copyTo(foreground, result); 
//~ 
//~ 
	//~ cv::Mat cont, threshold;
	//~ cont = foreground.clone();
	//~ cv::cvtColor(cont, cont, cv::COLOR_BGR2GRAY);
//~ 
	//~ cv::vector<cv::Vec4i> hierarchy;
	//~ cv::vector<cv::vector<cv::Point> > contours;
//~ 
	//~ cv::threshold(cont, cont, 128, 255, cv::THRESH_BINARY);
	//~ threshold = cont.clone();
//~ 
	//~ cv::findContours(cont, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);
//~ 
	//~ for (size_t i = 0; i < contours.size(); i++){
		//~ cv::drawContours(foreground, contours, i, cv::Scalar(0, 0, 255), 2);
	//~ }
//~ 
	//~ cv::imshow(IMG_PATH + "Contour", threshold);
	//~ //cv::imshow(IMG_PATH + "Video", frame);
	//~ // draw rectangle on original image
	//~ //cv::rectangle(image, rectangle, cv::Scalar(255, 255, 255), 1);
	//~ 
//~ 
	//~ imshow(IMG_PATH + "Foreground", foreground);	
}

//~ 
//~ bool isSkin(const Scalar color) {
    //~ Mat input = Mat(Size(1, 1), CV_8UC3, color);
    //~ Mat output;
//~ 
    //~ cvtColor(input, output, CV_BGR2YCrCb);
//~ 
    //~ Vec3b ycrcb = output.at<Vec3b>(0, 0);
    //~ return ((skinCrCbHist.at<uchar>(ycrcb[1], ycrcb[2]) > 0));
//~ }




int ImLabel_QG(Mat& imIn , int startlabel, Mat& imLabel)
{
	//dbi.WriteEnter("Entering ImLabel\n");
	imLabel.empty();
	int i=0, j=0, ii=0, jj=0, ilocal = 0, jlocal = 0;
	int w = imIn.cols, h = imIn.rows;
	int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
	int rayon = 1;
	int current_label = startlabel;
	// 	dbi.WriteInfo("debug 0.5\n");
	//std::queue<Node<uchar, int > > q;
	std::vector<Node<uchar, int > > q; q.clear();
	// 	dbi.WriteInfo("debug 1\n");
	// 	dbi.WriteInfo("here 1\n");
	for (j=0; j<h; j++){
		for (i=0; i<w; i++){
			if (imIn.at<uchar>(i,j) != 0){
				if (imLabel.at<int>(i,j) == 0){
					// 					dbi.WriteInfo("debug 1.5\n");
					Node<uchar, int> noeud(i,j,0,current_label);
					q.push_back(noeud);
					imLabel.at<int>(i,j) = current_label;
					// 					dbi.WriteInfo("debug 2\n");
					// 			dbi.WriteInfo("here 2\n");
					while (!q.empty())
					{
						Node<uchar, int> current_node = q.back();
						q.pop_back();
						ilocal = current_node.GetXPos();
						jlocal = current_node.GetYPos();

						borneInfx = std::max(0,ilocal-rayon);
						borneInfy = std::max(0,jlocal-rayon);
						borneMaxx = std::min(w, ilocal+rayon+1);
						borneMaxy = std::min(h, jlocal+rayon+1);

//						if (borneMaxx > w) dbi.WriteInfo("pb largeur\n");
//      				if (borneMaxy > h) dbi.WriteInfo("pb hauteur\n");

						for (ii=borneInfx; ii<borneMaxx; ii++){
							for (jj=borneInfy; jj<borneMaxy; jj++){
								if((imIn.at<uchar>(ii,jj) != 0) && (imLabel.at<int>(ii,jj) == 0)){
									Node<uchar, int> noeudlocal(ii,jj,0,current_label);
									q.push_back(noeudlocal);
									imLabel.at<int>(ii,jj) = current_label;
								}
							}
						}
					}
					// 			dbi.WriteInfo("here 3\n");
					current_label++;
					// 					dbi.WriteInfo("debug 3\n");
				}
			}
		}
	}
	//dbi.WriteOut("Leaving ImLabel\n");
	return (int(current_label-1));
}

// -----------------------------------------------------------------
//			 ESTIMATE FRAME RANGE TO DETECT HEAD PROFILES
// This function estimates the starting and the ending frames
// corresponding to head movement from left to right (or vice versa)
// The ultimate goal is to simplify the head profiles detection
//
// The function pass a array contains two values, which are the
// results of audio analysis. These values are the number of seconds
// in a video file when the person says "start" and "end" or similar
// -----------------------------------------------------------------
std::array<int, 2> getStartEnd( cv::VideoCapture& videoObj, const array<int, 2> timeRes )
{

	if (timeRes[0] < 0 || timeRes[1] < 0)
	{
		cerr << "Not valid video" << endl;
		exit(EXIT_FAILURE);
	}

	// Error in opening the video input
	if (!videoObj.isOpened())
	{
		cerr << "getStartEnd(): Unable to open video" << endl;
		exit(EXIT_FAILURE);
	}

	// Restart the video, just in case
	videoObj.set(CV_CAP_PROP_POS_AVI_RATIO, 0);

	//-- Convert seconds to frame numbers
	double frameRate = videoObj.get(CV_CAP_PROP_FPS);
	double nFrames = videoObj.get(CV_CAP_PROP_FRAME_COUNT);
	int startFrame = floor(timeRes[0] * frameRate);;
	int endFrame = floor(timeRes[1] * frameRate);

	// ***** TESTING ****
	// Video headers are sometimes corrupted or not set correctly
	if (frameRate >= 1000)
	{
		double temp = frameRate;
		frameRate = 30.0;
        startFrame = floor(timeRes[0] * frameRate);
		endFrame = floor(timeRes[1] * frameRate);
		nFrames = floor((nFrames/temp) * frameRate);
	}

	cout << "Frame rate = " << frameRate << endl
		 << "Starting frame = " << startFrame << endl
		 << "Ending frame = " << endFrame << endl
		 << "Number of frames = " << nFrames << endl;

	if (endFrame > nFrames)
	{
		cerr << "getStartEnd(): Final frame not valid" << endl;
		exit(EXIT_FAILURE);
	}

	array<int, 2> startEnd = {{startFrame, endFrame}};

	return startEnd;
}

// -----------------------------------------------------------------
//				CREATE VECTOR OF HEAD IMAGES FROM A VIDEO
// This function creates a vector of head images from a video input,
// following these steps:
// 1. Detect a face and get head bounding box
// 2. Enlarge bounding box
// 3. For each next frame, save original head image in a vector
//
// The function only processes the frames corresponding to head movement
// from left to right (already defined by getStartEnd() function)
//
// Parameters:
//	+ (input)	videoObj is the video object
//	+ (input)	startEnd contains the starting and ending frame numbers
//			 	corresponding to head movement
//	+ (output)	headImages contains the original head images
// -----------------------------------------------------------------
void video2HeadImages( cv::VideoCapture& videoObj, const std::array<int, 2>& startEnd, std::vector<cv::Mat>& headImages )
{
	// Error in opening the video input
	if (!videoObj.isOpened())
	{
		cerr << "video2HeadImages(): Unable to open video" << endl;
		exit(EXIT_FAILURE);
	}

	// Restart the video, just in case
	videoObj.set(CV_CAP_PROP_POS_AVI_RATIO, 0);

	Mat frame;
	Mat crop;
	Rect headBox;
	Rect largeHeadBox;
	int index = 0;

	// ------------------ VIDEO PROCESSING ------------------------
	for (int frameNumber = 1; frameNumber <= startEnd[1]; frameNumber++)
	{
		videoObj >> frame;

		if (!frame.empty() && !headBox.width)
		{
			// -- Detect face until one is found
			crop.empty();
			headBox = detectFace_dgil( frame, crop );

			if (!headBox.width)
			{
				cout << "Face not detected: " << frameNumber << endl;
				continue;
			}

			// Get number of columns of the frames
			double nCols = videoObj.get(CV_CAP_PROP_FRAME_WIDTH);

			// -- Enlarge bounding box to include head profile
			largeHeadBox = enlargeBBox( headBox, INCREASE_WIDTH, (int)nCols );
		}

		// -- Save original head image in a vector
		if (!frame.empty() && frameNumber >= startEnd[0])
		{
			// There is a special behaviour with Mat variables
			// It is necessary to have a separate image and copy it to the vector
			Mat tempHead = frame(largeHeadBox);
			tempHead.copyTo(headImages.at(index));
			index++;
		}
	}
}

// -----------------------------------------------------------------
//					REMOVE BACKGROUND
// Adapted from OpenCV Tutorials --> Video Analysis --> How to Use
// Background Subtraction Methods
// -----------------------------------------------------------------
cv::Ptr<BackgroundSubtractor> createBgSubtractor()
{
	// Create Background Subtractor object with MOG2 approach
	// Parameters are defined as constants in humanmodel.hpp
	Ptr<BackgroundSubtractor> bsMOG2;
	bsMOG2 = new BackgroundSubtractorMOG2(BS_HISTORY, BS_THRESHOLD, BS_SHADOW);

	return bsMOG2;
}

cv::Mat removeBackground( const cv::Mat& frame, cv::Ptr<BackgroundSubtractor> bgSub )
{
	// Foreground mask
	Mat fgMask;

	// Remove background and get foreground mask
	bgSub->operator()(frame, fgMask, LEARNING_RATE);

	// Apply opening to reduce noise
	int morphoType = MORPH_ELLIPSE;
	Mat element = getStructuringElement( morphoType, Size( 2*MORPHO_SIZE + 1, 2*MORPHO_SIZE + 1 ), Point( MORPHO_SIZE, MORPHO_SIZE ) );
	morphologyEx( fgMask, fgMask, 2, element ); // 2 indicates opening operation

	return fgMask;
}

// -----------------------------------------------------------------
//					ENLARGE BOUNDING BOX
// Enlarge width of a bounding box according to a percentage:
//
//       <-- original width -->
// ++++++----------------------++++++
// |    |                      |    |
// |    |                      |    |
// |    |                      |    |
// ++++++----------------------++++++
// added                       added
// <----------- new width ---------->
// -----------------------------------------------------------------
cv::Rect enlargeBBox( const cv::Rect& bBox, double increaseWidth, int nCols )
{
	Rect largeBBox;

	largeBBox.x = floor(bBox.x - bBox.width*increaseWidth/2);
	largeBBox.y = bBox.y;
	largeBBox.width = floor(bBox.width*(1 + increaseWidth));	// Width is increased
	largeBBox.height = bBox.height;

	if (largeBBox.x < 0)  
	{	
		largeBBox.x = 0;
	}

	if ((largeBBox.x + largeBBox.width) > nCols) 
	{
		largeBBox.width = nCols - largeBBox.x;
	}

	return largeBBox;
}

// -----------------------------------------------------------------
//					GET FRAME NUMBERS OF HEAD PROFILES
// To get the two maximum peaks that indicate the head profiles, an
// estimation of the frontal head is needed. This estimation is based
// on the minimum area of local minimums.
// One profile is determined as the maximum before the frame number
// indicating the frontal head. The other profile is gotten as the
// maximum after the frontal head.
// -----------------------------------------------------------------
vector<int> getPosProfiles( const vector<int>& maxPosFrame, const vector<int>& maxArea, int posMinArea )
{
	// Get the index of maximum value before minimum area
	Mat minAreaMask = (Mat)maxPosFrame < posMinArea; 			// Result: 255 & 0's
	int indexBefore = countNonZero(minAreaMask);				// Frame numbers are already organized
	auto indexLeftMax = distance( maxArea.begin(), max_element( maxArea.begin(), maxArea.begin() + indexBefore ) );

	// Get the index of maximum value after minimum area
	auto indexRightMax = distance( maxArea.begin(), max_element( maxArea.begin() + indexBefore, maxArea.end() ) );

	vector<int> profileFrames;
	profileFrames.push_back(maxPosFrame.at(indexLeftMax));
	profileFrames.push_back(0);
	profileFrames.push_back(maxPosFrame.at(indexRightMax));

	return profileFrames;
}

// -----------------------------------------------------------------
//  				GET THE HEAD IMAGE WITHOUT BACKGROUND
// This function just gets the head image without background, according
// to a mask already computed with the foreground extraction algorithm.
// So, two images are loaded: a head mask and a head image with background.
// - Options:
// If imgName is empty, the function reads the first two Mat variables.
// Otherwise, the function loads images from the computer.
// -----------------------------------------------------------------
cv::Mat getHeadImg( const cv::Mat& imgVar, const cv::Mat& maskVar, const int& frameNumber, string imgName )
{
	Mat img, mask;

	if (imgName.empty())
	{
		// Read image variables
		img = imgVar;
		mask = maskVar;
	}
	else
	{
		// If images are stored in the computer, read them
		img = imread(IMG_PATH + imgName + std::to_string(frameNumber) + ".png");
		mask = imread(IMG_PATH + "mask" + std::to_string(frameNumber) + ".png", CV_LOAD_IMAGE_GRAYSCALE);
	}

	if (img.empty() && mask.empty())
	{
		cerr << "getHeadImg(): Images # " << frameNumber << " not loaded" << endl;
		exit(EXIT_FAILURE);
	}

	if (mask.channels() != 1)
	{
		cerr << "getHeadImg(): Mask image # " << frameNumber << "not valid" << endl;
		exit(EXIT_FAILURE);
	}

	// If mask pixel values are white, set new pixel values to the original ones
	Mat mask3, headImg;
	cvtColor(mask, mask3, CV_GRAY2RGB);
	bitwise_and(img, mask3, headImg);

	return headImg;
}

// -----------------------------------------------------------------
//  						DETECT EAR
// Detect ear in order to determine if it's the right or the left profile
// Left or right ears depend on the cascade file loaded
// -----------------------------------------------------------------
bool detectEar( cv::Mat& headImg, std::string cascade_name )
{
	bool isEarDetected = false;
	CascadeClassifier ear_cascade;

	// Load cascade file provided by OpenCV
    if (!ear_cascade.load( cascade_name ))
    {
    	cerr << "detectEar(): **** Error loading ear cascade ****" << endl;
    	exit(EXIT_FAILURE);
    };

    std::vector<Rect> ears;
    Mat head_gray;

    if (headImg.channels() == 3)
    {
    	cvtColor( headImg, head_gray, COLOR_BGR2GRAY );
    }
    else
    {
    	head_gray = headImg;
    }

    equalizeHist( head_gray, head_gray );

    //-- Detect ear
    ear_cascade.detectMultiScale( head_gray, ears, 1.1, 2, 0|CASCADE_SCALE_IMAGE, Size(30, 30) );

    for (size_t i = 0; i < ears.size(); i++)
    {
    	Point center( ears[i].x + ears[i].width/2, ears[i].y + ears[i].height/2 );
    	ellipse( headImg, center, Size( ears[i].width/2, ears[i].height/2), 0, 0, 360, Scalar( 255, 0, 255 ), 4, 8, 0 );
    }

    if (ears.size())
    {
    	isEarDetected = true;
    	//imshow(IMG_PATH +  "Ear detected", headImg );
    	//waitKey();
    }

	return isEarDetected;
}

// -----------------------------------------------------------------
// 						IDENTIFY HEAD PROFILES
// This function identifies if it is the right or the left head profile
// using ear detection
// Left or right: view from an external subject
// -----------------------------------------------------------------
int identifyProfiles( cv::Mat& headImg1, cv::Mat& headImg2, std::vector<int>& headFrameNumbers )
{
	string leftCascadeName = HAARCASCADE_PATH + "haarcascade_mcs_leftear.xml";
	string rightCascadeName = HAARCASCADE_PATH + "haarcascade_mcs_rightear.xml";
	int idProfile = 0; // idea

	// Sometimes ear is not well detected
	// So, if just one is detected, we assume the other image contains the opposite ear
	// The final order of head frame numbers is: left profile, frontal view, right profile, frontal view
	if ( detectEar( headImg1, leftCascadeName ) ||  detectEar( headImg2, rightCascadeName ) )
	{
		idProfile = 1;
	}
	else if ( detectEar( headImg2, leftCascadeName ) || detectEar( headImg1, rightCascadeName ) )
	{
		vector<int> tempNumbers = headFrameNumbers;
		headFrameNumbers.at(LEFT) = headFrameNumbers.at(RIGHT); // left and right changed
		headFrameNumbers.at(RIGHT) = tempNumbers.at(LEFT);
		headFrameNumbers.at(FRONT1) = headFrameNumbers.at(FRONT2);
		headFrameNumbers.at(FRONT2) = tempNumbers.at(FRONT1);

		idProfile = 2;
	}
	else
	{
		cout <<  "identifyProfiles(): Ears not detected" << endl;	// TODO: what to do if ears are not detected?
		return idProfile;
	}

	cout << "Image # " << headFrameNumbers.at(LEFT) << ": left profile" << endl;
	cout << "Image # " << headFrameNumbers.at(FRONT1) << ": central view" << endl;
	cout << "Image # " << headFrameNumbers.at(RIGHT) << ": right profile" << endl;
	cout << "Image # " << headFrameNumbers.at(FRONT2) << ": central view" << endl;

	return idProfile;
}

// -----------------------------------------------------------------
//						DETECT HEAD PROFILES
// This function gets 4 key frame numbers from a video, corresponding
// to different positions of the head:
// left profile, frontal view, right profile and a second frontal view,
// respectively
//
// Basic steps are detailed below:
// 1. Detect face and get head bounding box
// 2. Enlarge bounding box
// 3. Get all local maxima and local minima
//	  3.1 Extract motion head mask (with background subtraction)
// 	  3.2 Compute current area of head
//    3.3 Extract head image and its mask
// 4. Get frame positions of head profiles and frontal views
// 5. Get corresponding head masks
//
// Parameters:
//	+ (input)	videoObj is the video object
//	+ (input)	startEnd contains the starting and ending frame numbers
//			 	corresponding to head movement
//	+ (output)	extHeads contains the head masks and key frame numbers
// ---------------------------------------------------------------------
void detectHeadProfiles( cv::VideoCapture& videoObj, const std::array<int, 2>& startEnd, HeadSequence& extHeads )
{

	// Error in opening the video input
	if (!videoObj.isOpened())
	{
		cerr << "detectHeadProfiles(): Unable to open video" << endl;
		exit(EXIT_FAILURE);
	}

	// Restart the video, just in case
	videoObj.set(CV_CAP_PROP_POS_AVI_RATIO, 0);

	Mat frame;
	Mat crop;
	Rect headBox;

	// For background subtraction
	Ptr<BackgroundSubtractor> bgSub = createBgSubtractor();

	// For head profile detection
	Rect largeHeadBox;
	int minAreaThreshold;
	Mat previousHead, motionMask;
	int previousArea, previousSlope;
	vector<int> maxPosFrame, maxArea;
	vector<int> minPosFrame, minArea;
	int posMinMin = 0;
	int frameNumberFace = 0;

	// For tests
	//ofstream myfile;
	//myfile.open ("Areas.csv");

	// ------------- VIDEO PROCESSING ------------------------
	for (int frameNumber = 1; frameNumber <= startEnd[1]; frameNumber++)
	{
		videoObj >> frame;

		// -- Detect face and initialize variables for profile detection
		if (!frame.empty() && !headBox.width)
		{
			crop.empty();
			headBox = detectFace_dgil(frame,crop);
			cout << "headBox = " << headBox << endl;
			if (!headBox.width)
			{
				cout << "Face not detected: " << frameNumber << endl;
				continue;
			}
			//imshow(IMG_PATH + "head",crop);
			frameNumberFace = frameNumber;

			// Get number of columns of the frames
			double nCols = videoObj.get(CV_CAP_PROP_FRAME_WIDTH);

			// Enlarge bounding box to include head profile
			largeHeadBox = enlargeBBox( headBox, INCREASE_WIDTH, (int)nCols );
			cout << largeHeadBox << endl;

			minAreaThreshold = MINBOX_PERCENTAGE * largeHeadBox.width * largeHeadBox.height;
			cout << "Area threshold = " << minAreaThreshold << endl;

			// Initialize previous values
            previousHead = frame(largeHeadBox);
            Mat previousMask = Mat::zeros(previousHead.rows,previousHead.cols, CV_16UC1);
			//previousMask = HeadForegroundExtraction( previousHead );
            //removeColorBackground( previousHead, previousMask );
            //removeUniformColorBackground(previousHead, previousMask);
            //removeUniformColorBackground_D2(previousHead, previousMask);
            removeUniformColorBackground_D3(previousHead, previousMask);
            //removeUniformColorBackground_Q3(previousHead, previousMask);

            previousArea = 0;
			previousSlope = 1;
		}

		// -- Get foreground/motion mask (head image)
		motionMask = removeBackground( frame(largeHeadBox), bgSub );

		// -- Get local maximums and minimums for the first seconds of the video
		// A local maximum may indicate a left or right profile
		// A local minimum may indicate a frontal head, and it helps to separate
		// left and right profiles
		if (!frame.empty() && frameNumber >= startEnd[0])
		{
			// Calculate area for the current frame and difference from the previous one
			int currentArea = countNonZero( motionMask );
			int slope = currentArea - previousArea;

			if (currentArea > minAreaThreshold)
			{
				// Get all possible local maximums
				//              *
				// 	      *           *
				// previous slope | current frame
				// Positive slope = ascending; negative slope = descending
				// --> the maximum is located at the previous frame
				if ((currentArea < previousArea) && (previousSlope > 0))
				{
					cout << "Local max at " << frameNumber - 1 << endl;
					maxPosFrame.push_back(frameNumber - 1);
					maxArea.push_back(previousArea);
				}

				// Get all possible local minimums
				if ((currentArea > previousArea) && (previousSlope < 0))
				{
					minPosFrame.push_back(frameNumber - 1);
					minArea.push_back(previousArea);
					cout << "Local min at " << frameNumber - 1 << endl;
				}
			}

			// Updating previous data
			previousArea = currentArea;
			previousSlope = slope;
			previousHead = frame(largeHeadBox);
            Mat previousMask= Mat::zeros(previousHead.rows,previousHead.cols, CV_8UC1);
            //previousMask = HeadForegroundExtraction( previousHead );
            //removeColorBackground( previousHead, previousMask );
            //removeUniformColorBackground(previousHead, previousMask);
            //removeUniformColorBackground_D2(previousHead, previousMask);
            removeUniformColorBackground_D3(previousHead, previousMask);
            //removeUniformColorBackground_Q3(previousHead, previousMask);

			// Save current head masks
			extHeads.headImages.push_back(previousMask);

			// For tests
			//myfile << frameNumber << "," << currentArea << endl;
		}

		//cout << frameNumber << endl;
	}

	// ------------------ GET KEY FRAME NUMBERS ----------------------
	// For tests
	//myfile.close();

	// -- Get the minimum of local minimums --> frontal head
	if (!minArea.empty())
	{
		posMinMin = distance( minArea.begin(), min_element( minArea.begin(), minArea.end() ) );
		cout << "posMinArea = " << minPosFrame.at(posMinMin) << endl;
	}
	else
	{
		// TODO: another idea to get the frontal head? --> the first image! but doesn't help to separate profiles
		cerr << "detectHeadProfiles(): Not possible to get the frontal head frame number" << endl;
		exit(EXIT_FAILURE);
	}

	// -- Get the frame numbers of head profiles
	if (maxArea.size() >= 2)
	{
		// For tests: display areas
		for (unsigned int k = 0; k < maxArea.size(); k++)
		{
			cout << "Max area # " << maxPosFrame.at(k) << " = " << maxArea.at(k) << endl;
		}

		// -- Get frame numbers corresponding to left and right profiles
		extHeads.headFrameNumbers = getPosProfiles( maxPosFrame, maxArea, minPosFrame.at(posMinMin) );

		// Add the positions of frontal views
		extHeads.headFrameNumbers.at(FRONT1) = minPosFrame.at(posMinMin);
		extHeads.headFrameNumbers.push_back(startEnd[1]-1);	// Assumption: the face at the end is in front of the camera

		// -- Choose corresponding head images
		for (unsigned int index = 0; index < extHeads.headFrameNumbers.size(); index++ )
		{
			int posVector;

			if ((frameNumberFace - 1) <= startEnd[0])
			{
				posVector = extHeads.headFrameNumbers.at(index) - startEnd[0] - 1;
			}
			else
			{
				posVector = extHeads.headFrameNumbers.at(index) - frameNumberFace;
			}

			if (posVector < 0) 
			{
				posVector = 0;
			}

			extHeads.headFrameNumbers.at(index) = posVector;
		}
	}
	else
	{
		cerr << "detectHeadProfiles(): Not possible to get the frame numbers of head profiles" << endl; // another idea?
		exit(EXIT_FAILURE);
	}
}

void sortProfiles( const std::vector<cv::Mat>& heads, std::vector<int>& headFrameNumbers )
{
	// -- Get images of head profiles
	Mat headProfile1 = heads.at(headFrameNumbers.at(LEFT));
	Mat headProfile2 = heads.at(headFrameNumbers.at(RIGHT));

	// -- Identify the type of head profile (left or right)
	// The initial order of of head frame numbers is: left or right profiles, and frontal views
	// The resulting order is: left profile, next frontal view, right profile, next frontal view
	// Idea: if necessary, create a struct with image, frame number and profile type
	int idProfile = identifyProfiles( headProfile1, headProfile2, headFrameNumbers );

	// For tests
	if (idProfile == 1)			// Order in frame numbers is the same
	{
		imwrite(IMG_PATH + "left" + std::to_string(headFrameNumbers.at(LEFT)) + ".png", headProfile1);
		imwrite(IMG_PATH + "right" + std::to_string(headFrameNumbers.at(RIGHT)) + ".png", headProfile2);
	}
	else if (idProfile == 2)	// Order in frame numbers has changed
	{
		imwrite(IMG_PATH + "left" + std::to_string(headFrameNumbers.at(LEFT)) + ".png", headProfile2);
		imwrite(IMG_PATH + "right" + std::to_string(headFrameNumbers.at(RIGHT)) + ".png", headProfile1);
	}

	// TODO: add imshow frontal views
}

void getSegmentedHeads( const std::vector<cv::Mat>& originalHeads, const std::vector<cv::Mat>& headMasks, std::vector<cv::Mat>& finalHeads )
{

	unsigned int smallerSize = (originalHeads.size() > headMasks.size()) ? headMasks.size() : originalHeads.size();

	for (unsigned int i = 0; i < smallerSize; i++)
	{
		finalHeads.push_back(getHeadImg( originalHeads.at(i), headMasks.at(i), i, "" ));
	}

}

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
}


void removeUniformColorBackground(Mat& colorImageInput, Mat& binaryFaceMask)
{
	//DETERMIN THE INPUT IMAGE IS LANDSCAPE OR PORTRAIT
	bool isLandscape;
	if (colorImageInput.rows > colorImageInput.cols){
		isLandscape = false;
	}
	else{
		isLandscape = true;
	}
	
	Mat colorImageInputHSV;
	cvtColor(colorImageInput, colorImageInputHSV, COLOR_BGR2HSV);

	//MEAN FILTER TO REMOVE THE NOISE IN THE IMAGE
	Mat blurredColorImage;
	blur(colorImageInput, blurredColorImage,Size(4,4),Point(-1,-1));
	imshow( "blurredColorImage", blurredColorImage);
    cout << "@11111111" << endl;
	//CONVER BGR TO HSV
	Mat blurredColorImageHSV;
	cvtColor(blurredColorImage,blurredColorImageHSV,COLOR_BGR2HSV);
	
	Mat BGInterestLeft,BGInterestRight;
	if(isLandscape){
		BGInterestLeft = blurredColorImageHSV(Range(0,blurredColorImage.rows/2),Range(0,blurredColorImage.cols/3));
		BGInterestRight = blurredColorImageHSV(Range(0,blurredColorImage.rows/2),Range(blurredColorImage.cols-blurredColorImage.cols/3,blurredColorImage.cols));
	}
	else{
		BGInterestLeft = blurredColorImage(Range(0,blurredColorImage.rows/2),Range(0,blurredColorImage.cols/5));
		BGInterestRight = blurredColorImage(Range(0,blurredColorImage.rows/2),Range(blurredColorImage.cols-blurredColorImage.cols/5,blurredColorImage.cols));
	}
    cout << "@222222222" << endl;
	imshow( "BGInterestLeft", BGInterestLeft);
	imshow( "BGInterestRight", BGInterestRight);
	
	int x,y;
	int count = 0;
	std::vector<int> sum; sum.resize(3,0);
	std::vector<int> sumSqure;  sumSqure.resize(3,0);
	std::vector<double> mean; mean.resize(3, 0);
	std::vector<double> stddev;  stddev.resize(3, 0);
	for (y = 0;y < BGInterestLeft.rows;y++)
	{
		for (x = 0;x < BGInterestLeft.cols;x++)
		{
			sum[0] = sum[0] + BGInterestLeft.at<Vec3b>(y, x).val[0];
			sum[1] = sum[1] + BGInterestLeft.at<Vec3b>(y, x).val[1];
			sum[2] = sum[2] + BGInterestLeft.at<Vec3b>(y, x).val[2];
			sumSqure[0] = sumSqure[0] + BGInterestLeft.at<Vec3b>(y, x).val[0] * BGInterestLeft.at<Vec3b>(y, x).val[0];
			sumSqure[1] = sumSqure[1] + BGInterestLeft.at<Vec3b>(y, x).val[1] * BGInterestLeft.at<Vec3b>(y, x).val[1];
			sumSqure[2] = sumSqure[2] + BGInterestLeft.at<Vec3b>(y, x).val[2] * BGInterestLeft.at<Vec3b>(y, x).val[2];
			count++;
		}
	} 
    cout << "@3333333333" << endl;
	for (y = 0; y < BGInterestRight.rows; y++)
	{
		for (x = 0; x < BGInterestRight.cols; x++)
		{
			sum[0] = sum[0] + BGInterestRight.at<Vec3b>(y, x).val[0];
			sum[1] = sum[1] + BGInterestRight.at<Vec3b>(y, x).val[1];
			sum[2] = sum[2] + BGInterestRight.at<Vec3b>(y, x).val[2];
			sumSqure[0] = sumSqure[0] + BGInterestRight.at<Vec3b>(y, x).val[0] * BGInterestRight.at<Vec3b>(y, x).val[0];
			sumSqure[1] = sumSqure[1] + BGInterestRight.at<Vec3b>(y, x).val[1] * BGInterestRight.at<Vec3b>(y, x).val[1];
			sumSqure[2] = sumSqure[2] + BGInterestRight.at<Vec3b>(y, x).val[2] * BGInterestRight.at<Vec3b>(y, x).val[2];
			count++;
		}
	}
    cout << "@44444444444" << endl;
	mean[0] = (double)(sum[0] / std::max(count, 1));
	mean[1] = (double)(sum[1] / std::max(count, 1));
	mean[2] = (double)(sum[2] / std::max(count, 1));
	stddev[0] = sqrt((double)((sumSqure[0] / std::max(count, 1)) - (mean[0] * mean[0])));
	stddev[1] = sqrt((double)((sumSqure[1] / std::max(count, 1)) - (mean[1] * mean[1])));
	stddev[2] = sqrt((double)((sumSqure[2] / std::max(count, 1)) - (mean[2] * mean[2])));


	//Generate binary mask of image with background removed
	Vec3b HSVPixel;
    Mat binaryFaceMask1(colorImageInput.rows,colorImageInput.cols,CV_8UC1);
    binaryFaceMask1.empty();
	for (int y = 0; y<colorImageInput.rows; y++)
	{
		for (int x = 0; x<colorImageInput.cols; x++)
		{
			//BGRPixel = colorImageInput.at<Vec3b>(y, x);
			HSVPixel = colorImageInputHSV.at<Vec3b>(y, x);
			
			//if (abs(HSVPixel.val[0] - mean[0]) <= (2 * stddev[0]) && abs(HSVPixel.val[1] - mean[1]) <= (2 * stddev[1]) && abs(HSVPixel.val[2] - mean[2]) <= (2*stddev[2]))
			if (abs(HSVPixel.val[1] - mean[1]) <= (2 * stddev[1]) && abs(HSVPixel.val[2] - mean[2]) <= (2 * stddev[2]))
			{
				binaryFaceMask1.at<uchar>(y, x) = 0;
			}
			else
			{
				binaryFaceMask1.at<uchar>(y, x) = 255;
			}
		}
	}
    cout << "@5555555555" << endl;
    Mat imLabel(colorImageInput.rows,colorImageInput.cols,CV_16UC1);
    imLabel.empty();
    ImLabelFlatZones_WithCriterion(binaryFaceMask1, binaryFaceMask1, Area, imLabel);
    imshow("lable",imLabel);
    cout << "@66666666666" << endl;
	std::pair<int,int> minmax = MeasMinMax(imLabel);

    cout << "@77777777777" << endl;
    //cvtColor(colorImageInput,binaryFaceMask,CV_RGB2GRAY);

    //binaryFaceMask.(colorImageInput.rows,colorImageInput.cols,CV_8U);
    for (int y = 0; y < colorImageInput.rows; y++)
	{
        for (int x = 0; x < colorImageInput.cols; x++)
		{
			if (imLabel.at<int>(y, x) == minmax.second)
			{
				binaryFaceMask.at<uchar>(y, x) = 255;
			}
			else
			{
				binaryFaceMask.at<uchar>(y, x) = 0;
			}
		}
	}

    cout << "@888888888888" << endl;

    //imshow("binaryFaceMask", binaryFaceMask);
}





void removeUniformColorBackground_D2(Mat& colorImageInput, Mat& binaryFaceMask)
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

    VFactor = 20.0;

//    cout<<"Mean Hue Value: " << float(meanHValue)<<endl;
//    cout<<"Mean Saturation Value: " << float(meanSValue)<<endl;
//    cout<<"Mean V Value: " << float(meanVValue)<<endl;
//    cout<<"Hue stddev: " << float(HValueStdDev)<<endl;
//    cout<<"Saturation stddev: " << float(SValueStdDev)<<endl;
//    cout<<"Values stddev: " << float(VValueStdDev)<<endl;

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
    kernelSize = double(binaryFaceMask.rows)*0.02;
	if (kernelSize%2 != 0 )
	   kernelSize = kernelSize+1;
   

    //Display kernel size
    cout<<"Kernel Size: "<<kernelSize<<endl;

    //Create structuring element
    element = getStructuringElement(MORPH_ELLIPSE,Size(kernelSize,kernelSize),Point(-1,-1));

    //Perform morphological closing on mask
    morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_OPEN, element ,Point(-1,-1),1,0,morphologyDefaultBorderValue());
    morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_CLOSE, element ,Point(-1,-1),1,0,morphologyDefaultBorderValue());

    /*
    Mat imLabel;
    //ImLabelFlatZones_WithCriterion(c, binaryFaceMask, Area, imLabel);
    ImLabel_QG(binaryFaceMask,1,imLabel);
    imshow("lable",imLabel);
    cout << "@66666666666" << endl;


    std::vector<int> histo;histo.resize(3000,0);
    for (int y = 0; y <imLabel.rows;y++)
    {
        for(int x = 0; x < imLabel.cols;x++)
        {
          histo[imLabel.at<int>(y,x)]++;
        }

    }
    int maximum = 0;
    int index = 0;
    for(int y = 0; y<histo.size();y++)
    {
        if(histo[y]>maximum)
            index = y;
    }

    //std::pair<int,int> minmax = MeasMinMax(imLabel);


    for (int y = 0; y < binaryFaceMask.rows; y++)
    {
        for (int x = 0; x < binaryFaceMask.cols; x++)
        {
            //if (imLabel.at<int>(y, x) == minmax.second)
            if (imLabel.at<int>(y,x) == index)
            {
                binaryFaceMask.at<uchar>(y, x) = 255;
            }
            else
            {
                binaryFaceMask.at<uchar>(y, x) = 0;
            }
        }
    }

*/



}




void ImLabelFlatZones_WithCriterion(const Mat & imFZ, const Mat & imData, const ZoneInfoType criterion, Mat & imLabel)
{
	int w = imFZ.cols, h = imFZ.rows;
	imLabel.empty();

	int i = 0, j = 0, ii = 0, jj = 0;
    Mat flag(h, w, CV_8UC1);

	flag.empty();


	int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
	int rayon = 1;
	int current_label = 0;
	int counter = 1;

	int maxival = std::numeric_limits<int>::max();
	int current_value = 0;

	std::list<std::pair<int, int> > regionpoints; if (!regionpoints.empty()) printf("The list regionpoints is not empty\n");
	std::list<std::pair<int, int> > neighborpoints; if (!neighborpoints.empty()) printf("The list neighborpoints is not empty\n");

	FZInfoNode<int> fzinode;

	std::pair<int, int> temppaire(0, 0);
	std::pair<int, int> local(0, 0);


	for (i = 0; i<w; i++){
		for (j = 0; j<h; j++){
			if (imFZ.at<uchar>(j, i) != 0){
				if (flag.at<uchar>(j, i) == 0){					
					flag.at<uchar>(j, i) = 1;
					current_label = imFZ.at<uchar>(j, i);
					fzinode.regionnumber = counter;
					counter++;

					temppaire.first = i; temppaire.second = j;
					regionpoints.push_back(temppaire);
					neighborpoints.push_back(temppaire);

					// 					imLabel(i,j) = fzinode.regionnumber;

					current_value = imData.at<uchar>(j, i);
					IsMinOrMax(i, fzinode.x_min, fzinode.x_max);
					IsMinOrMax(j, fzinode.y_min, fzinode.y_max);
					IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
					fzinode.area = fzinode.area + 1;
					fzinode.average += current_value;
					fzinode.volume += (long int)(maxival - current_value);

					while (!neighborpoints.empty()){
						local = neighborpoints.front();
						neighborpoints.pop_front();

						borneInfx = std::max(0, local.first - rayon);
						borneInfy = std::max(0, local.second - rayon);
						borneMaxx = std::min(w, local.first + rayon + 1);
						borneMaxy = std::min(h, local.second + rayon + 1);

						for (ii = borneInfx; ii<borneMaxx; ii++){
							for (jj = borneInfy; jj<borneMaxy; jj++){
								if ((imFZ.at<uchar>(jj, ii) == current_label) && (flag.at<uchar>(jj, ii) == 0)){
									if ((ii == i) && (jj == j)) printf("Pb: iteration sur le point central\n");

									temppaire.first = ii; temppaire.second = jj;
									regionpoints.push_back(temppaire);
									neighborpoints.push_back(temppaire);

									// 									imLabel(ii,jj) = fzinode.regionnumber;
									flag.at<uchar>(jj, ii) = 1;
									current_value = imData.at<uchar>(jj, ii);
									IsMinOrMax(ii, fzinode.x_min, fzinode.x_max);
									IsMinOrMax(jj, fzinode.y_min, fzinode.y_max);
									IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
									fzinode.area = fzinode.area + 1;
									fzinode.average += current_value;
									fzinode.stddev += double(current_value * current_value);
									fzinode.volume += (long int)(maxival - current_value);
								}
							}
						}
					}
					fzinode.average = double(fzinode.average) / double(fzinode.area);
					fzinode.stddev = std::sqrt(double(fzinode.stddev) / double(fzinode.area) - fzinode.average * fzinode.average);
					fzinode.dynamic = int((fzinode.val_max - fzinode.val_min));
					fzinode.volume = (long int)(fzinode.volume - fzinode.area*(maxival - fzinode.val_max));

					std::list<std::pair<int, int> >::iterator ittemp = regionpoints.begin();
					std::list<std::pair<int, int> >::iterator itend = regionpoints.end();

					switch (criterion)
					{
					case Width:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.x_max - fzinode.x_min);
						break;
					case Height:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.y_max - fzinode.y_min);
						break;
					case Val_min:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.val_min);
						break;
					case Average:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.average);
						break;
					case Dynamic:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.dynamic);
						break;
					case Area:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.area);
						break;
					case Volume:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.volume);
						break;
					case StdDev:
						for (; ittemp != itend; ittemp++)
							imLabel.at<int>(ittemp->second, ittemp->first) = int(fzinode.stddev);
						break;
					default:
						printf("Sorry, this feature is not available yet\n");
					}
					regionpoints.clear();
					neighborpoints.clear();
					fzinode.Reset();
				}
			}
		}
	}
}


void removeUniformColorBackground_D3(Mat& colorImageInput, Mat& binaryFaceMask)
{
	//Variable Declarations
	uchar meanRValue, meanGValue, meanBValue, RValueStdDev, GValueStdDev, BValueStdDev;
	uchar meanHValue, meanSValue, meanVValue;
	double HValueSum = 0, SValueSum = 0, VValueSum = 0;
    double HValueStdDev, SValueStdDev, VValueStdDev;
	double HValueVariance, SValueVariance, VValueVariance;
	Vec3b refPixel1, refPixel2, refPixel3, refPixel4, BGRPixel;
	Vec3b HSVRefPixel1, HSVRefPixel2, HSVRefPixel3, HSVRefPixel4, HSVRefPixel5, HSVRefPixel6, HSVPixel;
	int cornerSideLength = 0, windowWidth = 0, windowHeight = 0, numberHSVValues = 0;
	Mat blurredColorImageInput, HSVColorImageInput, HSVHistogram, HSVBlurredColorImageInput;
	int hBins = 500, sBins = 600;
	int channels[] = { 0, 1 };
	int histSizes[] = { hBins, sBins };
	float hueRange[] = { 0, 179 }, saturationRange[] = { 0, 255 };
	const float* ranges[] = { hueRange, saturationRange };
	float HFactor, SFactor, VFactor;
	int counter = 0;
	Mat element;
	int kernelSize = 0;

	// Declare variable 2D array to store HSV background pixel reference values
	uchar** refPixelHSVValues;

	//Convert image to HSV format
	cvtColor(colorImageInput, HSVColorImageInput, COLOR_BGR2HSV);

	// Calculate Histogram for HSV image
	calcHist(&HSVColorImageInput, 1, channels, Mat(), HSVHistogram, 2, histSizes, ranges, true, false);
	normalize(HSVHistogram, HSVHistogram, 0, 255, NORM_MINMAX, -1, Mat());

	// Extract pixels from two windows located in the top corners (2 corners) of the image (1/16 x 1/16) 
	windowWidth = HSVColorImageInput.cols / 16;
    windowHeight = HSVColorImageInput.rows / 2;

	//Total number of pixels used to evaluate background
	numberHSVValues = windowWidth*windowHeight * 2;

	//Allocate memory for each pixel to be read for background evaluation
	refPixelHSVValues = new uchar*[numberHSVValues];
	for (int i = 0; i < numberHSVValues; i++)
	{
		refPixelHSVValues[i] = new uchar[3];
	}

	//Read each pixel used for background evaluation
	//Read pixels on left side of image
	for (int x = 0; x<windowWidth; x++)
	{
		for (int y = 0; y<windowHeight; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
			refPixelHSVValues[counter][0] = HSVPixel.val[0];
			refPixelHSVValues[counter][1] = HSVPixel.val[1];
			refPixelHSVValues[counter][2] = HSVPixel.val[2];
			counter++;
		}
	}

	//Read pixels on right side of image
	for (int x = HSVColorImageInput.cols - 1; x>HSVColorImageInput.cols - 1 - windowWidth; x--)
	{
		for (int y = 0; y<windowHeight; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
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
	meanHValue = HValueSum / (numberHSVValues);
	meanSValue = SValueSum / (numberHSVValues);
	meanVValue = VValueSum / (numberHSVValues);

	//calculate variance of H,S,V values
	HValueVariance = 0;
	SValueVariance = 0;
	VValueVariance = 0;
	for (int i = 0; i < numberHSVValues; i++) //*************Major modification made here 2014-11-23 17h39**************
	{
		HValueVariance = HValueVariance + pow(double(refPixelHSVValues[i][0]) - double(meanHValue), 2.0);
		SValueVariance = SValueVariance + pow(double(refPixelHSVValues[i][1]) - double(meanSValue), 2.0);
		VValueVariance = VValueVariance + pow(double(refPixelHSVValues[i][2]) - double(meanVValue), 2.0);
	}
	HValueVariance = HValueVariance / (numberHSVValues);
	SValueVariance = SValueVariance / (numberHSVValues);
	VValueVariance = VValueVariance / (numberHSVValues);

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
	cvtColor(colorImageInput, binaryFaceMask, CV_RGB2GRAY);

	//Filtering factors (must be modified for each image!)
	if (HValueStdDev <2.0)
		HFactor = 18.0;
	else
	if (HValueStdDev <5.0)
		HFactor = 8.0;
	else
        HFactor = 3.0;
	if (SValueStdDev <2.0)
		SFactor = 22.0;
	else
	if (SValueStdDev <5.0)
		SFactor = 8.0;
	else
        SFactor = 3.0;

	VFactor = 1000.0;

//	cout << "Mean Hue Value: " << float(meanHValue) << endl;
//	cout << "Mean Saturation Value: " << float(meanSValue) << endl;
//	cout << "Mean V Value: " << float(meanVValue) << endl;
//	cout << "Hue stddev: " << float(HValueStdDev) << endl;
//	cout << "Saturation stddev: " << float(SValueStdDev) << endl;
//	cout << "Values stddev: " << float(VValueStdDev) << endl;

	//Generate binary mask of image with background removed
	for (int y = 0; y<colorImageInput.rows; y++)
	{
		for (int x = 0; x<colorImageInput.cols; x++)
		{
			BGRPixel = colorImageInput.at<Vec3b>(y, x);
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
			//if (abs(BGRPixel.val[0]-meanBValue)<(6*BValueStdDev) && abs(BGRPixel.val[1]-meanGValue)<(6*GValueStdDev) && abs(BGRPixel.val[2]-meanRValue)<(6*RValueStdDev))
			if (abs(HSVPixel.val[0] - meanHValue) <= (HFactor*HValueStdDev) && abs(HSVPixel.val[1] - meanSValue) <= (SFactor*SValueStdDev)) //&& abs(HSVPixel.val[2]-meanVValue)<=(VFactor*VValueStdDev))
			{
				binaryFaceMask.at<uchar>(y, x) = 0;
			}
			else
			{
				binaryFaceMask.at<uchar>(y, x) = 255;
			}
		}
	}
	//Display result Mask before morphological "closing"
	//namedWindow( "Result Mask before morphological closing", WINDOW_NORMAL );// Create a window for display.
	//imshow( "Result Mask before morphological closing", binaryFaceMask ); 

	//Perform morphological closing operation on mask
    kernelSize = double(binaryFaceMask.rows)*0.03;
	if (kernelSize % 2 != 0)
		kernelSize = kernelSize + 1;

	//Display kernel size
	cout << "Kernel Size: " << kernelSize << endl;

	//Create structuring element
	element = getStructuringElement(MORPH_ELLIPSE, Size(kernelSize, kernelSize), Point(-1, -1));

	//Perform morphological closing on mask
    morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_OPEN, element, Point(-1, -1), 1, 0, morphologyDefaultBorderValue());
    //element = getStructuringElement(MORPH_ELLIPSE, Size(2*kernelSize, 2*kernelSize), Point(-1, -1));

    morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_CLOSE, element, Point(-1, -1), 1, 0, morphologyDefaultBorderValue());

/*
    Mat imLabel(colorImageInput.rows,colorImageInput.cols,CV_16UC1);
    imLabel.empty();
    ImLabelFlatZones_WithCriterion(binaryFaceMask, binaryFaceMask, Area, imLabel);
    imshow("lable",imLabel);
    binaryFaceMask.empty();

    cout << "@66666666666" << endl;
    std::pair<int,int> minmax = MeasMinMax(imLabel);

    cout << "@77777777777" << endl;

    //binaryFaceMask.(colorImageInput.rows,colorImageInput.cols,CV_8U);
    for (int y = 0; y < colorImageInput.rows; y++)
    {
        for (int x = 0; x < colorImageInput.cols; x++)
        {
            if (imLabel.at<int>(y, x) == minmax.second)
            {
                binaryFaceMask.at<uchar>(y, x) = 255;
            }
            else
            {
                binaryFaceMask.at<uchar>(y, x) = 0;
            }
        }
    }

*/

}


void removeUniformColorBackground_Q3(Mat& colorImageInput, Mat& binaryFaceMask)
{
	//Variable Declarations
	uchar meanRValue, meanGValue, meanBValue;
	double RValueStdDev, GValueStdDev, BValueStdDev; 
	uchar meanHValue, meanSValue, meanVValue;
	int HValueSum = 0, SValueSum = 0, VValueSum = 0;
	double HValueStdDev, SValueStdDev, VValueStdDev;
	double HValueVariance, SValueVariance, VValueVariance;
	//Vec3b refPixel1, refPixel2, refPixel3, refPixel4, BGRPixel;
	Vec3b BGRPixel;
	//Vec3b HSVRefPixel1, HSVRefPixel2, HSVRefPixel3, HSVRefPixel4, HSVRefPixel5, HSVRefPixel6, HSVPixel;
	Vec3b HSVPixel;
	int cornerSideLength = 0, windowWidth = 0, windowHeight = 0, numberHSVValues = 0;
	Mat blurredColorImageInput, HSVColorImageInput, HSVHistogram, HSVBlurredColorImageInput;
	int hBins = 500, sBins = 600;
	int channels[] = { 0, 1 };
	int histSizes[] = { hBins, sBins };
	float hueRange[] = { 0, 179 }, saturationRange[] = { 0, 255 };
	const float* ranges[] = { hueRange, saturationRange };
	float HFactor, SFactor, VFactor;
	int counter = 0;
	Mat element;
	int kernelSize = 0;

	// Declare variable 2D array to store HSV background pixel reference values
	uchar** refPixelHSVValues;

	//Convert image to HSV format
	cvtColor(colorImageInput, HSVColorImageInput, COLOR_BGR2HSV);

	// Calculate Histogram for HSV image
	calcHist(&HSVColorImageInput, 1, channels, Mat(), HSVHistogram, 2, histSizes, ranges, true, false);
	normalize(HSVHistogram, HSVHistogram, 0, 255, NORM_MINMAX, -1, Mat());

	// Extract pixels from two windows located in the top corners (2 corners) of the image (1/16 x 1/16) 
	windowWidth = HSVColorImageInput.cols / 16;
	windowHeight = HSVColorImageInput.rows / 2;

	//Total number of pixels used to evaluate background
	numberHSVValues = windowWidth*windowHeight * 2;

	//Allocate memory for each pixel to be read for background evaluation
	refPixelHSVValues = new uchar*[numberHSVValues];
	for (int i = 0; i < numberHSVValues; i++)
	{
		refPixelHSVValues[i] = new uchar[3];
	}

	//Read each pixel used for background evaluation
	//Read pixels on left side of image
	for (int x = 0; x<windowWidth; x++)
	{
		for (int y = 0; y<windowHeight; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
			refPixelHSVValues[counter][0] = HSVPixel.val[0];
			refPixelHSVValues[counter][1] = HSVPixel.val[1];
			refPixelHSVValues[counter][2] = HSVPixel.val[2];
			counter++;
		}
	}

	//Read pixels on right side of image
	for (int x = HSVColorImageInput.cols - 1; x>HSVColorImageInput.cols - 1 - windowWidth; x--)
	{
		for (int y = 0; y<windowHeight; y++)
		{
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
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
	meanHValue = HValueSum / (numberHSVValues);
	meanSValue = SValueSum / (numberHSVValues);
	meanVValue = VValueSum / (numberHSVValues);

	//calculate variance of H,S,V values
	HValueVariance = 0;
	SValueVariance = 0;
	VValueVariance = 0;
	for (int i = 0; i < numberHSVValues; i++) //*************Major modification made here 2014-11-23 17h39**************
	{
		HValueVariance = HValueVariance + pow(double(refPixelHSVValues[i][0]) - double(meanHValue), 2.0);
		SValueVariance = SValueVariance + pow(double(refPixelHSVValues[i][1]) - double(meanSValue), 2.0);
		VValueVariance = VValueVariance + pow(double(refPixelHSVValues[i][2]) - double(meanVValue), 2.0);
	}
	HValueVariance = HValueVariance / (numberHSVValues);
	SValueVariance = SValueVariance / (numberHSVValues);
	VValueVariance = VValueVariance / (numberHSVValues);

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
	cvtColor(colorImageInput, binaryFaceMask, CV_RGB2GRAY);

	//Filtering factors (must be modified for each image!)
	if (HValueStdDev <2.0)
		HFactor = 18.0;
	else
	if (HValueStdDev <5.0)
		HFactor = 8.0;
	else
		HFactor = 4.0;
	if (SValueStdDev <2.0)
		SFactor = 22.0;
	else
	if (SValueStdDev <5.0)
		SFactor = 8.0;
	else
		SFactor = 4.0;

	VFactor = 1000.0;

	cout << "Mean Hue Value: " << float(meanHValue) << endl;
	cout << "Mean Saturation Value: " << float(meanSValue) << endl;
	cout << "Mean V Value: " << float(meanVValue) << endl;
	cout << "Hue stddev: " << float(HValueStdDev) << endl;
	cout << "Saturation stddev: " << float(SValueStdDev) << endl;
	cout << "Values stddev: " << float(VValueStdDev) << endl;

	//Generate binary mask of image with background removed
	for (int y = 0; y<colorImageInput.rows; y++)
	{
		for (int x = 0; x<colorImageInput.cols; x++)
		{
			BGRPixel = colorImageInput.at<Vec3b>(y, x);
			HSVPixel = HSVColorImageInput.at<Vec3b>(y, x);
			//if (abs(BGRPixel.val[0]-meanBValue)<(6*BValueStdDev) && abs(BGRPixel.val[1]-meanGValue)<(6*GValueStdDev) && abs(BGRPixel.val[2]-meanRValue)<(6*RValueStdDev))
			if (abs(HSVPixel.val[0] - meanHValue) <= (HFactor*HValueStdDev) && abs(HSVPixel.val[1] - meanSValue) <= (SFactor*SValueStdDev)) //&& abs(HSVPixel.val[2]-meanVValue)<=(VFactor*VValueStdDev))
			{
				binaryFaceMask.at<uchar>(y, x) = 0;
			}
			else
			{
				binaryFaceMask.at<uchar>(y, x) = 255;
			}
		}
	}
	//Display result Mask before morphological "closing"
	//namedWindow( "Result Mask before morphological closing", WINDOW_NORMAL );// Create a window for display.
	//imshow( "Result Mask before morphological closing", binaryFaceMask ); 
	
	Mat imLabel(colorImageInput.rows,colorImageInput.cols,CV_16UC1);
    imLabel.empty();
    ImLabelFlatZones_WithCriterion(binaryFaceMask, binaryFaceMask, Area, imLabel);
    imshow("lable",imLabel);
    binaryFaceMask.empty();
    
    cout << "@66666666666" << endl;
	std::pair<int,int> minmax = MeasMinMax(imLabel);

    cout << "@77777777777" << endl;   

    //binaryFaceMask.(colorImageInput.rows,colorImageInput.cols,CV_8U);
    for (int y = 0; y < colorImageInput.rows; y++)
	{
        for (int x = 0; x < colorImageInput.cols; x++)
		{
			if (imLabel.at<int>(y, x) == minmax.second)
			{
				binaryFaceMask.at<uchar>(y, x) = 255;
			}
			else
			{
				binaryFaceMask.at<uchar>(y, x) = 0;
			}
		}
	}

	//Perform morphological closing operation on mask
	kernelSize = double(binaryFaceMask.rows)*0.02;
	if (kernelSize % 2 != 0)
		kernelSize = kernelSize + 1;

	//Display kernel size
    //cout << "Kernel Size: " << kernelSize << endl;

	//Create structuring element
	element = getStructuringElement(MORPH_ELLIPSE, Size(kernelSize, kernelSize), Point(-1, -1));

	//Perform morphological closing on mask
	morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_OPEN, element, Point(-1, -1), 1, 0, morphologyDefaultBorderValue());
	morphologyEx(binaryFaceMask, binaryFaceMask, MORPH_CLOSE, element, Point(-1, -1), 1, 0, morphologyDefaultBorderValue());
}

