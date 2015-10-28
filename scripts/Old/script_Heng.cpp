#include "opencv2/core/core.hpp"
#include "opencv/cv.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>

#include "../include/humanmodel.hpp"

using namespace cv;
using namespace std;


int main( int argc, char** argv )
{
    if( argc != 2)
    {
     cout <<" Usage: display_image ImageToLoadAndDisplay" << endl;
     return -1;
    }

    Mat image;
    image = imread(argv[1], CV_LOAD_IMAGE_COLOR);   // Read the file

    if(! image.data )                              // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }
    
    //For testing CountBrightPixels
    //vector<int> counts = CountBrightPixels(image, ROW);
    


    //For testing FindBlobd
    cv::Mat output = cv::Mat::zeros(image.size(), CV_8UC3);
    std::vector < std::vector<cv::Point2i > > blobs;

    //Grayscale matrix
    Mat gray_image(image.size(), CV_8U);

    //Convert RGB to Gray
    cvtColor(image, gray_image, CV_BGR2GRAY);

    //Binary image
    Mat binary_image;

    //Apply thresholding
    threshold(gray_image, binary_image, 0, 1, THRESH_BINARY);   //TODO


    FindBlobs(binary_image, blobs);
    
    // Randomy color the blobs
    for(size_t i=0; i < blobs.size(); i++) {
        unsigned char r = 255 * (rand()/(1.0 + RAND_MAX));
        unsigned char g = 255 * (rand()/(1.0 + RAND_MAX));
        unsigned char b = 255 * (rand()/(1.0 + RAND_MAX));

        for(size_t j=0; j < blobs[i].size(); j++) {
            int x = blobs[i][j].x;
            int y = blobs[i][j].y;

            output.at<cv::Vec3b>(y,x)[0] = b;
            output.at<cv::Vec3b>(y,x)[1] = g;
            output.at<cv::Vec3b>(y,x)[2] = r;
        }
    }

    imshow("output image", output);
    

    waitKey(0);                                          // Wait for a keystroke in the window
    return 0;
}



