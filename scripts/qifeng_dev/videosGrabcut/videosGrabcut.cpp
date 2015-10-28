
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream> // for standard I/O
#include <string>   // for strings
#include <iomanip>  // for controlling float print precision
#include <sstream>  // string to number conversion
#include <boost/lexical_cast.hpp>



#include "../../../include/Include.hpp"


using namespace std;
using namespace cv;
using namespace boost::filesystem;

void rotate(const cv::Mat& image, int degrees, cv::Mat & imOut )
{
        cv::Mat res;
        switch (degrees % 360) {
            case 0:
                res = image;
                break;
            case 90:
                res = image.t();
                cv::flip(res, res, 1);
                break;
            case 180:
                cv::flip(image, res, -1);
                break;
            case 270:
                res = image.t();
                cv::flip(res, res, 0);
                break;
            default:
                cv::Mat r = cv::getRotationMatrix2D({image.cols/2.0F, image.rows/2.0F}, degrees, 1.0);
                int len = std::max(image.cols, image.rows);
                cv::warpAffine(image, res, r, cv::Size(len, len));
                break; //image size will change
        }
        imOut = res;
    }




// Function main
int main(int argc, char * argv[])
{
	if(argc!=3){
			std::cout << "Usage: " << argv[0] << " infile  outfile"<<std::endl;
	}
		
	path dir_path(argv[1]);
	path dir_path_out(argv[2]);
	path dir_path_model(argv[3]);
	
	path im_name;	
	path im_path;
	path result_path;
	path intermediate_path;
	
	boost::filesystem::directory_iterator end_itr;
	for ( boost::filesystem::directory_iterator itr( dir_path ); itr != end_itr; ++itr){
		if ( !is_directory(itr->status()) )
		{
			im_name = itr->path().filename();
			im_path = dir_path/im_name;
			result_path = dir_path_out/im_name.stem();
			std::cout << im_path.string().c_str() << std::endl;
			std::vector<ParamT> param; param.clear();
			std::string path_model = dir_path_model.string();
			
			intermediate_path = dir_path_out/im_name.stem()/"/";
			if (!exists(intermediate_path))
				create_directories(intermediate_path);
			dbi.SetPath(intermediate_path);
					
			std::vector<DipColorImage<uchar> > frames;
			DipVideo videoin(im_path.string().c_str());

			videoin.decodeFrames(frames);
			int ex = static_cast<int>(videoin.get(CV_CAP_PROP_FOURCC));
			int w  = videoin.getFrameWidth();
			int h  = videoin.getFrameHeight();
			int frameNumber = frames.size();
			
			DipImage<uchar> imMask(w,h);imMask.init(0);
			DipImage<uchar> imMaskInitial(w,h);imMaskInitial.init(0);
			maskDetection(frames[0],path_model,imMask);
			ImCopy(imMask,imMaskInitial);
				
			
				

			int count = 0;
			std::vector<DipColorImage<uchar> > framesFG;
			for(int p = 1;p < frameNumber; p++){
				if(!frames[p].empty()){
					DipColorImage<uchar> imFGColor(w,h);imFGColor.init(0);					
					DipColorImage<uchar> imColorPrevious(w,h);imColorPrevious.init(0);					
					DipColorImage<uchar> imColor(w,h);imColor.init(0);					
					DipImage<uchar> imMaskAll(w,h);imMaskAll.init(255);
					
					//Add motion mask
					//~ ImColorMeanFilter(frames[p],imMaskAll,int(5),imColor);
					//~ ImColorMeanFilter(frames[p-1],imMaskAll,int(5),imColorPrevious);

					DipImage<uchar> imfgMask(w,h);imfgMask.init(0);
					imageMatting_3(frames[p],frames[p-1],frames[0],imMask,imMaskInitial,imfgMask);
					//~ imageMatting_2(imColor,imColorPrevious,imMask,imfgMask);
					dbi.WriteImage(imfgMask, "imfgMask.png");
					
					//~ ImCopy(imfgMask,imMask);
					
					//extract the largest one
					DipImage<double> imFGD(w,h); imFGD.init(0);
					DipImage<uchar> imFGLarget(w,h); imFGLarget.init(0);
					DipImage<int> imFGLabel(w,h); imFGLabel.init(0);
					ImLabelFlatZones_WithCriterion(imfgMask,imFGD,Area,imFGLabel);
					std::pair<double,double> minmax = MeasMinMax(imFGLabel,imfgMask);
					ImThreshold(imFGLabel, int(minmax.second), int(minmax.second), uchar(255),uchar(0),imFGLarget);
					dbi.WriteImage(imFGLarget, "imFGLarget.png");
					ImCopy(imFGLarget,imMask);
					
					
					
					//ImClose
					param.clear();
					param.push_back(int(w/50));
					DipImage<uchar> imFGLargetClose(w,h);imFGLargetClose.init(0);
					ImClose(imFGLarget,Fast,param,imFGLargetClose);
					dbi.WriteImage(imFGLargetClose,"imFGLargetClose.png");
					ImCopy(imFGLargetClose,imFGLarget);

																
					//fill the hole
					DipImage<uchar> imInInvert(w,h); imInInvert.init(0);
					DipImage<double> imInInvertD(w,h); imInInvertD.init(0);
					ImInvert(imFGLarget,imInInvert);
					ImCopy(imInInvert,imInInvertD);
					
					DipImage<int> imLabel(w,h); imLabel.init(0);
					ImLabelFlatZones_WithCriterion(imInInvert,imInInvertD,Area,imLabel);


					DipImage<uchar> imBG(w,h); imBG.init(0);
					DipImage<uchar> imBGFinal(w,h); imBGFinal.init(0);
					ImThreshold(imLabel, int(0),int(w*h/20), uchar(255),uchar(0),imBG);
					dbi.WriteImage(imBG,"imBG.png");
					
					ImAddImage(imFGLarget,Clip,imBG,imBGFinal);
					dbi.WriteImage(imBGFinal,"imBGFinal.png");				
				

					
					DipColorImage<uchar> imFGColorBGR(w,h);imFGColorBGR.init(0);
						
					
					ImCopy(frames[p],imBGFinal,imFGColor);
					cvtColor(imFGColor,imFGColorBGR,CV_RGB2BGR);
					framesFG.push_back(imFGColorBGR);
					dbi.WriteImage(imFGColor, "imFGColor" + boost::lexical_cast<string>(p) + ".png");
					
					count = count +1;
					
				}
			}
			
			
			const std::string NAME = result_path.string() +".avi";
			VideoWriter outputVideo; 
			Size S = Size((int) videoin.get(CV_CAP_PROP_FRAME_WIDTH),(int) videoin.get(CV_CAP_PROP_FRAME_HEIGHT));   
			
			outputVideo.open(NAME, CV_FOURCC('D','I','V','X'), videoin.get(CV_CAP_PROP_FPS), S, true);

			if (!outputVideo.isOpened())
			{
				cout  << "Could not open the output video for write: "  << endl;
				return -1;
			}
			for(int p = 1;p < frameNumber; p++){ //Show the image captured in the window and repeat
       
				if (frames[p].empty()) break;   // check if at end			
				outputVideo << framesFG[p-1];
			}
			outputVideo.release();
			
		}
		
	}


	return 0;
}
