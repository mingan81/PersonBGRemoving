
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream> // for standard I/O
#include <string>   // for strings
#include <iomanip>  // for controlling float print precision
#include <sstream>  // string to number conversion
#include <boost/lexical_cast.hpp>



#include "../include/Include.hpp"



using namespace std;
using namespace cv;
using namespace boost::filesystem;



// Function main
int main(int argc, char * argv[])
{
	if(argc!=4){
			std::cout << "Usage: " << argv[0] << " infile pathmodel outfile"<<std::endl;
	}
	
	
	path dir_path(argv[1]);
	path dir_path_out(argv[3]);
	path dir_path_model(argv[2]);
	
	
	std::string path_model = dir_path_model.string();

	std::cout << dir_path.string().c_str() << std::endl;
	std::vector<ParamT> param; param.clear();
			
	std::vector<DipColorImage<uchar> > frames;
	DipVideo videoin(dir_path.string().c_str());

	videoin.decodeFrames(frames);
	int ex = static_cast<int>(videoin.get(CV_CAP_PROP_FOURCC));
	int w  = videoin.getFrameWidth();
	int h  = videoin.getFrameHeight();
	int frameNumber = frames.size();
			
	
	
	DipImage<uchar> imMask(w,h);imMask.init(0);
	maskDetection(frames[0],path_model,imMask);
	
	int count = 0;
	std::vector<DipColorImage<uchar> > framesFG;
	for(int p = 1;p < frameNumber; p++){
		if(!frames[p].empty()){
			std::cout << p << endl;
			DipColorImage<uchar> imFGColor(w,h);imFGColor.init(0);					
			
			//@todo Add motion mask
			//~ ImMeanFilter(

			DipImage<uchar> imfgMask(w,h);imfgMask.init(0);
			//~ imageMatting_2(frames[p],frames[p-1],imMask,imfgMask);
			imageMatting_4(frames[p],frames[p-1],imMask,imfgMask);
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
			//~ ImCopy(imFGLarget,imMask);
			
			
			
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
			ImCopy(imBGFinal,imMask);

			
			DipColorImage<uchar> imFGColorBGR(w,h);imFGColorBGR.init(0);
				
			
			ImCopy(frames[p],imBGFinal,imFGColor);
			cvtColor(imFGColor,imFGColorBGR,CV_RGB2BGR);
			framesFG.push_back(imFGColorBGR);
			dbi.WriteImage(imFGColor, "imFGColor" + boost::lexical_cast<string>(p) + ".png");
			
			count = count +1;
			
		}
	}
		
	
	const std::string NAME = dir_path_out.string() +"/output.avi";
	VideoWriter outputVideo; 
	Size S = Size((int) videoin.get(CV_CAP_PROP_FRAME_WIDTH),(int) videoin.get(CV_CAP_PROP_FRAME_HEIGHT));   
	
	outputVideo.open(NAME, CV_FOURCC('D','I','V','X'), videoin.get(CV_CAP_PROP_FPS), S, true);

	if (!outputVideo.isOpened())
	{
		cout  << "Could not open the output video for write: "  << endl;
		return -1;
	}
	for(int p = 0;p < frameNumber-1; p++){ //Show the image captured in the window and repeat

		if (frames[p+1].empty()) break;   // check if at end			
		outputVideo << framesFG[p];
	}
	outputVideo.release();
	
	
	return 0;
	
}


