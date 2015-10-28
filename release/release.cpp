
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

inline bool exists (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}


// Function main
int main(int argc, char * argv[])
{
	if(argc!=4){
			std::cout << "Usage: " << argv[0] << " path filename pathmodel"<<std::endl;
	}
	
	std::string path(argv[1]);
	std::string filename(argv[2]);
	std::string path_model(argv[3]);
	std::vector<ParamT> param; param.clear();

	DipColorImage<uchar> imFistFrame;
	std::string path_image = path + filename + "0001.png";
	readDipImage(path_image.c_str(),imFistFrame);
	int w  = imFistFrame.width();
	int h  = imFistFrame.height();
	
	std::cout << "w =" << w<<endl;

	DipImage<uchar> imMask(w,h);imMask.init(0);
	DipImage<uchar> imFirstFrameMask(w,h);imFirstFrameMask.init(0);
	DipColorImage<uchar> imPrivous(w,h);imPrivous.init(0);
	DipColorImage<uchar> imCurrent(w,h);imCurrent.init(0);
	DipColorImage<uchar> imFirstFrameResult(w,h);imFirstFrameResult.init(0);
	
	maskDetection(imFistFrame,path_model,imMask);
	ImCopy(imFistFrame,imMask,imFirstFrameResult);
	ImCopy(imMask,imFirstFrameMask);
	
	std::string path_image_out;
	path_image_out = path + filename + "_out0001.png";
	writeDipImage(imFirstFrameResult,path_image_out.c_str());

	//~ std::vector<DipColorImage<uchar> > frames;frames.clear();
	ImCopy(imFistFrame,imCurrent);
	bool isEndImage = false;
	int countImage = 2;
	while(!isEndImage){
		imPrivous.init(0);
		ImCopy(imCurrent,imPrivous);
		
		if(countImage<10){
			path_image = path + filename + "000" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else if(countImage<100){
			path_image = path + filename + "00" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else if(countImage<1000){
			path_image = path + filename + "0" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else{
			path_image = path + filename + boost::lexical_cast<string>(countImage) + ".png";
		}
		imCurrent.init(0);
		readDipImage(path_image.c_str(),imCurrent);
		
		DipColorImage<uchar> imFGColor(w,h);imFGColor.init(0);					
			
		//@todo Add motion mask


		DipImage<uchar> imfgMask(w,h);imfgMask.init(0);
		imageMatting_2(imCurrent,imPrivous,imMask,imfgMask);
		//~ imageMatting_3(imCurrent,imPrivous,imFistFrame,imMask,imFirstFrameMask,imfgMask);

		//~ imageMatting_4(imCurrent,imPrivous,imMask,imfgMask);
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
		//~ ImCopy(imBGFinal,imMask);

		
		//~ DipColorImage<uchar> imFGColorBGR(w,h);imFGColorBGR.init(0);
			
		
		ImCopy(imCurrent,imBGFinal,imFGColor);
		//~ cvtColor(imFGColor,imFGColorBGR,CV_RGB2BGR);
		//~ framesFG.push_back(imFGColorBGR);
		//~ framesFG.push_back(imFGColor);
		//~ dbi.WriteImage(imFGColor, "imFGColor" + boost::lexical_cast<string>(p) + ".png");
		

		if(countImage<10){
			path_image_out = path + filename + "_out000" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else if(countImage<100){
			path_image_out = path + filename + "_out00" + boost::lexical_cast<string>(countImage) + ".png";		
		}
		else if(countImage<1000){
			path_image_out = path + filename + "_out0" + boost::lexical_cast<string>(countImage) + ".png";		
		}
		else{
			path_image_out = path + filename + "_out" + boost::lexical_cast<string>(countImage) + ".png";
		}		
		
		
		writeDipImage(imFGColor,path_image_out.c_str());
		
		countImage = countImage + 1;
		if(countImage<10){
			path_image = path + filename + "000" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else if(countImage<100){
			path_image = path + filename + "00" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else if(countImage<1000){
			path_image = path + filename + "0" + boost::lexical_cast<string>(countImage) + ".png";
		}
		else{
			path_image = path + filename + boost::lexical_cast<string>(countImage) + ".png";			
		}
		//~ std::cout << countImage << endl;
				
		if(!exists(path_image)){
			std::cout << path_image << endl;
			isEndImage = true;
			countImage = countImage -1;
		}
	}
	
	


			
	

	
	
	return 0;
	
}


