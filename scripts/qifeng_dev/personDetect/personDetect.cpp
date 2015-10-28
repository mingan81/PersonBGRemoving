

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream> // for standard I/O
#include <string>   // for strings
#include <iomanip>  // for controlling float print precision
#include <sstream>  // string to number conversion
#include <boost/lexical_cast.hpp>



#include "../../../include/Include.hpp"

using namespace boost::filesystem;



// Function main
int main(int argc, char * argv[])
{
	if(argc!=4){
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
			std::cout << im_path.string().c_str() << std::endl;
			std::vector<ParamT> param; param.clear();
			
			intermediate_path = dir_path_out/im_name.stem()/"/";
			if (!exists(intermediate_path))
				create_directories(intermediate_path);
			dbi.SetPath(intermediate_path);
			
			std::string path_model = dir_path_model.string();
			
			DipColorImage<uchar> frame;
			readDipImage(im_path.string().c_str(),frame);
			
			int w  = frame.width();
			int h  = frame.height();
			DipImage<uchar> imMask(w,h);imMask.init(0);
			maskDetection(frame,path_model,imMask);
			
			dbi.WriteImage(imMask, "imMaskBeforeMatting.png");
			DipColorImage<uchar> imFGColor(w,h);imFGColor.init(0);
			//~ ImCopy(frame,imMask,imFGColor);
			//~ dbi.WriteImage(imFGColor, "imFGColor.png");

			DipImage<uchar> imfgMask(w,h);imfgMask.init(0);
			imageMatting_initial(frame,imMask,imfgMask);
			dbi.WriteImage(imfgMask, "imfgMask.png");
			
						
			//extract the largest one
			DipImage<double> imFGD(w,h); imFGD.init(0);
			DipImage<uchar> imFGLarget(w,h); imFGLarget.init(0);
			DipImage<int> imFGLabel(w,h); imFGLabel.init(0);
			ImLabelFlatZones_WithCriterion(imfgMask,imFGD,Area,imFGLabel);
			std::pair<double,double> minmax = MeasMinMax(imFGLabel,imfgMask);
			ImThreshold(imFGLabel, int(minmax.second), int(minmax.second), uchar(255),uchar(0),imFGLarget);
			dbi.WriteImage(imFGLarget, "imFGLarget.png");
														
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
			
			
			ImCopy(frame,imBGFinal,imFGColor);
			dbi.WriteImage(imFGColor, "imFGColor.png");
										
		}
		
	}


	return 0;
}

