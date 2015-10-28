

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream> // for standard I/O
#include <string>   // for strings
#include <iomanip>  // for controlling float print precision
#include <sstream>  // string to number conversion
#include <boost/lexical_cast.hpp>



#include "../../../include/Include.hpp"



// Function main
int main(int argc, char * argv[])
{
	if(argc!=3){
			std::cout << "Usage: " << argv[0] << " infile  outfile"<<std::endl;
	}
	
	path dir_path(argv[1]);
	path dir_path_out(argv[2]);
	
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
					
			std::vector<DipColorImage<uchar> > frames;
			DipVideo videoin(im_path.string().c_str());

			videoin.decodeFrames(frames);
			
			int w  = videoin.getFrameWidth();
			int h  = videoin.getFrameHeight();
			
			
			std::vector<DipImage<uchar> > framesV;
			int frameNumber = frames.size();
			for(int i = 1;i < frameNumber; i++){
				if(!frames[i].empty()){
					DipColorImage<uchar> ImTempPrevious(w,h);ImTempPrevious.init(0);
					DipColorImage<uchar> ImTempCurrent(w,h);ImTempCurrent.init(0);
					cvtColor(frames[i-1],ImTempPrevious,CV_RGB2HSV);
					cvtColor(frames[i],ImTempCurrent,CV_RGB2HSV);
					DipImage<uchar> ImTempH(w,h);ImTempH.init(0);
					DipImage<uchar> ImTempS(w,h);ImTempS.init(0);
					DipImage<uchar> ImTempVPrevious(w,h);ImTempVPrevious.init(0);
					DipImage<uchar> ImTempVCurrent(w,h);ImTempVCurrent.init(0);
					DipImage<uchar> ImMotion(w,h);ImMotion.init(0);
					ImColorBandSeparation(ImTempPrevious,ImTempH,ImTempS,ImTempVPrevious);
					ImColorBandSeparation(ImTempCurrent,ImTempH,ImTempS,ImTempVCurrent);
					for(int n = 0; n < w; ++n){
							for(int m=0;m < h;++m){
									if(abs(ImTempVPrevious(n,m)-ImTempVCurrent(n,m))>5)ImMotion(n,m) = 255;
							}
					}
					framesV.push_back(ImMotion);
					dbi.WriteImage(ImMotion, boost::lexical_cast<string>(i) +"ImMotion.png");
				}
			}
			
			
			
			
			
			
			
			
			
		}
		
	}


	return 0;
}

