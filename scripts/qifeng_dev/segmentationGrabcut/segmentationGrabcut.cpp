

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
	if(argc!=4){
			std::cout << "Usage: " << argv[0] << " infile infile(general Mask)  outfile"<<std::endl;
	}
	
	path dir_path(argv[1]);
	path dir_path_out(argv[3]);
	path dir_path_mask(argv[2]);
	
	path im_name;	
	path im_path;
	path mask_path;
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
			
			mask_path = dir_path_mask/im_name.stem();
			std::cout << mask_path.string().c_str() << std::endl;
			intermediate_path = dir_path_out/im_name.stem()/"/";
			if (!exists(intermediate_path))
				create_directories(intermediate_path);
			dbi.SetPath(intermediate_path);
					
			std::vector<DipColorImage<uchar> > frames;
			DipVideo videoin(im_path.string().c_str());
			string pathmask = mask_path.string() + ".png";
			DipImage<uchar> imGMask;
			readDipImage(pathmask.c_str(),imGMask);
			dbi.WriteImage(imGMask,"imGMask.png");

			videoin.decodeFrames(frames);
			
			int w  = videoin.getFrameWidth();
			int h  = videoin.getFrameHeight();
			
			
			DipImage<uchar> imGrabMask(w,h);imGrabMask.init(0);
			for(int i = 0; i< w;++i){
				for(int j = 0; j< h; ++j){
					if(imGMask(i,j) == 0) {imGrabMask(i,j) = cv::GC_BGD;}
					else{
							if( i> w/2 - 10 && i < w/2 + 10 && j > h/2) imGrabMask(i,j) = cv::GC_FGD;
							else imGrabMask(i,j) = cv::GC_PR_FGD;
					}
				}
			}

			int frameNumber = frames.size();
			for(int p = 0;p < frameNumber; p++){
				if(!frames[p].empty()){
					cv::Mat bgModel,fgModel;
					//~ DipImage<double> fgModel(w,h);
					DipImage<uchar> imGrabMaskTmp(w,h);
					cv::Rect rectangle(1,1,h-1,w-1);
					ImCopy(imGrabMask,imGrabMaskTmp);
					//~ dbi.WriteImage(imGrabMaskTmp, boost::lexical_cast<string>(p) +"imGrabMaskTmp.png");
					grabCut(frames[p], imGrabMaskTmp,rectangle,bgModel,fgModel,7, cv::GC_INIT_WITH_MASK); 
					
					DipImage<uchar> imFG(w,h);imFG.init(0);
					for(int i = 0; i< w;++i){
						for(int j = 0; j< h; ++j){
							if(imGrabMaskTmp(i,j) == GC_FGD || imGrabMaskTmp(i,j) == 3 ) {imFG(i,j) = 255;}
					
						}
					}
					//~ dbi.WriteImage(imGrabMaskTmp, boost::lexical_cast<string>(p) +"imGrabMaskTmpf.png");
					dbi.WriteImage(imFG, boost::lexical_cast<string>(p) +"imFG.png");
					
					
				}
			}
			


			
			
			
			
			
		}
		
	}


	return 0;
}

