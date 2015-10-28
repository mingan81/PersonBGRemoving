

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream> // for standard I/O
#include <string>   // for strings
#include <iomanip>  // for controlling float print precision
#include <sstream>  // string to number conversion
#include <boost/lexical_cast.hpp>

#include <opencv2/video/background_segm.hpp>


#include "../../../include/Include.hpp"
const int BS_HISTORY = 10;
const float BS_THRESHOLD = 16.0;
const bool BS_SHADOW = true;
const double LEARNING_RATE = -1.0;
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

template <typename T> void removeBackground( const DipColorImage<T> & frame, DipImage<T> & imMaskOut,DipColorImage<T> & imMaskBG, cv::Ptr<BackgroundSubtractor> bgSub )
{
	// Remove background and get foreground mask
	bgSub->operator()(frame, imMaskOut, LEARNING_RATE);
	bgSub->getBackgroundImage(imMaskBG);
	//~ // Apply opening to reduce noise
	//~ int morphoType = MORPH_ELLIPSE;
	//~ Mat element = getStructuringElement( morphoType, Size( 2*MORPHO_SIZE + 1, 2*MORPHO_SIZE + 1 ), Point( MORPHO_SIZE, MORPHO_SIZE ) );
	//~ morphologyEx( fgMask, fgMask, 2, element ); // 2 indicates opening operation
}

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
			
			
			Ptr<BackgroundSubtractor> bgSub = createBgSubtractor();
			
			std::vector<DipImage<uchar> > framesMotion;
			int frameNumber = frames.size();
			for(int i = 0;i < frameNumber; i++){
				if(!frames[i].empty()){
					DipImage<uchar> ImTemp(w,h);ImTemp.init(0);
					DipColorImage<uchar> ImTempBG(w,h);ImTempBG.init(0);
					// -- Get foreground/motion mask (head image)
					removeBackground(frames[i],ImTemp,ImTempBG, bgSub );				
					framesMotion.push_back(ImTemp);
					dbi.WriteImage(ImTempBG, boost::lexical_cast<string>(i) +"ImTempBG.png");
				}
			}
			
		}
		
	}


	return 0;
}

