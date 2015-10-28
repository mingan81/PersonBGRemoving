

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
			DipImage<int> ImEdgeAcc(w,h);ImEdgeAcc.init(0);
			
			DipImage<uchar> ImEdgesPrevious(w,h);ImEdgesPrevious.init(0);
			DipImage<uchar> ImEdgesCurrent(w,h);ImEdgesCurrent.init(0);
			
			int frameNumber = frames.size();
			for(int i = 0;i < frameNumber; i++){
				if(!frames[i].empty()){
					ImEdgesCurrent.init(0);
					DipImage<uchar> ImGreen(w,h);ImGreen.init(0);
					DipImage<uchar> ImRed(w,h);ImRed.init(0);
					DipImage<uchar> ImBlue(w,h);ImBlue.init(0);
					DipImage<uchar> ImGrayscale(w,h);ImGrayscale.init(0);					
					ImColorBandSeparation(frames[i],ImRed,ImGreen,ImBlue);
					cvtColor(frames[i],ImGrayscale,CV_BGR2GRAY);
					blur( ImGrayscale, ImGrayscale, Size(3,3) );
					DipImage<uchar> edges(w,h);edges.init(0);
					DipImage<int> ImEdgesAdd(w,h);ImEdgesAdd.init(0);
					DipImage<uchar> _img(w,h);_img.init(0);
					double otsu_thresh_val = cv::threshold(ImGrayscale, _img, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
					double high_thresh_val  = otsu_thresh_val*0.5,lower_thresh_val = otsu_thresh_val * 0.1;
					Canny(ImGrayscale, edges,lower_thresh_val,high_thresh_val);
					dbi.WriteImage(edges, boost::lexical_cast<string>(i) +"edges.png");
					ImCopy(edges,ImEdgesCurrent);
					ImAddImage(ImEdgesCurrent,NoClip,ImEdgesPrevious,ImEdgesAdd);
					DipImage<int> ImEdgesDiff(w,h);ImEdgesDiff.init(0);

					ImDivConst(ImEdgesAdd,int(2),ImEdgesDiff);
					//~ dbi.WriteImage(ImEdgesDiff, boost::lexical_cast<string>(i) +"ImEdgesDiff.png");
					ImCopy(ImEdgesCurrent,ImEdgesPrevious);
					
					ImAddImage(ImEdgeAcc,NoClip,edges,ImEdgeAcc);
				}
			}
			DipImage<uchar> ImEdgebased(w,h);ImEdgebased.init(0);
			dbi.WriteImage(ImEdgeAcc, "ImEdgeAcc.png");
			
			ImThreshold(ImEdgeAcc,int(0.1*frameNumber*255), int(frameNumber*255*2/3),uchar(255),uchar(0),ImEdgebased);
			dbi.WriteImage(ImEdgebased, "ImEdgebased.png");
			
			
							

			
			
		}
		
	}


	return 0;
}

