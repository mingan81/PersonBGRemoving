

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
			
			
			DipImage<uchar> imIn;
			readDipImage(im_path.string().c_str(),imIn);
					
			int w  = imIn.width();
			int h  = imIn.height();;
			
			DipImage<uchar> imInInvert(w,h); imInInvert.init(0);
			DipImage<double> imInInvertD(w,h); imInInvertD.init(0);
			ImInvert(imIn,imInInvert);
			ImCopy(imInInvert,imInInvertD);
			
			DipImage<int> imLabel(w,h); imLabel.init(0);
			ImLabelFlatZones_WithCriterion(imInInvert,imInInvertD,Area,imLabel);


			DipImage<uchar> imBG(w,h); imBG.init(0);
			DipImage<uchar> imBGFinal(w,h); imBGFinal.init(0);
			ImThreshold(imLabel, int(0),int(w*h/20), uchar(255),uchar(0),imBG);
			
			ImAddImage(imIn,Clip,imBG,imBGFinal);
			dbi.WriteImage(imBGFinal,"imBGFinal.png");			
			
			
		}
		
	}


	return 0;
}

