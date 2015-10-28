

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
			
			std::string genstring = intermediate_path.string() + "generalData.txt";
			fstream generaldatafilestr;
			generaldatafilestr.open(genstring.c_str(), ios_base::out);
			
			std::vector<DipColorImage<uchar> > frames;
			DipVideo videoin(im_path.string().c_str());

			videoin.decodeFrames(frames);
			
			int w  = videoin.getFrameWidth();
			int h  = videoin.getFrameHeight();
			
			
			std::vector<DipColorImage<uchar> > framesHSV;
			int frameNumber = frames.size();
			for(int i = 0;i < frameNumber; i++){
				if(!frames[i].empty()){
					DipColorImage<uchar> ImTemp(w,h);ImTemp.init(0);
					cvtColor(frames[i],ImTemp,CV_RGB2HSV);
					framesHSV.push_back(ImTemp);
					//~ dbi.WriteImage(frames[i], boost::lexical_cast<string>(i) +"frames.png");
					//~ DipImage<uchar> ImGreen(w,h);ImGreen.init(0);
					//~ DipImage<uchar> ImRed(w,h);ImRed.init(0);
					//~ DipImage<uchar> ImBlue(w,h);ImBlue.init(0);
					//~ ImColorBandSeparation(frames[i],ImRed,ImGreen,ImBlue);
					//~ DipImage<uchar> edges(w,h);edges.init(0);
					//~ Canny(ImGreen, edges,1,3);
					//~ dbi.WriteImage(edges, boost::lexical_cast<string>(i) +"edges.png");
				}
			}
			
			//~ for(int i = 0; i < w; ++i){
				//~ for(int j=0;j < h; ++j){
					//~ generaldatafilestr << "(" << i << "_" << j << ")" <<",";
				//~ }
			//~ }
			//~ generaldatafilestr <<endl;
			//~ for(int p = 0;p < frameNumber; p++){
				//~ for(int i = 0; i < w; ++i){
					//~ for(int j=0;j < h; ++j){
						//~ generaldatafilestr << int(framesHSV[p](i,j)[2]) <<",";
					//~ }
				//~ }
				//~ generaldatafilestr << endl;
			//~ }
			//~ 
			//~ generaldatafilestr.flush();
			//~ generaldatafilestr.close();
			//dbi.WriteImage(frames[50], "frames.png");

			
			DipImage<uchar> BGestimation(w,h); BGestimation.init(0);
			DipImage<double> imStddev(w,h); imStddev.init(0);
			detectBGEstimation(frames,BGestimation,imStddev);
			dbi.WriteImage(BGestimation, "BGestimation.png");
			dbi.WriteImage(imStddev, "imStddev.png");
			
			
			DipImage<uchar> imStddevUchar(w,h); imStddevUchar.init(0);
			for(int i = 0; i < w; ++i){
				for(int j=0;j < h; ++j){
					imStddevUchar(i,j) = uchar(imStddev(i,j)/3);
				}
			}
			
			
			
			//Ostu method(didnot work)
			//~ DipImage<unsigned char> img_bw(w,h); img_bw.init(0);
			//~ cv::threshold(imStddevUchar, img_bw, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
			//~ dbi.WriteImage(img_bw, "img_bw.png");
			//custom thresholding(didnot work)
			//~ DipImage<uchar> imMask(w,h);imMask.init(255);
			//~ DipImage<uchar> imOut(w,h);imOut.init(255);
			//~ ImAutoThresh(imStddev,imMask,imOut);
			//~ dbi.WriteImage(imOut, "imOut.png");
			
			
			DipImage<uchar> imForeGround(w,h);imForeGround.init(0);
			DipImage<uchar> imMask(w,h);imMask.init(255);
			std::vector<double> meanStddev = MeasStatistics(imStddev,imMask,all);
			dbi.WriteInfo("%f\n",meanStddev[0]);
			dbi.WriteInfo("%f\n",meanStddev[1]);
			dbi.WriteInfo("%f\n",meanStddev[2]);
			std::pair<double,double> minMax = MeasMinMax(imStddev,imMask);
			dbi.WriteInfo("%f\n",minMax.first);
			dbi.WriteInfo("%f\n",minMax.second);
			
			double factor = 0.1 * meanStddev[0]/3;
			double lowBorder = meanStddev[0];
			bool isGeneralBG = false;
			

			DipImage<uchar> imFGOpen(w,h);imFGOpen.init(0);
			DipImage<uchar> imFGClose(w,h);imFGClose.init(0);
			DipImage<int> imLabel(w,h);imLabel.init(0);
			DipImage<uchar> imBG(w,h);imBG.init(0);

			while(!isGeneralBG){
				ImThreshold(imStddev, lowBorder,minMax.second, uchar(255),uchar(0),imForeGround);
				param.clear();
				param.push_back(int(w/100));
				ImOpen(imForeGround,Fast,param,imFGOpen);
				param.clear();
				param.push_back(int(w/50));
				ImClose(imFGOpen,Fast,param,imFGClose);
				DipImage<double> imFGCloseDouble(w,h);imFGCloseDouble.init(0);
				ImCopy(imFGClose,imFGCloseDouble);
				ImLabelFlatZones_WithCriterion(imFGClose,imFGCloseDouble,Area,imLabel);
				std::pair<double,double> minmaxImLabel = MeasMinMax(imLabel,imMask);
				
				ImThreshold(imLabel, int(minmaxImLabel.second),int(minmaxImLabel.second), uchar(255),uchar(0),imBG);
				if( double(ImMeasVolume(imBG)/255) < 0.4 * w * h){
					lowBorder = lowBorder - factor;
				}
				else {isGeneralBG = true;}

			
			}
			dbi.WriteImage(imFGOpen, "imFGOpen.png");
			dbi.WriteImage(imFGClose, "imFGClose.png");
			dbi.WriteImage(imForeGround, "imForeGround.png");
			dbi.WriteImage(imLabel, "imLabel.png");
			dbi.WriteImage(imBG, "imBG.png");
			
			//~ std::vector<DipColorImage<uchar> > framesOut; framesOut.resize(frameNumber);
			//~ for(int i = 0;i < frameNumber; i++){
				//~ if(!frames[i].empty())
					//~ getForeGround(frames[i],imBG,framesOut[i]);
					//~ dbi.WriteImage(framesOut[i], boost::lexical_cast<string>(i) +"frames.png");
			//~ }
			
			
			
			
		}
		
	}


	return 0;
}

