

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
			result_path = dir_path_out/im_name.stem();
			std::cout << im_path.string().c_str() << std::endl;
			std::vector<ParamT> param; param.clear();
			
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
			
			DipImage<double> imStddev(w,h); imStddev.init(0);
			{
				std::vector<DipImage<uchar> > framesV;
				
				for(int i = 0;i < frameNumber; i++){
					if(!frames[i].empty()){
						DipColorImage<uchar> ImTemp(w,h);ImTemp.init(0);
						cvtColor(frames[i],ImTemp,CV_RGB2HSV);
						DipImage<uchar> ImTempH(w,h);ImTempH.init(0);
						DipImage<uchar> ImTempS(w,h);ImTempS.init(0);
						DipImage<uchar> ImTempV(w,h);ImTempV.init(0);
						ImColorBandSeparation(ImTemp,ImTempH,ImTempS,ImTempV);
						framesV.push_back(ImTempV);
					}
				}
				


				
				DipImage<uchar> BGestimation(w,h); BGestimation.init(0);
				
				detectBGEstimation(framesV,imStddev);			
				dbi.WriteImage(imStddev, "imStddev.png");
			}
			
			DipImage<double> imStddevUchar(w,h); imStddevUchar.init(0);
			{
				DipImage<uchar> imMask(w,h);imMask.init(255);
				std::pair<double,double> minMax = MeasMinMax(imStddev,imMask);
				dbi.WriteInfo("%f\n",minMax.first);
				dbi.WriteInfo("%f\n",minMax.second);
				ImDivConst(imStddev,double(255/minMax.second),imStddevUchar);
				dbi.WriteImage(imStddevUchar, "imStddevUchar.png");
			}
			
			
			DipImage<uchar> imBG(w,h);imBG.init(0);
			{

			
				DipImage<uchar> imForeGround(w,h);imForeGround.init(0);
				DipImage<uchar> imMask(w,h);imMask.init(255);
				std::vector<double> meanStddev = MeasStatistics(imStddevUchar,imMask,all);
				dbi.WriteInfo("%f\n",meanStddev[0]);
				dbi.WriteInfo("%f\n",meanStddev[1]);
				dbi.WriteInfo("%f\n",meanStddev[2]);
				std::pair<double,double> minMax = MeasMinMax(imStddevUchar,imMask);
				dbi.WriteInfo("%f\n",minMax.first);
				dbi.WriteInfo("%f\n",minMax.second);
				
				double factor = 0.01 * meanStddev[0]/3;
				double lowBorder = meanStddev[0];
				bool isGeneralBG = false;
				

				DipImage<uchar> imFGOpen(w,h);imFGOpen.init(0);
				DipImage<uchar> imFGClose(w,h);imFGClose.init(0);
				DipImage<int> imLabel(w,h);imLabel.init(0);
				

				while(!isGeneralBG){
					ImThreshold(imStddevUchar, lowBorder,minMax.second, uchar(255),uchar(0),imForeGround);
					param.clear();
					param.push_back(int(w/100));
					ImOpen(imForeGround,Fast,param,imFGOpen);
					param.clear();
					param.push_back(int(w/50));
					ImClose(imFGOpen,Fast,param,imFGClose);
					DipImage<double> imFGCloseDouble(w,h);imFGCloseDouble.init(0);
					ImCopy(imFGClose,imFGCloseDouble);
					ImLabelFlatZones_WithCriterion(imFGClose,imFGCloseDouble,Area,imLabel);
					//~ DipImage<double> imFGOpenDouble(w,h);imFGOpenDouble.init(0);
					//~ ImCopy(imFGOpen,imFGOpenDouble);
					//~ ImLabelFlatZones_WithCriterion(imFGOpen,imFGOpenDouble,Area,imLabel);
					
					
					
					std::pair<double,double> minmaxImLabel = MeasMinMax(imLabel,imMask);
					
					ImThreshold(imLabel, int(minmaxImLabel.second),int(minmaxImLabel.second), uchar(255),uchar(0),imBG);
					if( double(ImMeasVolume(imBG)/255) < 0.3 * w * h){
						lowBorder = lowBorder - factor;
					}
					else {isGeneralBG = true;}

				
				}
				dbi.WriteImage(imFGOpen, "imFGOpen.png");
				dbi.WriteImage(imFGClose, "imFGClose.png");
				dbi.WriteImage(imForeGround, "imForeGround.png");
				dbi.WriteImage(imLabel, "imLabel.png");
				dbi.WriteImage(imBG, "imBG.png");
			}
			
			DipImage<uchar> imBGFinal(w,h); imBGFinal.init(0);
			{
				DipImage<uchar> imInInvert(w,h); imInInvert.init(0);
				DipImage<double> imInInvertD(w,h); imInInvertD.init(0);
				ImInvert(imBG,imInInvert);
				ImCopy(imInInvert,imInInvertD);
				
				DipImage<int> imLabel2(w,h); imLabel2.init(0);
				ImLabelFlatZones_WithCriterion(imInInvert,imInInvertD,Area,imLabel2);


				DipImage<uchar> imBG2(w,h); imBG2.init(0);
				
				ImThreshold(imLabel2, int(0),int(w*h/20), uchar(255),uchar(0),imBG2);
				
				ImAddImage(imBG,Clip,imBG2,imBGFinal);
				dbi.WriteImage(imBGFinal,"imBGFinal.png");
			}	
			
			DipImage<uchar> imBGFinalDilate(w,h); imBGFinalDilate.init(0);
			param.clear();
			param.push_back(int(w/25));
			ImDilate(imBGFinal,Fast,param,imBGFinalDilate);
			dbi.WriteImage(imBGFinalDilate,"imBGFinalDilate.png");
			
			
			DipImage<uchar> imGrabMask(w,h);imGrabMask.init(0);
			for(int i = 0; i< w;++i){
				for(int j = 0; j< h; ++j){
					if(imBGFinalDilate(i,j) == 0) {imGrabMask(i,j) = cv::GC_BGD;}
					else{
							if( i> w/2 - 10 && i < w/2 + 10 && j > h/2) imGrabMask(i,j) = cv::GC_FGD;
							else imGrabMask(i,j) = cv::GC_PR_FGD;
					}
				}
			}
			
			std::vector<DipColorImage<uchar> > framesFG;
			for(int p = 0;p < frameNumber; p++){
				if(!frames[p].empty()){
					if(p==0){
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
						dbi.WriteImage(imGrabMaskTmp, boost::lexical_cast<string>(p) +"imGrabMaskTmpf.png");
						dbi.WriteImage(imFG, boost::lexical_cast<string>(p) +"imFG.png");
						DipColorImage<uchar> imFGColor(w,h);imFGColor.init(0);
						DipColorImage<uchar> imFGColorBGR(w,h);imFGColorBGR.init(0);
						
						//EXTRACT THE LARGEST ONE
						DipImage<double> imFGDouble(w,h);imFGDouble.init(0);
						DipImage<int> imFGLabel(w,h);imFGLabel.init(0);
						DipImage<uchar> imMask(w,h);imMask.init(255);
						DipImage<uchar> imFGLarge(w,h);imFGLarge.init(255);
						DipImage<uchar> imFGFinal(w,h);imFGFinal.init(255);
						ImCopy(imFG,imFGDouble);
						ImLabelFlatZones_WithCriterion(imFG,imFGDouble,Area,imFGLabel);
						std::pair<double,double> minmaxImFGLabel = MeasMinMax(imFGLabel,imMask);					
						ImThreshold(imFGLabel, int(minmaxImFGLabel.second),int(minmaxImFGLabel.second), uchar(255),uchar(0),imFGLarge);
						
						
						
						param.clear();
						param.push_back(int(10));
						ImDilate(imFGLarge,Fast,param,imFGFinal);
						dbi.WriteImage(imFGFinal, boost::lexical_cast<string>(p) +"imFGFinal.png");

						
						DipImage<uchar> ImGrayscale(w,h);ImGrayscale.init(0);					
						DipImage<uchar> ImGrayscaleBlur(w,h);ImGrayscaleBlur.init(0);					
						DipImage<uchar> ImEdgeinMask(w,h);ImGrayscale.init(0);					
						cvtColor(frames[p],ImGrayscale,CV_BGR2GRAY);
						DipImage<uchar> edges(w,h);edges.init(0);
						DipImage<uchar> imEdge(w,h);imEdge.init(0);
						DipImage<uchar> _img(w,h);_img.init(0);
						blur( ImGrayscale, ImGrayscaleBlur, Size(3,3) );
						double otsu_thresh_val = cv::threshold(ImGrayscaleBlur, _img, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
						double high_thresh_val  = otsu_thresh_val* 0.5,lower_thresh_val = otsu_thresh_val * 0.1;
						Canny(ImGrayscaleBlur, edges,lower_thresh_val,high_thresh_val);
						
						dbi.WriteImage(edges,"edges.png");
						ImCopy(edges,imFGFinal,ImEdgeinMask);
						dbi.WriteImage(ImEdgeinMask,"ImEdgeinMask.png");
						
						ImThreshold(ImEdgeinMask,uchar(255), uchar(255),uchar(255),uchar(0),imEdge);
						dbi.WriteImage(imEdge,"imEdge.png");
						
						vector<vector<Point> > contours;
						vector<Vec4i> hierarchy;

						/// Detect edges using canny 
						/// Find contours
						findContours( imEdge, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_TC89_KCOS, Point(0, 0) );
						
						
						DipImage<uchar> drawing(w,h);drawing.init(0);					

						for( int i = 0; i< contours.size(); i++ )
						{
							Scalar color = Scalar(255,255,255);
							drawContours( drawing, contours, i, color, CV_FILLED, 8);
						}
						
						dbi.WriteImage(drawing,"drawing.png");
						
									
						
						
						ImCopy(frames[p],drawing,imFGColor);
						cvtColor(imFGColor,imFGColorBGR,CV_RGB2BGR);
						framesFG.push_back(imFGColorBGR);
						dbi.WriteImage(imFGColor, boost::lexical_cast<string>(p) +"imFGColor.png");
					}
					else{
					}
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
			for(int p = 0;p < frameNumber; p++){ //Show the image captured in the window and repeat
       
				if (frames[p].empty()) break;   // check if at end			
				outputVideo << framesFG[p];
			}
			outputVideo.release();
			
		}
		
	}


	return 0;
}
