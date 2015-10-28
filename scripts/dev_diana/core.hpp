 /** 
  * @file core.hpp
  * @author Qifeng Gan
  * @brief Core for Porjects  
  * In Core.hpp are enum, class, structs, etc... that are usefull for other hpp containts. 
  */

#ifndef CORE_HPP
#define CORE_HPP
 
#include <math.h>
#include <stdlib.h>
//#include <stddef.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <functional>
#include <queue>
#include <map>
#include <limits>
#include <list>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdarg.h>
#include <ctime>
#include <time.h>
//#include <assert.h>
//#include <typeinfo> 
 
#include "boost/filesystem.hpp"
#include "boost/random.hpp"
#include "boost/tuple/tuple_io.hpp"
#include <boost/thread/thread.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/array.hpp> 


#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>


#ifdef _TRACES_
#define DEBUG 1
#else
#define DEBUG 0
#endif

#ifdef DEBUG
namespace bf = boost::filesystem;
#endif
using namespace cv;

/*!
* DipImage inherits from opencv Mat_
*/
template <typename T> class DipImage : public cv::Mat_<T>{
	public:

		/*!
		* Default constructor
		*/
		DipImage() :cv::Mat_<T>(){};
		DipImage(const int w, const int h) :cv::Mat_<T>(h, w){};

		/*!
		* Accessors
		*/
		int width() const { return Mat_<T>::cols; };
		int height() const { return Mat_<T>::rows; };
		T const& operator() (int x, int y) const { return Mat_<T>::operator()(y, x); };
		T& operator()(int x, int y){ return Mat_<T>::operator()(y, x); }

		/*!
		* Scalar initialisation
		*/
		template<typename Tval>
		void init(const Tval val){ (*this).setTo(cv::Scalar(T(val))); }

};

/*!
* DipColor inherits from opencv Vec
*/
template <typename T> class DipColor : public cv::Vec<T, 3>{
	public:

		/**
		* Default constructor
		*/
		DipColor() :cv::Vec<T, 3>(){};
		DipColor(T r, T g, T b) :cv::Vec<T, 3>(r, g, b){};
		DipColor(T c) :cv::Vec<T, 3>(c, c, c){};

		T red() const { return((*this)[0]); }
		T green() const { return((*this)[1]); }
		T blue() const { return((*this)[2]); }

};

/*!
* DipColorImage inherits from opencv Mat_
*/
template <typename T> class DipColorImage : public cv::Mat_<Vec<T, 3> >{
	public:

		/**
		* Default constructor
		*/
		DipColorImage() :cv::Mat_<Vec<T, 3> >(){};
		DipColorImage(const int w, const int h) :cv::Mat_<Vec<T, 3> >(h, w){};

		/*!
		* Accessors
		*/
		int width() const { return Mat_<Vec<T, 3> >::cols; };
		int height() const { return Mat_<Vec<T, 3> >::rows; };
		Vec<T, 3>& operator()(int x, int y){ return Mat_<Vec<T, 3> >::operator()(y, x); }
		Vec<T, 3> const & operator()(int x, int y) const { return Mat_<Vec<T, 3> >::operator()(y, x); }


		/*!
		* Scalar initialisation
		*/
		void init(const DipColor<T> val){ Mat_<Vec<T, 3> > M(Mat_<Vec<T, 3> >::rows, Mat_<Vec<T, 3> >::cols, val); M.copyTo((*this)); }
		void init(const T val){ DipColor<T> val3(val); init(val3); }
};

/*!
 * video class inherits from VideoCapture (opencv)
 */
 
class video: public cv::VideoCapture{
	/*!
	 * constructor
	*/
	 video():cv::VideoCapture(){};
	 //~ video(const string& filename): cv::VideoCapture(const string& filename){};
	 video(const string& filename);
	 //~ video(int device): cv::VideoCapture(int device){};
	 video(int device);
	 
	 bool decode2frames(const string& filename);
	 bool decode2frames(int device);
	 bool isOpen();
	 void release();
	 void loadFrames(std::vector<DipColorImage<uchar> > &  images);
	 
	 private:
		std::vector<DipColorImage<uchar> > frames;
	 
	 
	
}


/*!
 * DipVideo class inherits from VideoCapture (OpenCV)
 */
class DipVideo : public cv::VideoCapture
{
public:
	DipVideo();
	DipVideo(const std::string&);
	DipVideo(int device);
	~DipVideo();

	std::string getFilename();
	int getCamera();
	double getFrameWidth();
	double getFrameHeight();
	double getFPS();
	double getNumFrames();

	void checkInit();
	void restartVideo();
	std::vector<DipColorImage<uchar> > decodeFrames();

private:
	std::string filename;
	int camera;
};


/*!
* @brief Import DIP Color Image from given path
* @param impath : path to the image
* @param[out] I : output image
*/
template<typename T> void readDipImage(const char * impath, DipColorImage<T> & I){
	Mat_<Vec<T, 3> > J = imread(impath, 1);
	cvtColor(J, I, CV_BGR2RGB);
};

/*!
* @brief Import DIP Image from given path
* @param impath : path to the image
* @param[out] I : output image
*/
template<typename T> void readDipImage(const char * impath, DipImage<T> & I)
{
	Mat_<T> J = imread(impath, 0);
	J.copyTo(I);
};

/*!
* @brief Export DIP Color Image to given path
* @param I : image to be written
* @param impath : path to the image
*/
template<typename T> void writeDipImage(const DipColorImage<T> & I, const char * impath)
{
	Mat_<Vec<T, 3> > J;
	imwrite(impath, I);
	cvtColor(I, J, CV_RGB2BGR);
	imwrite(impath, J);
};

/*!
* @brief Export DIP Image to given path
* @param I : image to be written
* @param impath : path to the image
*/
template<typename T> void writeDipImage(const DipImage<T> & I, const char * impath)
{
	imwrite(impath, I);
};
/*!
* @brief Export DIP Color Image to given path
* @param I : image to be written
* @param impath : path to the image
*/
template<typename T> void writeDipImage(const DipColorImage<T> & I, const char * impath, const char * compression)
{
	writeDipImage(I, impath);
};

/*!
* @brief Export DIP Image to given path
* @param I : image to be written
* @param impath : path to the image
*/
template<typename T> void writeDipImage(const DipImage<T> & I, const char * impath, const char * compression)
{
	writeDipImage(I, impath);
};


template <typename T> int sign(T val)
{
	if (val == T(0))
		return int(0);
	else if (T(0) < val)
		return int(1);
	else
		return int(-1);
}

 
/**
  * TimeInfo
  */
struct TimeInfo{
	public:
	/**
	 * constructor
	 */
	TimeInfo(){}
	/**
	  * Start
	 */
	void Start(){
		#if DEBUG
			tstart = clock_t(clock());
		#endif // DEBUG
	}

	/**
	  * Stop
	 */
	void Stop(){
		#if DEBUG
			tstop = clock();
			ti = (double(tstop) - double(tstart)) / CLOCKS_PER_SEC;
			printf("Processing time: %f seconds\n", ti);
		#endif // DEBUG
	}

	/**
	  * Stop
	  *
	  * @param level
	 */
	void Stop(const int level){
		#if DEBUG
			tstop = clock();
			ti = (double(tstop) - double(tstart)) / CLOCKS_PER_SEC;
			for (int i = 0; i< level; i++)
				printf("\t");
			printf("Processing time: %f seconds - level = %d\n", ti, level);
		#endif // DEBUG
	}

	private:
		clock_t tstart, tstop; 	/**< data members */
		double ti;				 	/**< data members */
};
 
 
/**
  * DebugInfo
 */
struct DebugInfo{
	public:
		/**
		 * default constructor
		 */
		DebugInfo(){
			level = 0;
			ttil.clear();
			EnterPath();
		}

		/**
		  * constructor
		  * @param _level
		 */
		DebugInfo(int _level){
			level = _level;
			ttil.clear();
			EnterPath();
		}

		int GetLevel(){return(level);} 				/**< GetLevel @return ? */
		void SetLevel(int _level){level = _level;} 	/**< SetLevel @param _level ? */

		/**
		  * GetPath
		  *
		  * @return
		 */
		std::string GetPath(){
			return(pathname);
		}

		/**
		  * SetPath
		  *
		  * @param _pathname
		 */
		void SetPath(std::string _pathname){
			pathname = _pathname;
		}

		/**
		  * SetPath
		  *
		  * @param _pathname
		 */
		void SetPath(bf::path _pathname){
			pathname = _pathname.string();
			// printf("Pathname = %s\n", pathname);
		}

		/**
		  * WriteInfo
		  *
		  * @param mask
		 */
		void WriteInfo(const char * mask, ...){
			#if DEBUG
				for (int i=0; i<level; i++)
					printf("\t");
				va_list vl;
				va_start(vl, mask);
				vprintf(mask, vl);
				va_end(vl);
			 #endif // DEBUG
		}

		/**
		  * WriteEnter
		  *
		  * @param mask
		*/
		void WriteEnter(const char * mask, ...){
			#if DEBUG
				if (level == 0)
					printf("\n");
				for (int i=0; i<level; i++)
					printf("\t");
				va_list vl;
				va_start(vl, mask);
				vprintf(mask, vl);
				va_end(vl);
				level++;

				TimeInfo ttis;
				ttis.Start();
				ttil.push_back(ttis);

			#endif // DEBUG
		}

		/**
		  * WriteOut
		  *
		  * @param mask
		 */
		void WriteOut(const char * mask, ...){
			#if DEBUG
				level--;

				ttil.back().Stop(level);
				ttil.pop_back();

				for (int i=0; i<level; i++)
					printf("\t");
				va_list vl;
				va_start(vl, mask);
				vprintf(mask, vl);
				va_end(vl);


			#endif // DEBUG
		}

		/**
		 * WriteImage
		 *
		 * @param imIn
		 * @param filename
		 * @param compression
		*/
		template<typename Tin> void WriteImage(const DipImage<Tin> & imIn, const std::string filename, const char * compression){
			#if DEBUG
				std::string temp(pathname);
				temp.append(filename);
				writeDipImage(imIn, temp.c_str(), compression);
			#endif
		}

		/**
		 * WriteImage
		 *
		 * @param imIn
		 * @param filename
		*/
		template<typename Tin> void WriteImage(const DipImage<Tin> & imIn, const std::string filename){
			#if DEBUG
				WriteImage(imIn, filename, "97");
			#endif
		}

		/**
		 * WriteImage
		 *
		 * @param imIn
		 * @param filename
		 * @param compression
		*/
		template<typename Tin> void WriteImage(const DipColorImage<Tin> & imIn, const std::string filename, const char * compression){
			#if DEBUG
				std::string temp(pathname);
				temp.append(filename);				
				writeDipImage(imIn, temp.c_str(), compression);
			#endif
		}

		/**
		 * WriteImage
		 *
		 * @param imIn
		 * @param filename
		*/
		template<typename Tin> void WriteImage(const DipColorImage<Tin> & imIn, const std::string filename){
			#if DEBUG
				WriteImage(imIn, filename, "97");
			#endif
		}

	private:
		int level;							/**< level of ??? */
		std::string pathname;				/**< ??? */

		clock_t tstart;					/**< ??? */
		clock_t tstop;						/**< ??? */
		double ti;							/**< ??? */
		std::vector<TimeInfo> ttil;		/**< ??? */

		/**
		  * EnterPath
		 */
		void EnterPath(){
			#if DEBUG
				printf("Please enter the path: ");

				getline(std::cin, pathname, '\n');
				printf("pathname = %s\n", pathname.c_str());
				printf("you are at\n");
				if (pathname == "")
					pathname = "/tmp";
				if (pathname == "\n")
					pathname = "/tmp";
				if (pathname.c_str()[pathname.length()-1] != '/')
					pathname.append("/");

				bf::path chemin(pathname);
				printf("Path = %s\n", pathname.c_str());
				if (!bf::exists( chemin ) )
					bf::create_directories( chemin );

			#endif
		}
};

extern DebugInfo dbi;


/*! @struct FZInfoNode
*
*/
template<typename T> struct FZInfoNode
{
	typedef T TypeTmp;

	int regionnumber;
	int x_min;
	int x_max;
	int y_min;
	int y_max;
	T val_min;
	T val_max;
	int area;
	double average;
	double stddev;
	long int volume;
	int nbptnoz;
	T dynamic;


	FZInfoNode()// Default constructor
	{
		x_max = y_max = area = regionnumber = 0;
		val_max = 0;
		volume = 0;
		average = 0;
		stddev = 0;
		dynamic = 0;
		nbptnoz = 0;
		x_min = y_min = std::numeric_limits<int>::max();
		val_min = std::numeric_limits<TypeTmp>::max();
	}

	FZInfoNode(int _regionnumber, int _x_min, int _x_max, int _y_min, int _y_max, T _val_min, T _val_max, int _area, double _average, double _stddev, long int _volume, T _dynamic, int _nbptnoz)
	{
		regionnumber = _regionnumber;
		x_min = _x_min; x_max = _x_max;
		y_min = _y_min; y_max = _y_max;
		val_min = _val_min; val_max = _val_max;
		area = _area;
		average = _average;
		stddev = _stddev;
		volume = _volume;
		dynamic = _dynamic;
		nbptnoz = _nbptnoz;
	}

	FZInfoNode<T>(const FZInfoNode<T> & fzin)
	{
		regionnumber = fzin.regionnumber;
		x_min = fzin.x_min;
		x_max = fzin.x_max;
		y_min = fzin.y_min;
		y_max = fzin.y_max;
		val_min = fzin.val_min;
		val_max = fzin.val_max;
		area = fzin.area;
		average = fzin.average;
		stddev = fzin.stddev;
		volume = fzin.volume;
		dynamic = fzin.dynamic;
		nbptnoz = fzin.nbptnoz;
	}

	void operator=(FZInfoNode fzin)
	{
		regionnumber = fzin.regionnumber;
		x_min = fzin.x_min;
		x_max = fzin.x_max;
		y_min = fzin.y_min;
		y_max = fzin.y_max;
		val_min = fzin.val_min;
		val_max = fzin.val_max;
		area = fzin.area;
		average = fzin.average;
		stddev = fzin.stddev;
		volume = fzin.volume;
		dynamic = fzin.dynamic;
		nbptnoz = fzin.nbptnoz;
	}

	void Reset()
	{
		x_max = y_max = area = regionnumber = nbptnoz = 0;
		val_max = 0;
		average = 0;
		stddev = 0;
		volume = 0;
		dynamic = 0;
		x_min = y_min = std::numeric_limits<int>::max();
		val_min = std::numeric_limits<TypeTmp>::max();
	}

	void DisplayNode()
	{
		std::cout << "Region number : " << regionnumber << std::endl;
		std::cout << "x_min : " << x_min << std::endl;
		std::cout << "x_max : " << x_max << std::endl;
		std::cout << "y_min : " << y_min << std::endl;
		std::cout << "y_max : " << y_max << std::endl;
		std::cout << "val_min : " << int(val_min) << std::endl;
		std::cout << "val_max : " << int(val_max) << std::endl;
		std::cout << "area : " << area << std::endl;
		std::cout << "average : " << double(average) << std::endl;
		std::cout << "stddev : " << double(stddev) << std::endl;
		std::cout << "volume : " << int(volume) << std::endl;
		std::cout << "dynamic : " << int(dynamic) << std::endl;
		std::cout << "nbptnoz : " << int(nbptnoz) << std::endl;
	}
};



/*! @struct Node
  * This structure is a container dedicated to connected components information gathering.
  */
 template<typename Tvalue, typename Tlabel> struct Node
 {
	 typedef Tvalue TypeValue;
	 typedef Tlabel TypeLabel;

	 public:		
		Node(){pos_x = pos_y = value = label = 0;}
		Node(int _x, int _y, TypeValue _value, TypeLabel _label){pos_x = _x; pos_y = _y; value = _value; label = _label;}
		void operator=(Node noeud){
			pos_x = noeud.pos_x;
			pos_y = noeud.pos_y;
			value = (TypeValue)noeud.value;
			label = (TypeLabel)noeud.label;
		}
		void SetPos(int _x, int _y) {pos_x = _x; pos_y = _y;};
		void SetXPos(int _x) {pos_x = _x;};
		void SetYPos(int _y) {pos_y = _y;};
		std::pair<int, int> GetPos() const {return (std::pair<int, int>(pos_x,pos_y));};
		int GetXPos() const {return pos_x;};
		int GetYPos() const {return pos_y;};
		void SetValue(TypeValue _value){value = _value;};
		TypeValue GetValue() const {return value;};
		void SetValues(int _x, int _y, TypeValue _value, TypeLabel _label){pos_x = _x; pos_y = _y; value = _value; label = _label;};
		void SetLabel(TypeLabel _label){label = _label;};
		TypeLabel GetLabel() const {return label;};
		void Reset(){pos_x = pos_y = value = label = 0;};
	 private:
		int pos_x;
		int pos_y;
		TypeValue value;
		TypeLabel label;
 };






#endif
