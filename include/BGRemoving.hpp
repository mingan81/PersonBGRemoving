/*! @file BGRemoving.hpp
 *  @brief BGRemoving
 */
#ifndef BGREMOVING_HPP
#define BGREMOVING_HPP

#include <list>
#include <vector>
#include <algorithm>  
using namespace std;
using namespace cv;
using namespace boost::filesystem;



/*! 
 * @brief Makes the comparison op between each imIn pixel and the value seuil. If true, the result is the corresponding pixel of image input or the value input. If false, the result is valNo. 
 * @param imIn : input image
 * @param op : operator
 * @param method : sis or sss or ssi
 * @param parameter : depends on methods
 * @param imOut : output image
 * @todo when new method input, modify brief
 */
template<typename Tin, typename Tout> void ImCompare(const DipImage<Tin> & imIn, const OperationType op, imCompareMethod method, std::vector<ParamT> parameter, DipImage<Tout> & imOut){
	
	dbi.WriteEnter("Entering ImCompare\n");	
	DipImage<unsigned char> imMask(imIn.width(), imIn.height()); imMask.init(255);
	ImCompare(imIn, imMask, op, method, parameter, imOut);
	dbi.WriteOut("Leaving ImCompare\n");
}


/*! 
 * @brief Makes the comparison op between each imIn pixel and the value seuil. If true, the result is the corresponding pixel of image input or the value input. If false, the result is valNo. 
 * @param imIn : input image
 * @param imMask : mask image
 * @param op : operator
 * @param method : sis or sss or ssi
 * @param parameter : depends on methods
 * @param imOut : output image
 * @todo when new method input, modify brief
 */
template<typename Tin, typename Tmask, typename Tout> void ImCompare(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const OperationType op, imCompareMethod method, std::vector<ParamT> parameter, DipImage<Tout> & imOut){
	dbi.WriteEnter("Entering ImCompare\n");
	int w = imIn.width(), h = imIn.height();
	
	switch(method){
		case sis:{
			if(parameter.size() !=3){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(Tin) && parameter[1].type() != typeid(DipImage<Tout>) && parameter[2].type() != typeid(Tout) ){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			Tin seuil = boost::get<Tin>(parameter[0]);
			DipImage<Tout> imValYes(w,h); imValYes.init(0);
			imValYes = boost::get<DipImage<Tout> > (parameter[1]);
			//ImCopy(boost::get<DipImage<Tout> > (parameter[1]), imValYes);
			Tout valNo = boost::get<Tout>(parameter[2]);
			int i, j;
			for (i=0; i<w; i++){
				for (j=0; j<h; j++){
					if (imMask(i,j)!=0)
					switch(op){
						case Equal:
							imOut(i,j) = (imIn(i,j) == seuil) ? imValYes(i,j) : valNo;
							break;
						case SupEqual:
							imOut(i,j) = (imIn(i,j) >= seuil) ? imValYes(i,j) : valNo;
							break;
						case InfEqual:
							imOut(i,j) = (imIn(i,j) <= seuil) ? imValYes(i,j) : valNo;
							break;
						case Sup:
							imOut(i,j) = (imIn(i,j) > seuil) ? imValYes(i,j) : valNo;
							break;
						case Inf:
							imOut(i,j) = (imIn(i,j) < seuil) ? imValYes(i,j) : valNo;
							break;
						default:
							printf("Wrong operator\n");
					}
				}
			}
		}
		break;
		case sss:{
			if(parameter.size() !=3){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(Tin) && parameter[1].type() != typeid(Tout) && parameter[2].type() != typeid(Tout) ){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			Tin seuil = boost::get<Tin>(parameter[0]);
			Tout valYes = boost::get<Tout>(parameter[1]);
			Tout valNo = boost::get<Tout>(parameter[2]);
			int i, j;
				for (i=0; i<w; i++){
					for (j=0; j<h; j++){
						if (imMask(i,j)!=0)
						switch(op){
							case Equal:
								imOut(i,j) = (imIn(i,j) == seuil) ? valYes : valNo;
								break;
							case SupEqual:
								imOut(i,j) = (imIn(i,j) >= seuil) ? valYes : valNo;
								break;
							case InfEqual:
								imOut(i,j) = (imIn(i,j) <= seuil) ? valYes : valNo;
								break;
							case Sup:
								imOut(i,j) = (imIn(i,j) > seuil) ? valYes : valNo;
								break;
							case Inf:
								imOut(i,j) = (imIn(i,j) < seuil) ? valYes : valNo;
								break;
							default:
								printf("Wrong operator\n");
						}
					}
				}
		}
		break;
		case ssi:{
			if(parameter.size() !=3){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(Tin) && parameter[1].type() != typeid(Tout) && parameter[2].type() != typeid(DipImage<Tout>) ){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			Tin seuil = boost::get<Tin>(parameter[0]);
			Tout valYes = boost::get<Tout>(parameter[1]);
			DipImage<Tout> imValNo(w,h); imValNo.init(0);
			imValNo = boost::get<DipImage<Tout> > (parameter[2]);
			int i, j;
			for (i=0; i<w; i++){
				for (j=0; j<h; j++){
					if (imMask(i,j)!=0)
					switch(op){
						case Equal:
							imOut(i,j) = (imIn(i,j) == seuil) ? valYes : imValNo(i,j);
							break;
						case SupEqual:
							imOut(i,j) = (imIn(i,j) >= seuil) ? valYes : imValNo(i,j);
							break;
						case InfEqual:
							imOut(i,j) = (imIn(i,j) <= seuil) ? valYes : imValNo(i,j);
							break;
						case Sup:
							imOut(i,j) = (imIn(i,j) > seuil) ? valYes : imValNo(i,j);
							break;
						case Inf:
							imOut(i,j) = (imIn(i,j) < seuil) ? valYes : imValNo(i,j);
							break;
						default:
							printf("Wrong operator\n");
					}
				}
			}
		}
		break;
		case iii:{
			if(parameter.size() !=3){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(DipImage<Tin>) && parameter[1].type() != typeid(DipImage<Tin>) && parameter[2].type() != typeid(DipImage<Tin>) ){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			DipImage<Tin> imCompare(w,h); imCompare.init(0);
			imCompare = boost::get<DipImage<Tin> > (parameter[0]);
			DipImage<Tin> imValYes(w,h); imValYes.init(0);
			imValYes = boost::get<DipImage<Tin> > (parameter[1]);
			DipImage<Tin> imValNo(w,h); imValNo.init(0);
			imValNo = boost::get<DipImage<Tin> > (parameter[2]);
			int i,j;
			for (i=0; i<w; i++){
				for (j=0; j<h; j++){
					if (imMask(i,j)!=0){
						switch(op)
						{
							case Equal:
								imOut(i,j) = (imIn(i,j) == imCompare(i,j)) ? Tout(imValYes(i,j)) : Tout(imValNo(i,j));
								break;
							case SupEqual:
								imOut(i,j) = (imIn(i,j) >= imCompare(i,j)) ? Tout(imValYes(i,j)) : Tout(imValNo(i,j));
								break;
							case InfEqual:
								imOut(i,j) = (imIn(i,j) <= imCompare(i,j)) ? Tout(imValYes(i,j)) : Tout(imValNo(i,j));
								break;
							case Sup:
								imOut(i,j) = (imIn(i,j) > imCompare(i,j)) ? Tout(imValYes(i,j)) : Tout(imValNo(i,j));
								break;
							case Inf:
								imOut(i,j) = (imIn(i,j) < imCompare(i,j)) ? Tout(imValYes(i,j)) : Tout(imValNo(i,j));
								break;
							default:
								printf("Wrong operator\n");
						}
					}
				}
			}
		}
		break;
	}
	dbi.WriteOut("Leaving ImCompare\n");
}


/*! 
 * @brief RGB image convert to grayscale image based on luminosity
 * @param imIn
 * @param imOut
 */
template<typename T> void RGBToGrayscale_Luminosity(const DipColorImage<T> & imIn, DipImage<T> & imOut)
{
  const int w = imIn.width(), h = imIn.height();
  const T * p_in1 = &imIn(0,0)[0];
  const T * p_in2 = &imIn(0,0)[1];
  const T * p_in3 = &imIn(0,0)[2];
  const T * p_in1end = p_in1 + 3 * w * h;
  T * p_out = &imOut(0,0);

  for(; p_in1 != p_in1end; p_in1+=3, p_in2+=3, p_in3+=3, p_out++)
  {
    *p_out = T(0.21 * *p_in1 + 0.71 * *p_in2 + 0.07 * *p_in3);
  }
}



/*! 
 * @brief Generate binary image from grayscale image by returnin the max value for all non-zero pixels
 * @param imIn : input grayscale image
 * @param [out] imOut : output image
 */
template<typename Tin, typename Tout> void ImBinarisation(const DipImage<Tin> & imIn, DipImage<Tout> & imOut)
{
	// dbi.WriteEnter("Entering ImBinarisation\n");
	const int w = imIn.width(), h = imIn.height();
	Tout maxiout = std::numeric_limits<Tout>::max();
	const Tin * p_in = &imIn(0,0);
	const Tin * p_inend = p_in + w*h;
	Tout * p_out = &imOut(0,0);

	for (; p_in != p_inend; p_in++, p_out++){
		*p_out = (*p_in == 0) ? 0 : maxiout;
	}
	// dbi.WriteOut("Leaving ImBinarisation\n");
}

template<typename T> T ListeAutoThresh(const std::list<T> & ValueListe)
{
	dbi.WriteEnter("\nEntering GetAutoThreshold \n");
	std::list<T> Liste = ValueListe;
	Liste.sort();
	//~ std::list<T> UniqueList;
	//~ UniqueList.clear();
	//~ UniqueList = Liste;
	//~ UniqueList.unique();
	
	//~ T NewThresh((UniqueList.front()+UniqueList.back())/2);
	T NewThresh((Liste.front()+Liste.back())/2);
	//~ dbi.WriteInfo("Min Max class : [ %d - %d ]\n",UniqueList.front(),UniqueList.back());
	//~ dbi.WriteInfo("Min Max class : [ %d - %d ]\n",Liste.front(),Liste.back());
	//~ dbi.WriteInfo("NewThresh = %d ", NewThresh);
	T Thresh(0);
	double sumf1, sumf2;
	int count1, count2;
	
	typename std::list<T>::iterator it = Liste.begin(), itend = Liste.end();
	do
	{
		Thresh = NewThresh;
		sumf1 = 0; count1 = 0;
		sumf2 = 0; count2 = 0;
		for (it = Liste.begin(); it!=itend; it++ )		
			if ((*it)<Thresh)
			{
				sumf1+= (*it);
				count1++;
			}
			else if ((*it)>=Thresh)
			{
				sumf2+= (*it);
				count2++;
			}
		NewThresh = (T)(((sumf1/count1)+(sumf2/count2))/2);
		//~ dbi.WriteInfo("%d ", NewThresh);
	}
	while (NewThresh!=Thresh);
	dbi.WriteInfo("Automatic Threshold : %d\n", NewThresh);
	dbi.WriteOut("Leaving GetAutoThreshold \n");
	return(NewThresh);
}

template <typename Tin, typename Tmask, typename Tout> void ImLocalAutoThresh(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("\nEntering ImLocalAutoThresh \n");
	int w = imIn.width(), h = imIn.height();
	imOut.init(0);
	Neighborhood<int> nh;
	//~ Disc2DNListGenerator(200, nh);
	Square2DNListGenerator(100, nh);
	Neighborhood<int>::iterator it = nh.begin(), itend = nh.end();
	
	DipImage<Tin> imAutoThresh(w,h); imAutoThresh.init(0);
	uint NTotal= ImMeasVolume(imMask)/255;
	uint inc(0);
	Tin AutoThresh(0);
	std::list<Tin> ValueList, UniqueList;
	
	// For each pixel in the ROI
	int i, j, x, y;
	for (j=0; j<h; j++)
		for (i=0; i<w; i++)
			if ((imMask(i,j)!=0))
			{
				ValueList.clear();
				// Go through neighborhood anf fill intensity list
				for (it = nh.begin(); it!=itend; it++)
				{
					x = i+it->first;
					y = j+it->second;
					if ((x>=0)&&(x<w)&&(y>=0)&&(y<h))
						if (imMask(x,y)!=0)
							ValueList.push_back(imIn(x,y));
				}
								
				// Compute automatic threshold
				inc++;
				dbi.WriteInfo("On en est à %f /100 des pixels.\n", double(inc)*100/NTotal);
				AutoThresh = ListeAutoThresh(ValueList);
				imAutoThresh(i,j) = AutoThresh;
				
				// Threshold current pixel value
				if (imIn(i,j)<(AutoThresh)) 
						imOut(i,j) = std::numeric_limits<Tout>::max();
			}
	dbi.WriteImage(imAutoThresh,"imAutoThresh.png");
	dbi.WriteOut("Leaving ImLocalAutoThresh \n");
}


template <typename Tin, typename Tmask, typename Tout> void ImAutoThresh(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("\nEntering ImLocalAutoThresh \n");
	int w = imIn.width(), h = imIn.height();
	imOut.init(0);
	
	uint inc(0);
	Tin AutoThresh(0);
	std::list<Tin> ValueList, UniqueList;
	
	// For each pixel in the ROI
	int i, j, x, y;
	for (j=0; j<h; j++){
		for (i=0; i<w; i++){
			if ((imMask(i,j)!=0))
			{
				ValueList.push_back(imIn(x,y));
			}
		}
	}
	AutoThresh = ListeAutoThresh(ValueList);
	
	// Threshold current pixel value
	if (imIn(i,j)>(AutoThresh)){ 
		imOut(i,j) = std::numeric_limits<Tout>::max();
	}
	dbi.WriteOut("Leaving ImLocalAutoThresh \n");
}


/*! 
 * @brief Performs labeling with respect to a given criterion
 * @param imFZ : image to be labeled. Background = 0
 * @param imMask : mask image 
 * @param method : criterion to be assigned to output image
 * @param parameter : data image to compute criterion or labelling starts from input value
 * @param imLabel : label image. Background = 0
 */
template<typename Tflatzones, typename Tlabel> int ImLabel(const DipImage<Tflatzones> & imFZ,  const DipImage<unsigned char> & imMask, ZoneInfoType method, std::vector<ParamT> parameter, DipImage<Tlabel> & imLabel)
{
	//~ dbi.WriteEnter("Entering ImLabel WithMask\n");
	int w = imFZ.width(), h = imFZ.height();


	imLabel.init(0);
// 	dbi.WriteInfo("Fin init imLabel\n");

	if(method == iteration){
		if(parameter.size() !=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
		}
		if(parameter[0].type() != typeid(TLabel)){
			//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
		}
		
		ParamT parametertype = parameter[0];
		TLabel current_label = boost::get<TLabel>(parametertype);
		
		int i=0, j=0, ii=0, jj=0, ilocal = 0, jlocal = 0;
		int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		int rayon = 1;
		// 	dbi.WriteInfo("debug 0.5\n");
		std::queue<Node<Tflatzones, TLabel> > q;
		// 	dbi.WriteInfo("debug 1\n");
		// 	dbi.WriteInfo("here 1\n");
		for (j=0; j<h; j++){
			for (i=0; i<w; i++){
				if (imFZ(i,j) != 0){
					if(imMask(i,j) != 0){
						if (imLabel(i,j) == Tlabel(0)){
							// dbi.WriteInfo("debug 1.5\n");
							Node<Tflatzones, TLabel> noeud(i,j,Tflatzones(0),current_label);
							q.push(noeud);
							imLabel(i,j) = Tlabel(current_label);
							// 					dbi.WriteInfo("debug 2\n");
							// 			dbi.WriteInfo("here 2\n");
							while (!q.empty()){
								Node<Tflatzones, TLabel> current_node = q.front();
								q.pop();
								ilocal = current_node.GetXPos();
								jlocal = current_node.GetYPos();
								borneInfx = std::max(0,ilocal-rayon);

								borneInfy = std::max(0,jlocal-rayon);
								borneMaxx = std::min(w, ilocal+rayon+1);
								borneMaxy = std::min(h, jlocal+rayon+1);
								if (borneMaxx > w) dbi.WriteInfo("pb largeur\n");
								if (borneMaxy > h) dbi.WriteInfo("pb hauteur\n");
								for (ii=borneInfx; ii<borneMaxx; ii++){
									for (jj=borneInfy; jj<borneMaxy; jj++){
										if((imFZ(ii,jj) != 0) && (imLabel(ii,jj) == Tlabel(0))){
											Node<Tflatzones, TLabel> noeudlocal(ii,jj,0,current_label);
											q.push(noeudlocal);
											imLabel(ii,jj) = Tlabel(current_label);
										}
									}
								}
							}
							// 			dbi.WriteInfo("here 3\n");
							current_label++;
							// 					dbi.WriteInfo("debug 3\n");
						}
					}
				}
			}
		}
		//~ dbi.WriteOut("Leaving ImLabel WithMask\n");
		return (int(current_label-1));
	}
	else if(method == IterationWithNeighbour){
		if(parameter.size() !=2){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
		}
		if(parameter[0].type() != typeid(int) && parameter[1].type() != typeid(Neighborhood<int>) ){
			//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
		}
		
		int firstlabel = boost::get<int>(parameter[0]);
		std::list<std::pair<int, int> > nh = boost::get<std::list<std::pair<int, int> > >(parameter[1]);
		
		DipImage<uchar> imTmp(w,h); imTmp.init(0);
		uchar * p_tmpstart= &imTmp(0,0);
		uchar * p_tmp = &imTmp(0,0);
		uchar * p_tmpnh = &imTmp(0,0);

		const Tflatzones * p_instart = &imFZ(0,0);
		const Tflatzones * p_inx = &imFZ(0,0);
		const Tflatzones * p_iny = &imFZ(0,0);
		const Tflatzones * p_intmp = &imFZ(0,0);
		const Tflatzones * p_inend = p_instart + w*h;
		const Tflatzones * p_in = p_instart;
		const Tflatzones * p_inlinestart = p_instart;
		const Tflatzones * p_inlineend = p_instart + w;
	
		const unsigned char * p_maskstart = &imMask(0,0);
		const unsigned char * p_mask = &imMask(0,0);
		const unsigned char * p_masktmp = p_mask;

		Tlabel * p_outstart = &imLabel(0,0);
		Tlabel * p_out = &imLabel(0,0);
		// Tout * p_outfin = &imLabel(0,0) + w*h;

		int absoffset = 0;

		int currentlabel = firstlabel;

		std::list<int> centers; centers.clear();
		std::list<int> neigh; neigh.clear();
		std::list<int>::iterator itl, itlend;

		std::list<std::pair<int, int> >::iterator itnh, itnhend = nh.end();

		for (;p_in != p_inend; p_in++, p_tmp++, p_mask++){ // Pour tous les points de l'image
			if ((*p_mask) != 0){ // Si le point appartient au masque
				if ((*p_in) != 0){ // Si le point est non nul
					if ((*p_tmp) == 0){ // Si le point n'est pas marque
						absoffset = p_in - p_instart; // on calcule l'offset
						// if (absoffset == 0){ dbi.WriteInfo("Offset = %d\n", absoffset); }
						centers.push_back(absoffset); // on empile le point
						*p_tmp = 10; // on marque le point
						while(!centers.empty()){ // tant que la file de points de centre n'est pas vide
							itl = centers.begin(); // on recale les iterateurs
							itlend = centers.end();
							for (;itl != itlend; itl++){ // pour tous les points de la liste
								absoffset = *itl;
								// dbi.WriteInfo("\toffset = %d\n", absoffset);
								p_out = p_outstart + absoffset;
								// if (absoffset >= 9){
								//	dbi.WriteInfo("p_out en dehors de l'image - offset = %d - maxoffset = %d - w * h = %d\n", absoffset, p_outfin-p_outstart, w*h);}
								// else
								*p_out = currentlabel;
								p_intmp = p_instart + absoffset;
								//dbi.WriteInfo("offset = %d\n", absoffset);
								p_inlinestart = p_instart + int(w * floor(double(absoffset) / double(w)));
								// dbi.WriteInfo("offset = %d\n", (p_inlinestart-p_instart));
								p_inlineend = p_inlinestart + w;
								itnh = nh.begin();
								for (;itnh != itnhend; itnh++){
									// dbi.WriteInfo("Offset du voisinage = (%d,%d)\n", (*itnh).first, (*itnh).second);
									p_inx = p_intmp + (*itnh).first;
									// dbi.WriteInfo("\tp_inx = %d\n", p_inx - p_instart);
									if ((p_inx >= p_inlinestart) && (p_inx < p_inlineend)){
										p_iny = p_inx + w * (*itnh).second;
										if (p_iny < p_inend){
											absoffset = p_iny - p_instart;
											p_masktmp = p_maskstart + absoffset;
											p_tmpnh = p_tmpstart + absoffset;
											if ((*p_masktmp) != 0){
												if ((*p_iny) != 0){
													if ((*p_tmpnh) == 0){
														neigh.push_back(absoffset);
														// dbi.WriteInfo("\t\toffset neigh = %d\n", absoffsettmp);
														*p_tmpnh = 10;
													}
												}
											}
										}
									}
								}
							}
							centers.clear();
							centers.swap(neigh);
							neigh.clear();
						}
						currentlabel++;
					}
				}
			}
		}

	// 	dbi.WriteInfo("Valeur de retour = %d\n", currentlabel);
	//	dbi.WriteImage(imOut, "imLabel.png");
	//	dbi.WriteImage(imTmp, "imTmp.png");
	//	dbi.WriteImage(imIn, "imOrig.png");
	//~ dbi.WriteOut("Leaving ImLabel WithMask\n");
	return(currentlabel-1);
	}
	else 
	{
		if(parameter.size() !=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
		}
		if(parameter[0].type() != typeid(DipImage<double>)){
			//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
		}
		
		ParamT parametertype = parameter[0];
		DipImage<double> imData = boost::get<DipImage<double> >(parametertype);
	
		int i=0, j=0, ii=0, jj=0;
		DipImage<unsigned char> flag(w,h);
		// 	dbi.WriteInfo("Debut init de flag\n");
		flag.init(0);
		// 	dbi.WriteInfo("Init de flag terminee\n");

		int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		int rayon = 1;
		Tflatzones current_label = 0;
		int counter = 1;

		double maxival = std::numeric_limits<double>::max();
		double current_value = 0;

		std::list<std::pair<int, int> > regionpoints; if(!regionpoints.empty()) printf("The list regionpoints is not empty\n");
		std::list<std::pair<int, int> > neighborpoints; if(!neighborpoints.empty()) printf("The list neighborpoints is not empty\n");

		FZInfoNode<double> fzinode;

		std::pair<int, int> temppaire(0,0);
		std::pair<int, int> local(0,0);


		for (i=0; i<w; i++){
			for (j=0; j<h; j++){
			if (imMask(i,j) != 0){
				if (imFZ(i,j) != 0){
					if (flag(i,j) == 0){
// 						dbi.WriteInfo("area = %d\n", fzinode.area);
						flag(i,j) = 1;
						current_label = imFZ(i,j);
						fzinode.regionnumber = counter;
						counter++;

						temppaire.first = i; temppaire.second = j;
						regionpoints.push_back(temppaire);
						neighborpoints.push_back(temppaire);

	// 					imLabel(i,j) = fzinode.regionnumber;

						current_value = imData(i,j);
						IsMinOrMax(i, fzinode.x_min, fzinode.x_max);
						IsMinOrMax(j, fzinode.y_min, fzinode.y_max);
						IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
						fzinode.area = fzinode.area + 1;
						fzinode.average += current_value;
						fzinode.volume += (long int)(maxival - current_value);
						if (current_value != 0)
							fzinode.nbptnoz = 1;
						else
							fzinode.nbptnoz = 0;

						while(!neighborpoints.empty()){
							local = neighborpoints.front();
							neighborpoints.pop_front();

							borneInfx = std::max(0,local.first - rayon);
							borneInfy = std::max(0,local.second - rayon);
							borneMaxx = std::min(w, local.first + rayon + 1);
							borneMaxy = std::min(h, local.second + rayon + 1);

							for (ii=borneInfx; ii<borneMaxx; ii++){
								for (jj=borneInfy; jj<borneMaxy; jj++){
									if (imMask(ii,jj) != 0){
										if((imFZ(ii,jj) == current_label) && (flag(ii,jj) == 0)){
											if ((ii==i)&&(jj==j)) printf("Pb: iteration sur le point central\n");

											temppaire.first = ii; temppaire.second = jj;
											regionpoints.push_back(temppaire);
											neighborpoints.push_back(temppaire);

// 											imLabel(ii,jj) = fzinode.regionnumber;
											flag(ii,jj) = 1;
											current_value = imData(ii,jj);
											IsMinOrMax(ii, fzinode.x_min, fzinode.x_max);
											IsMinOrMax(jj, fzinode.y_min, fzinode.y_max);
											IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
											fzinode.area = fzinode.area + 1;
											fzinode.average += current_value;
											fzinode.stddev += double(current_value * current_value);
											fzinode.volume += (long int)(maxival - current_value);
											if (current_value != 0)
												fzinode.nbptnoz = fzinode.nbptnoz + 1;
										}
									}
								}
							}
						}
						fzinode.average = double(fzinode.average) / double(fzinode.area);
						fzinode.stddev = std::sqrt(double(fzinode.stddev) / double(fzinode.area) - fzinode.average * fzinode.average);
						fzinode.dynamic = int((fzinode.val_max - fzinode.val_min));
						fzinode.volume = (long int)(fzinode.volume - fzinode.area*(maxival - fzinode.val_max));

						std::list<std::pair<int, int> >::iterator ittemp = regionpoints.begin();
						std::list<std::pair<int, int> >::iterator itend = regionpoints.end();

						// dbi.WriteInfo("Valeur de la moyenne = %g\n", double(fzinode.average));

						switch (method){
							case Width:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.x_max - fzinode.x_min);
								break;
							case Height:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.y_max - fzinode.y_min);
								break;
							case Val_min:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.val_min);
								break;
							case Val_max:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.val_max);
								break;
							case Average:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.average);
								break;
							case Dynamic:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.dynamic);
								break;
							case Area:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.area);
								break;
							case Volume:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.volume);
								break;
							case StdDev:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.stddev);
								break;
							case NbPtNoZ:
								for(;ittemp != itend; ittemp++)
									imLabel(ittemp->first, ittemp->second) = Tlabel(fzinode.nbptnoz);
								break;
							default:
								printf("Sorry, this feature is not available yet\n");
						}
						regionpoints.clear();
						neighborpoints.clear();
						fzinode.Reset();
						}
					}
				}
			}
		}
		//~ dbi.WriteOut("Leaving ImLabel WithMask\n");
		return 0;
	}
	
	
	
}

/*! 
 * @brief Performs labeling with respect to a given criterion without mask
 * @param imFZ : image to be labeled. Background = 0 
 * @param method : criterion to be assigned to output image
 * @param parameter : data image to compute criterion or labelling starts from input value
 * @param imLabel : label image. Background = 0
 */
template<typename Tflatzones, typename Tout> int ImLabel(const DipImage<Tflatzones> & imFZ, ZoneInfoType method, std::vector<ParamT> parameter, DipImage<Tout> & imLabel)
{
	dbi.WriteEnter("Entering ImLabel\n");
	DipImage<unsigned char> imMask(imFZ.width(), imFZ.height()); imMask.init(255);
	int result;
	result = ImLabel(imFZ,imMask,method,parameter,imLabel);	
	dbi.WriteOut("Leaving ImLabel\n");
	return result;
}




/*! @fn void ImLabelFlatZones_WithCriterion(const BasicImage<Tflatzones> & imFZ, const BasicImage<Tdata> & imData, const ZoneInfoType criterion, BasicImage<Tout> & imLabel)
 * @brief Performs flat zones labeling with respect to a given criterion
 * @param imFZ : image to be labeled. Background = 0
 * @param imData : data image to compute criterion
 * @param criterion : criterion to be assigned to output image
 * @param imLabel : label image. Background = 0
 */
template<typename Tflatzones, typename Tdata, typename Tout> void ImLabelFlatZones_WithCriterion(const DipImage<Tflatzones> & imFZ, const DipImage<Tdata> & imData, const ZoneInfoType criterion, DipImage<Tout> & imLabel)
{
	dbi.WriteEnter("Entering ImLabelFlatZones_WithCriterion\n");
	int w = imFZ.width(), h = imFZ.height();
// 	dbi.WriteInfo("Largeur = %d\n", w);
// 	dbi.WriteInfo("Hauteur = %d\n", h);
// 	dbi.WriteInfo("Critere = %d\n", criterion);

	imLabel.init(0);
// 	dbi.WriteInfo("Fin init imLabel\n");

	int i=0, j=0, ii=0, jj=0;
	DipImage<unsigned char> flag(w,h);
// 	dbi.WriteInfo("Debut init de flag\n");
	flag.init(0);
// 	dbi.WriteInfo("Init de flag terminee\n");

	int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
	int rayon = 1;
	Tflatzones current_label = 0;
	int counter = 1;

	Tdata maxival = std::numeric_limits<Tdata>::max();
	Tdata current_value = 0;

	std::list<std::pair<int, int> > regionpoints; if(!regionpoints.empty()) printf("The list regionpoints is not empty\n");
	std::list<std::pair<int, int> > neighborpoints; if(!neighborpoints.empty()) printf("The list neighborpoints is not empty\n");

	FZInfoNode<Tdata> fzinode;

	std::pair<int, int> temppaire(0,0);
	std::pair<int, int> local(0,0);


	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if (imFZ(i,j) != 0){
				if (flag(i,j) == 0){
 					dbi.WriteInfo("area = %d\n", fzinode.area);
					flag(i,j) = 1;
					current_label = imFZ(i,j);
					fzinode.regionnumber = counter;
					counter++;

					temppaire.first = i; temppaire.second = j;
					regionpoints.push_back(temppaire);
					neighborpoints.push_back(temppaire);

// 					imLabel(i,j) = fzinode.regionnumber;

					current_value = imData(i,j);
					IsMinOrMax(i, fzinode.x_min, fzinode.x_max);
					IsMinOrMax(j, fzinode.y_min, fzinode.y_max);
					IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
					fzinode.area = fzinode.area + 1;
					fzinode.average += current_value;
					fzinode.volume += (long int)(maxival - current_value);

					while(!neighborpoints.empty()){
						local = neighborpoints.front();
						neighborpoints.pop_front();

						borneInfx = std::max(0,local.first - rayon);
						borneInfy = std::max(0,local.second - rayon);
						borneMaxx = std::min(w, local.first + rayon + 1);
						borneMaxy = std::min(h, local.second + rayon + 1);

						for (ii=borneInfx; ii<borneMaxx; ii++){
							for (jj=borneInfy; jj<borneMaxy; jj++){
								if((imFZ(ii,jj) == current_label) && (flag(ii,jj) == 0)){
									if ((ii==i)&&(jj==j)) printf("Pb: iteration sur le point central\n");

									temppaire.first = ii; temppaire.second = jj;
									regionpoints.push_back(temppaire);
									neighborpoints.push_back(temppaire);

// 									imLabel(ii,jj) = fzinode.regionnumber;
									flag(ii,jj) = 1;
									current_value = imData(ii,jj);
									IsMinOrMax(ii, fzinode.x_min, fzinode.x_max);
									IsMinOrMax(jj, fzinode.y_min, fzinode.y_max);
									IsMinOrMax(current_value, fzinode.val_min, fzinode.val_max);
									fzinode.area = fzinode.area + 1;
									fzinode.average += current_value;
									fzinode.stddev += double(current_value * current_value);
									fzinode.volume += (long int)(maxival - current_value);
								}
							}
						}
						//fzinode.DisplayNode();
					}
					fzinode.average = double(fzinode.average) / double(fzinode.area);
					fzinode.stddev = std::sqrt(double(fzinode.stddev) / double(fzinode.area) - fzinode.average * fzinode.average);
					fzinode.dynamic = Tdata((fzinode.val_max - fzinode.val_min));
					fzinode.volume = (long int)(fzinode.volume - fzinode.area*(maxival - fzinode.val_max));

					std::list<std::pair<int, int> >::iterator ittemp = regionpoints.begin();
					std::list<std::pair<int, int> >::iterator itend = regionpoints.end();

					switch (criterion)
					{
						case Width:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.x_max - fzinode.x_min);
							break;
						case Height:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.y_max - fzinode.y_min);
							break;
						case Val_min:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.val_min);
							break;
						case Average:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.average);
							break;
						case Dynamic:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.dynamic);
							break;
						case Area:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.area);
							break;
						case Volume:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.volume);
							break;
						case StdDev:
							for(;ittemp != itend; ittemp++)
								imLabel(ittemp->first, ittemp->second) = Tout(fzinode.stddev);
							break;
						default:
							printf("Sorry, this feature is not available yet\n");
					}
					dbi.WriteInfo("area = %d\n", fzinode.area);
					regionpoints.clear();
					neighborpoints.clear();
					fzinode.Reset();
				}
			}
		}
	}
	dbi.WriteInfo("Nombre de regions traitees = %d\n", counter-1);
	dbi.WriteOut("Leaving ImLabelFlatZones_WithCriterion\n");
}


/**
 * ImNormalisation_WithMask
 * @param imIn
 * @param imMask
 * @param imOut
 */
template<typename Tin,typename Tmask,typename Tout> void ImNormalisation(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImNormalisation_WithMask\n");
	const int w = imIn.width(), h = imIn.height();
	
	std::vector<double> meanStddev = MeasStatistics(imIn,imMask,all);

	int i,j;
	for(i=0;i<w;i++){
		for(j=0;j<h;j++){
			if(imMask(i,j)!=0){
				imOut(i,j) = Tout(((double)imIn(i,j) - meanStddev[0]) / meanStddev[1]);
				if(abs(imOut(i,j)) < std::numeric_limits<Tout>::epsilon()) imOut(i,j) = Tout(0); 
			}
		}
	}
	dbi.WriteOut("Leaving ImNormalisation_WithMask\n");
}

/*! 
 * @brief Computes the volume of the input image restrained to the given mask. volume = Sum(Grey level)
 * @param imIn : input image
 * @param imMask : mask information image
 * @return the resulting volume
 */
template<typename Tin, typename Tmask> uint ImMeasVolume(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask)
{
	
		const int w = imIn.width(), h = imIn.height();
		const Tin * p_imIn = &imIn(0,0);
		const Tin * p_imInEnd = p_imIn + w*h;
		const Tmask * p_imMask = &imMask(0,0);
		uint volume = 0;
    
		for (;p_imIn != p_imInEnd; p_imIn++, p_imMask++){
			if ((*p_imMask) != 0){
				volume += uint(*p_imIn);
			}
		}    
	
		return (volume);
}



/*! 
 * @brief Computes the volume of the input image. volume = Sum(Grey level)
 * @param imIn : input image
 * @return : the resulting volume
 */
template<typename T> uint ImMeasVolume(const DipImage<T> & imIn)
{
 
    DipImage<uchar> imMask(imIn.width(), imIn.height());
    imMask.init(255);
    return ImMeasVolume(imIn, imMask); 
}



/*! 
 * @brief Computes the minimum and maximum values in the input image with respect to the mask
 * @param imIn : input image
 * @param imMask : mask image
 * @return a pair containing respectively the minimum and the maximum values
 */
template<typename T, typename Tmask> std::pair<T, T> MeasMinMax(const DipImage<T> & imIn, const DipImage<Tmask> & imMask)
{
	dbi.WriteEnter("Entering MeasMinMax toto\n");
	// 	minmax.first = std::numeric_limits<T>::max(); minmax.second = std::numeric_limits<T>::min();
	T mini = std::numeric_limits<T>::max();
	T maxi;
	if (std::numeric_limits<T>::min() != 0){
		maxi = -abs(std::numeric_limits<T>::max());
	}
	else{
		maxi = 0;
	}

	dbi.WriteInfo("mini = %e\n", double(mini));
	dbi.WriteInfo("maxi = %e\n", double(maxi));


	// dbi.WriteImage(imMask, "MeasMinMax_WithMask_imMask.png");
	const int w = imIn.width(), h = imIn.height();
	const T * p_in = &imIn(0,0);
	const T * p_inend = p_in + w*h;
	const Tmask * p_mask = &imMask(0,0);

	int counter = 0;

	for (;p_in != p_inend; p_in++, p_mask++)
	{
		if (*p_mask != 0)
		{
			if (*p_in < mini) mini = *p_in;
			// if (*p_in < 0)
			//	dbi.WriteInfo("value = %e\n", double(*p_in));
			if (*p_in > maxi) maxi = *p_in;
			counter++;
		}
	}
	dbi.WriteInfo("counter = %d\n", counter);
	dbi.WriteInfo("mini = %e\n", double(mini));
	dbi.WriteInfo("maxi = %e\n", double(maxi));
	std::pair<T,T> minmax;
	minmax.first = mini;
	minmax.second = maxi;
	dbi.WriteOut("Leaving MeasMinMax\n");
	return minmax;
}

/*! 
 * @brief Computes the minimum and maximum values in the input image 
 * @param imIn : input image
 * @return a pair containing respectively the minimum and the maximum values
 */
template<typename T> std::pair<T, T> MeasMinMax(DipImage<T> & imIn)
{
	DipImage<unsigned char> imMask(imIn.width(), imIn.height()); imMask.init(255);
	return MeasMinMax(imIn, imMask);
}


template<typename T> void ImColorBandSeparation(const DipColorImage<T> & imIn, DipImage<T> & imCanal1, DipImage<T> & imCanal2, DipImage<T> & imCanal3)
{
  const T * p_in1 = &imIn(0,0)[0];
  const T * p_in1end = p_in1 + 3 * imIn.width() * imIn.height();
  const T * p_in2 = p_in1 + 1;
  const T * p_in3 = p_in1 + 2;
  T * p_canal1 = &imCanal1(0,0);
  T * p_canal2 = &imCanal2(0,0);
  T * p_canal3 = &imCanal3(0,0);

  for (; p_in1 != p_in1end; p_in1+=3, p_in2+=3, p_in3+=3, p_canal1++, p_canal2++, p_canal3++)
  {
    *p_canal1 = *p_in1;
    *p_canal2 = *p_in2;
    *p_canal3 = *p_in3;
  }
}

/*! @fn template<typename Tin, typename Tmask, typename Tin2, typename Tout> void ImAddImage(const BasicImage<Tin> & imIn, const BasicImage<Tmask> & imMask, const ClipMethod method, const BasicImage<Tin2> & imIn2, BasicImage<Tout> & imOut)
 * @brief Operates the summation of imIn1 and imIn2. 
 * @param imIn1 : input image
 * @param imMask : mask image
 * @param method : either Clip or NoClip (with Clip: Truncation is done if result is superior to max value).
 * @param imIn2 : input image (must have the same type than imIn1)
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tmask, typename Tin2, typename Tout> void ImAddImage(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const ClipMethod method, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImAddImage with a Mask\n");
	switch (method)
	{
		case NoClip:
		{
			const Tin * p_in1 = &imIn(0,0);
			const Tin * p_in1end = &imIn(0,0) + imIn.width()*imIn.height();
			const Tin2 * p_in2 = &imIn2(0,0);
			const Tmask * p_mask = &imMask(0,0);
			Tout * p_out = &imOut(0,0);
			
			for (; p_in1 != p_in1end; p_in1++, p_in2++, p_mask++, p_out++){
				if (*p_mask != 0){
						*p_out = (Tout)((double)*p_in1 + (double)*p_in2);
				}
				else{
					*p_out = (Tout) *p_in1;
				}
			}
		break;
		}
		case Clip:
		{
			
				
			const Tin * p_in1 = &imIn(0,0);
			const Tin * p_in1end = &imIn(0,0) + imIn.width()*imIn.height();
			const Tin2 * p_in2 = &imIn2(0,0);
			const Tmask * p_mask = &imMask(0,0);
			Tout * p_out = &imOut(0,0);
			
			for (; p_in1 != p_in1end; p_in1++, p_in2++, p_mask++, p_out++){
				if (*p_mask != 0){
					*p_out = (Tout)std::min(double(*p_in1) + double(*p_in2), double(std::numeric_limits<Tout>::max()));
				}
				else{
					*p_out = (Tout) *p_in1;
				}
			}

			break;
		}
	}
	dbi.WriteOut("Leaving ImAddImage with a Mask.\n");
}

/*! @fn template<typename Tin, typename Tin2, typename Tout> void ImAddImage(const BasicImage<Tin> & imIn, const ClipMethod method, const BasicImage<Tin2> & imIn2, BasicImage<Tout> & imOut)
 * @brief Operates the summation of imIn1 and imIn2. 
 * @param imIn1 : input image
 * @param method : either Clip or NoClip (with Clip: Truncation is done if result is superior to max value).
 * @param imIn2 : input image (must have the same type than imIn1)
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tin2, typename Tout> void ImAddImage(const DipImage<Tin> & imIn, const ClipMethod method, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImAddImage\n");
	int w = imIn.width(), h = imIn.height();
	DipImage<unsigned char> imMask(w,h); imMask.init(255);
	ImAddImage(imIn, imMask, method, imIn2, imOut);
	dbi.WriteOut("Leaving ImAddImage.\n");
}

/*! @fn template<typename Tin, typename Tmask, typename Tconst, typename Tout> void ImDivConst(const BasicImage<Tin> & imIn,  const BasicImage<Tmask> & imMask,const Tconst constante, BasicImage<Tout> & imOut)
 * @brief Operates the division of a imIn1 by a constant. 
 * @param imIn1 : input image
 * @param imMask : mask image
 * @param constante : scalar value to be divided to the input image
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tmask, typename Tconst, typename Tout> void ImDivConst(const DipImage<Tin> & imIn,  const DipImage<Tmask> & imMask,const Tconst constante, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImDivConst with a Mask\n");
	
	int w = imIn.width(), h = imIn.height();
	DipImage<Tconst> imIn2(w,h); imIn2.init(constante);
	ImDivImage(imIn, imMask, imIn2, imOut);
	
	dbi.WriteOut("Leaving ImDivConst  with a Mask.\n");
}

/*! @fn template<typename Tin, typename Tconst, typename Tout> void ImDivConst(const BasicImage<Tin> & imIn, const Tconst constante, BasicImage<Tout> & imOut)
 * @brief Operates the division of a imIn1 by a constant. 
 * @param imIn1 : input image
 * @param constante : scalar value to be divided to the input image
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tconst, typename Tout> void ImDivConst(const DipImage<Tin> & imIn, const Tconst constante, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImDivConst\n");
	
	int w = imIn.width(), h = imIn.height();
	DipImage<Tconst> imIn2(w,h); imIn2.init(constante);
	ImDivImage(imIn, imIn2, imOut);
	
	dbi.WriteOut("Leaving ImDivConst.\n");
}


/*! @fn template<typename Tin, typename Tmask, typename Tin2, typename Tout> void ImDivImage(const BasicImage<Tin> & imIn, const BasicImage<Tmask> & imMask, const BasicImage<Tin2> & imIn2, BasicImage<Tout> & imOut)
 * @brief Operates the division of a imIn1 by imIn2. 
 * @param imIn1 : input image
 * @param imMask : mask image
 * @param imIn2 : input image (must have the same type than imIn1)
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tmask, typename Tin2, typename Tout> void ImDivImage(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImDivImage with a Mask\n");
	//~ try
	//~ {
		//~ vigra_assert(imIn.height()==imIn2.height(), "ImDivImage: Images must have the same height");
		//~ vigra_assert(imIn.width()==imIn2.width(), "ImDivImage: Images must have the same width");
	
		const Tin * p_in = &imIn(0,0);
		const Tin2 * p_in2 = &imIn2(0,0);
		const Tin * p_inend = p_in + imIn.width() * imIn.height();
		const Tmask * p_mask = &imMask(0,0);
		Tout * p_out = &imOut(0,0);

		for (; p_in != p_inend; p_in++, p_in2++, p_mask++, p_out++){
			if (*p_mask != 0){
				*p_out = Tout(double(*p_in) / double(*p_in2));
			}
		}
	//~ }
	//~ catch (vigra::StdException & e)
	//~ {
		//~ // catch any errors that might have occured and print their reason
		//~ std::cout << e.what() << std::endl;
	//~ }
	dbi.WriteOut("Leaving ImDivImage with a Mask.\n");
}

/*! @fn template<typename Tin, typename Tin2, typename Tout> void ImDivImage(const BasicImage<Tin> & imIn, const BasicImage<Tin2> & imIn2, BasicImage<Tout> & imOut)
 * @brief Operates the division of a imIn1 by imIn2. 
 * @param imIn1 : input image
 * @param imIn2 : input image (must have the same type than imIn1)
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tin2, typename Tout> void ImDivImage(const DipImage<Tin> & imIn, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImDivImage\n");
	
	int w = imIn.width(), h = imIn.height();
	DipImage<unsigned char> imMask(w,h); imMask.init(255);
	ImDivImage(imIn, imMask, imIn2, imOut);
	
	dbi.WriteOut("Leaving ImDivImage.\n");
}


/*! 
 * @brief Computes the statistics of the input image with respect to the mask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method: input method (either dip_mean, stddev, variance or all)
 * @return a vector containing all the desired statistics
 * @todo modify it later
 */
template<typename Tin, typename Tmask> std::vector<double> MeasStatistics(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask,statisticsMethod method){
	dbi.WriteEnter("Entering MeasStatistics\n");
	int w = imIn.width(), h = imIn.height();
	int i=0, j=0;
	
	double temp = 0;
	double tempsquare = 0;
	int counter = 0;
	std::vector<double> result;result.clear();
	
	Tin valeur = 0;
	
	double volume = ImMeasVolume(imMask);
	dbi.WriteInfo("nb of points = %d\n", int(volume/255));
	
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if (imMask(i,j) != 0){
				valeur = imIn(i,j);
				temp += double(valeur);
				tempsquare += double(valeur * valeur);
				counter +=1;
			}
		}
	}
	if (counter == 0){
		switch(method){
			case dip_mean:{
				result.push_back(0);
			}
			break;
			case stddev:{
				result.push_back(0);
			}
			break;
			case variance:{
				result.push_back(0);
			}
			break;
			case all:{
				result.push_back(0);
				result.push_back(0);
				result.push_back(0);
			}
			break;
		}
		dbi.WriteOut("Leaving MeasStatistics\n");
		return result;
	}
	else{
		double moyenne = temp / counter;
		double varianceV = abs((tempsquare/counter) - (moyenne * moyenne));
		double stddevV = sqrt(varianceV);
		dbi.WriteInfo("compteur = %d\n", counter);
		switch(method){
			case dip_mean:{
				result.push_back(moyenne);
			}
			break;
			case stddev:{
				result.push_back(stddevV);
			}
			break;
			case variance:{
				result.push_back(varianceV);
			}
			break;
			case all:{
				result.push_back(moyenne);
				result.push_back(stddevV);
				result.push_back(varianceV);
			}
			break;
		}
		dbi.WriteOut("Leaving MeasStatistics\n");
		return result;
	}
}

/*! 
 * @brief Assigns val_yes to imOut point if corresponding imIn point value is in range [thresh_min, thresh_max]. Assigns val_no if not.
 * @param imIn : input image
 * @param thresh_min : low threshold boundary
 * @param thresh_max : high threshold boundary
 * @param val_yes : value assigned to imOut if imIn value belongs to interval [thresh_min, thresh_max]
 * @param val_no : alternative value
 * @param imOut : output image
 */
template<typename Tin, typename Tout> void ImThreshold(const DipImage<Tin> & imIn, const Tin thresh_min, const Tin thresh_max, const Tout val_yes, const Tout val_no, DipImage<Tout> & imOut)
{
  dbi.WriteEnter("Entering ImThreshold\n");
  const Tin * p_in = &imIn(0,0);
  const Tin * p_inend = p_in + imIn.width() * imIn.height() - 1;
  Tout * p_out = &imOut(0,0);
  
  for (;p_in != p_inend; p_in++, p_out++){
      *p_out = ((*p_in <= thresh_max) && (*p_in >= thresh_min)) ?  val_yes : val_no;
  }
  dbi.WriteOut("Leaving ImThreshold\n");
}



void detectBGEstimation(std::vector<DipColorImage<uchar> > & frames, DipImage<uchar> & BGEstimation,DipImage<double> & stddev){
	dbi.WriteEnter("Enter detectBGEstimation");
	int frameNumber = frames.size();
	int w = BGEstimation.width();
	int h = BGEstimation.height();
	
	dbi.WriteInfo("frameNumber = %d",frameNumber);
	
	DipImage<uchar> imIn1(w,h); imIn1.init(0);
	DipImage<uchar> imIn2(w,h); imIn2.init(0);
	DipImage<uchar> imIn3(w,h); imIn3.init(0);
	
	
	DipImage<double> stddev1(w,h); stddev1.init(0);
	DipImage<double> stddev2(w,h); stddev2.init(0);
	DipImage<double> stddev3(w,h); stddev3.init(0);
	DipImage<double> mean1(w,h); mean1.init(0);
	DipImage<double> mean2(w,h); mean2.init(0);
	DipImage<double> mean3(w,h); mean3.init(0);
	DipImage<double> sum1(w,h); sum1.init(0);
	DipImage<double> sum2(w,h); sum2.init(0);
	DipImage<double> sum3(w,h); sum3.init(0);
	DipImage<double> squareSum1(w,h); squareSum1.init(0);
	DipImage<double> squareSum2(w,h); squareSum2.init(0);
	DipImage<double> squareSum3(w,h); squareSum3.init(0);
	int count = 0;

	
	
	for(int i = 0;i < frameNumber; i++){
		if(!frames[i].empty()){
			count = count +1;
			imIn1.init(0);
			imIn2.init(0);
			imIn3.init(0);
			ImColorBandSeparation(frames.at(i),imIn1,imIn2,imIn3);
			//~ dbi.WriteImage(imIn1, "imIn1.png");
			//~ dbi.WriteImage(imIn2, "imIn2.png");
			//~ dbi.WriteImage(imIn3, "imIn3.png");
			for(int x = 0;x < w;x++){
				for(int y = 0;y < h;y++){
					sum1(x,y) = sum1(x,y) + double(imIn1(x,y));
					sum2(x,y) = sum2(x,y) + double(imIn2(x,y));
					sum3(x,y) = sum3(x,y) + double(imIn3(x,y));
					squareSum1(x,y) = squareSum1(x,y) + double(imIn1(x,y)) * double(imIn1(x,y));
					squareSum2(x,y) = squareSum2(x,y) + double(imIn2(x,y)) * double(imIn2(x,y));
					squareSum3(x,y) = squareSum3(x,y) + double(imIn3(x,y)) * double(imIn3(x,y));
					
				}
			}
		}
	}
	
	for(int x = 0;x < w;x++){
		for(int y = 0;y < h;y++){
					
			mean1(x,y) = sum1(x,y)/count;
			mean2(x,y) = sum2(x,y)/count;
			mean3(x,y) = sum3(x,y)/count;
			
			//~ dbi.WriteInfo("mean1 = %f",mean1(x,y));
			//~ dbi.WriteInfo("mean2 = %f",mean2(x,y));
			//~ dbi.WriteInfo("mean3 = %f\n",mean3(x,y));	
			stddev1(x,y) = sqrt(squareSum1(x,y)/count - mean1(x,y) *mean1(x,y));
			stddev2(x,y) = sqrt(squareSum2(x,y)/count - mean2(x,y) *mean2(x,y));
			stddev3(x,y) = sqrt(squareSum3(x,y)/count - mean3(x,y) *mean3(x,y));
			
			
		}
	}
	
	
	for(int x = 0;x < w;x++){
		for(int y = 0;y < h;y++){
			//~ dbi.WriteInfo("stddev1 = %f",stddev1(x,y));
			//~ dbi.WriteInfo("stddev2 = %f",stddev2(x,y));
			//~ dbi.WriteInfo("stddev3 = %f\n",stddev3(x,y));
			if(stddev1(x,y) < 5 && stddev2(x,y) < 5 && stddev3(x,y) < 5){
				BGEstimation(x,y) = uchar(255);
			}
		}
	}
	
	
	dbi.WriteImage(stddev1, "stddev1.png");
	dbi.WriteImage(stddev2, "stddev2.png");
	dbi.WriteImage(stddev3, "stddev3.png");
	
	DipImage<double> stddev12(w,h); stddev12.init(0);
	ImAddImage(stddev1,NoClip,stddev2,stddev12);

	ImAddImage(stddev12,NoClip,stddev3,stddev);

	
	dbi.WriteOut("Leave detectBGEstimation"); 
	
	
}


void detectBGEstimation(std::vector<DipImage<uchar> > & frames, DipImage<double> & stddev1){
	dbi.WriteEnter("Enter detectBGEstimation");
	int frameNumber = frames.size();
	int w = stddev1.width();
	int h = stddev1.height();
	
	dbi.WriteInfo("frameNumber = %d",frameNumber);
	
	
	
	//~ DipImage<double> stddev1(w,h); stddev1.init(0);	
	DipImage<double> mean1(w,h); mean1.init(0);	
	DipImage<double> sum1(w,h); sum1.init(0);	
	DipImage<double> squareSum1(w,h); squareSum1.init(0);	
	int count = 0;

	
	
	for(int i = 0;i < frameNumber; i++){
		if(!frames[i].empty()){
			count = count +1;		
			for(int x = 0;x < w;x++){
				for(int y = 0;y < h;y++){
					sum1(x,y) = sum1(x,y) + double(frames[i](x,y));					
					squareSum1(x,y) = squareSum1(x,y) + double(frames[i](x,y)) * double(frames[i](x,y));				
					
				}
			}
		}
	}
	
	for(int x = 0;x < w;x++){
		for(int y = 0;y < h;y++){
					
			mean1(x,y) = sum1(x,y)/count;		
			stddev1(x,y) = sqrt(squareSum1(x,y)/count - mean1(x,y) *mean1(x,y));			
			
			
		}
	}
	
	
	
	dbi.WriteImage(stddev1, "stddevspecial.png");	
	
	
	
	dbi.WriteOut("Leave detectBGEstimation"); 
	
	
}

template<typename Tin, typename Tmask, typename Tin2, typename Tout> void ImSubImage(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const ClipMethod method, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImSubImage with a Mask\n");
	switch (method)
	{
		case NoClip:
		{
			//~ try
			//~ {
				//~ vigra_assert(imIn.height()==imIn2.height(), "ImSubImage: Images must have the same height");
				//~ vigra_assert(imIn.width()==imIn2.width(), "ImSubImage: Images must have the same width");

				const Tin * p_in1 = &imIn(0,0);
				const Tin * p_in1end = &imIn(0,0) + imIn.width()*imIn.height();
				const Tin2 * p_in2 = &imIn2(0,0);
				const Tmask * p_mask = &imMask(0,0);
				Tout * p_out = &imOut(0,0);
				
				for (; p_in1 != p_in1end; p_in1++, p_in2++, p_mask++, p_out++){
					if (*p_mask != 0){
							*p_out = (Tout)((double)*p_in1 - (double)*p_in2);
					}
					else{
						*p_out = (Tout)0;
						//~ *p_out = (Tout)*p_in1;
					}
				}
			//~ }
			//~ catch (vigra::StdException & e)
			//~ {
				//~ // catch any errors that might have occured and print their reason
				//~ std::cout << e.what() << std::endl;
			//~ }
			break;
		}
		case Clip:
		{
			//~ try
			//~ {
				//~ vigra_assert(imIn.height()==imIn2.height(), "ImSubImage: Images must have the same height");
				//~ vigra_assert(imIn.width()==imIn2.width(), "ImSubImage: Images must have the same width");
				
				const Tin * p_in1 = &imIn(0,0);
				const Tin * p_in1end = &imIn(0,0) + imIn.width()*imIn.height();
				const Tin2 * p_in2 = &imIn2(0,0);
				const Tmask * p_mask = &imMask(0,0);
				Tout * p_out = &imOut(0,0);
				
				for (; p_in1 != p_in1end; p_in1++, p_in2++, p_mask++, p_out++){
					if (*p_mask != 0){
						*p_out = (Tout)std::max(double(*p_in1) - double(*p_in2), 0.0);
					}
					else{
						*p_out = (Tout)0;
						//~ *p_out = (Tout)*p_in1;
					}
				}
			//~ }
			//~ catch (vigra::StdException & e)
			//~ {
				//~ // catch any errors that might have occured and print their reason
				//~ std::cout << e.what() << std::endl;
			//~ }
			break;
		}
	}
	dbi.WriteOut("Leaving ImSubImage with a Mask.\n");
}

/*! @fn template<typename Tin, typename Tin2, typename Tout> void ImSubImage(const BasicImage<Tin> & imIn, const ClipMethod method, const BasicImage<Tin2> & imIn2, BasicImage<Tout> & imOut)
 * @brief Operates the subtraction of imIn2 from imIn1. 
 * @param imIn1 : input image
 * @param method : either Clip or NoClip (with Clip: Truncation is done if result is negative).
 * @param imIn2 : input image (must have the same type than imIn1)
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tin2, typename Tout> void ImSubImage(const DipImage<Tin> & imIn, const ClipMethod method, const DipImage<Tin2> & imIn2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImSubImage\n");
	int w = imIn.width(), h = imIn.height();
	DipImage<unsigned char> imMask(w,h); imMask.init(255);
	ImSubImage(imIn, imMask, method, imIn2, imOut);
	dbi.WriteOut("Leaving ImSubImage.\n");
}


/**
 *
 * @param imIn
 * @param imMask
 * @param firstlabel
 * @param nh
 * @param imOut
 * @return
 */
template<typename Tin, typename Tmask, typename Tout> int ImLabel_WithMask(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const int firstlabel, Neighborhood<int> nh, DipImage<Tout> & imOut)
{
// 	dbi.WriteEnter("Entering ImLabel_WithMask\n");
	const int w = imIn.width(), h = imIn.height();


// 	dbi.WriteInfo("Taille du voisinage = %d\n", nh.size());

	DipImage<uchar> imTmp(w,h); imTmp.init(0);
	uchar * p_tmpstart= &imTmp(0,0);
	uchar * p_tmp = &imTmp(0,0);
	uchar * p_tmpnh = &imTmp(0,0);

	const Tin * p_instart = &imIn(0,0);
	const Tin * p_inx = &imIn(0,0);
	const Tin * p_iny = &imIn(0,0);
	const Tin * p_intmp = &imIn(0,0);
	const Tin * p_inend = p_instart + w*h;
	const Tin * p_in = p_instart;
	const Tin * p_inlinestart = p_instart;
	const Tin * p_inlineend = p_instart + w;

	const Tmask * p_maskstart = &imMask(0,0);
	const Tmask * p_mask = &imMask(0,0);
	const Tmask * p_masktmp = p_mask;

	Tout * p_outstart = &imOut(0,0);
	Tout * p_out = &imOut(0,0);
	// Tout * p_outfin = &imOut(0,0) + w*h;

	int absoffset = 0;

	int currentlabel = firstlabel;

	std::list<int> centers; centers.clear();
	std::list<int> neigh; neigh.clear();
	std::list<int>::iterator itl, itlend;

	Neighborhood<int>::iterator itnh, itnhend = nh.end();

	for (;p_in != p_inend; p_in++, p_tmp++, p_mask++){ // Pour tous les points de l'image
		if ((*p_mask) != 0){ // Si le point appartient au masque
			if ((*p_in) != 0){ // Si le point est non nul
				if ((*p_tmp) == 0){ // Si le point n'est pas marque
					absoffset = p_in - p_instart; // on calcule l'offset
					// if (absoffset == 0){ dbi.WriteInfo("Offset = %d\n", absoffset); }
					centers.push_back(absoffset); // on empile le point
					*p_tmp = 10; // on marque le point
					while(!centers.empty()){ // tant que la file de points de centre n'est pas vide
						itl = centers.begin(); // on recale les iterateurs
						itlend = centers.end();
						for (;itl != itlend; itl++){ // pour tous les points de la liste
							absoffset = *itl;
							// dbi.WriteInfo("\toffset = %d\n", absoffset);
							p_out = p_outstart + absoffset;
							// if (absoffset >= 9){
							//	dbi.WriteInfo("p_out en dehors de l'image - offset = %d - maxoffset = %d - w * h = %d\n", absoffset, p_outfin-p_outstart, w*h);}
							// else
							*p_out = currentlabel;
							p_intmp = p_instart + absoffset;
							//dbi.WriteInfo("offset = %d\n", absoffset);
							p_inlinestart = p_instart + int(w * floor(double(absoffset) / double(w)));
							// dbi.WriteInfo("offset = %d\n", (p_inlinestart-p_instart));
							p_inlineend = p_inlinestart + w;
							itnh = nh.begin();
							for (;itnh != itnhend; itnh++){
								// dbi.WriteInfo("Offset du voisinage = (%d,%d)\n", (*itnh).first, (*itnh).second);
								p_inx = p_intmp + (*itnh).first;
								// dbi.WriteInfo("\tp_inx = %d\n", p_inx - p_instart);
								if ((p_inx >= p_inlinestart) && (p_inx < p_inlineend)){
									p_iny = p_inx + w * (*itnh).second;
									if (p_iny < p_inend){
										absoffset = p_iny - p_instart;
										p_masktmp = p_maskstart + absoffset;
										p_tmpnh = p_tmpstart + absoffset;
										if ((*p_masktmp) != 0){
											if ((*p_iny) != 0){
												if ((*p_tmpnh) == 0){
													neigh.push_back(absoffset);
													// dbi.WriteInfo("\t\toffset neigh = %d\n", absoffsettmp);
													*p_tmpnh = 10;
												}
											}
										}
									}
								}
							}
						}
						centers.clear();
						centers.swap(neigh);
						neigh.clear();
					}
					currentlabel++;
				}
			}
		}
	}
	return(currentlabel-1);
}




/*! @fn template<typename Tin, typename Tmask, typename Tconst, typename Tout> void ImAddConst(const BasicImage<Tin> & imIn, const BasicImage<Tmask> & imMask, const ClipMethod method, const Tconst constante, BasicImage<Tout> & imOut)
 * @brief Operates the addition of a constante to imIn. 
 * @param imIn1 : input image
 * @param imMask : mask image
 * @param method : either Clip or NoClip (with Clip: Truncation is done if result is superior to max value).
 * @param constante : scalar value to be added to the input image
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tmask, typename Tconst, typename Tout> void ImAddConst(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const ClipMethod method, const Tconst constante, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImAddConst with a Mask\n");
	
	int w = imIn.width(), h = imIn.height();
	DipImage<Tconst> imIn2(w,h); imIn2.init(constante);
	ImAddImage(imIn, imMask, method, imIn2, imOut);
	
	dbi.WriteOut("Leaving ImAddConst with a Mask.\n");
	
}

/*! @fn template<typename Tin, typename Tconst, typename Tout> void ImAddConst(const BasicImage<Tin> & imIn, const ClipMethod method, const Tconst constante, BasicImage<Tout> & imOut)
 * @brief Operates the addition of a constante to imIn. 
 * @param imIn1 : input image
 * @param method : either Clip or NoClip (with Clip: Truncation is done if result is superior to max value).
 * @param constante : scalar value to be added to the input image
 * @param [out] imOut : resulting image
 */
template<typename Tin, typename Tconst, typename Tout> void ImAddConst(const DipImage<Tin> & imIn, const ClipMethod method, const Tconst constante, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImAddConst\n");
	int w = imIn.width(), h = imIn.height();
	DipImage<unsigned char> imMask(w,h); imMask.init(255);
	ImAddConst(imIn, imMask, method, constante, imOut);
	dbi.WriteOut("Leaving ImAddConst.\n");
}


/*! 
 * @brief Merges three single channel images into a color image. Each input image is considered as a channel of the output
 * @param imCanal1 : first channel image
 * @param imCanal2 : second channel image
 * @param imCanal3 : third channel image
 * @param imIn : output color image
 */
template<typename T> void ImColorBandMerge(const DipImage<T> & imCanal1, const DipImage<T> & imCanal2, const DipImage<T> & imCanal3, DipColorImage<T> & imOut)
{
  dbi.WriteEnter("Entering ImColorBandMerge...\n");
  const T * p_canal1 = &imCanal1(0,0);
  const T * p_canal2 = &imCanal2(0,0);
  const T * p_canal3 = &imCanal3(0,0);

  T * p_out1 = &imOut(0,0)[0];
  T * p_out1end = p_out1 + 3 * imOut.width() * imOut.height();
  T * p_out2 = p_out1 + 1;
  T * p_out3 = p_out1 + 2;

  for (; p_out1 != p_out1end; p_out1+=3, p_out2+=3, p_out3+=3, p_canal1++, p_canal2++, p_canal3++)
  {
    *p_out1 = *p_canal1;
    *p_out2 = *p_canal2;
    *p_out3 = *p_canal3;
  }
  dbi.WriteOut("Leaving ImColorBandMerge.\n");
}


/** 
 * @brief Computes the mean filter of the input image with different methods
 * @param imIn : input image
 * @param imMask : mask image
 * @param op  :  method selection
 * @param parameter  :  v1,v2: integer radius input. neighborLimits: integer radius input and number of neighborLimits input
 * @param imOut : output image
 */
template<typename Tin, typename Tmask, typename Tout > void ImMeanFilter(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask,meanFilter op, std::vector<ParamT> parameter, DipImage<Tout> & imOut){
	
	switch(op){
		case v1:{
			dbi.WriteEnter("Entering ImMeanFilter v1\n");
			if(parameter.size() !=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter.front().type() != typeid(int)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}			
			ParamT parametertype = parameter.front();
			int rayon = boost::get<int>(parametertype);
	
			clock_t start, stop;
			//~ double time = 0.0;
	
			//~ start = clock();
	
			const int neighborhoodwidth = 2*rayon+1;
			int w = imIn.width(), h = imIn.height();
			int i = 0, j = 0, ii = 0, jj = 0;
			int borneInfx = 0, borneSupx = rayon + 1;
	
			std::vector<double> accuYMemory; accuYMemory.resize(rayon+1, 0.0);
			std::vector<double> compteurYMemory; compteurYMemory.resize(rayon+1, 0.0);
	
			double nbpt= 0.0;
			double nbpttemp= 0.0;
			double sommetemp = 0.0;
			double somme = 0.0;
			double compteurligneY = 0.0;
			double sommeligneY = 0.0;
			
			double valtemp = 0.0;
			
			int index = 0;
			
			dbi.WriteInfo("Valeur du rayon = %d\n", rayon);
			
			//~ stop = clock();
			//~ time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
			//~ std::cout << "Step 1 = " << time << std::endl;
			
			//~ start = clock();
			// Initialisation
			dbi.WriteInfo("Debut de l'initialisation\n");
			for (jj=0; jj<rayon+1; jj++){
				for (ii=0; ii<rayon+1; ii++){
					if (imMask(ii,jj) != 0){
						valtemp = double(imIn(ii,jj));
						accuYMemory[jj] = accuYMemory[jj] + valtemp; // Somme des valeurs horizontalement
						compteurYMemory[jj] = compteurYMemory[jj] + 1.0;
						somme += valtemp;
						nbpt += 1.0;
					}
				}
			}
	
			dbi.WriteInfo("Fin de l'initialisation\n");
			//~ stop = clock();
			//~ time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
			//~ std::cout << "Step 2 = " << time << std::endl;
			
			//~ start = clock();
			std::vector<double> compteurY; compteurY.resize(rayon+1, 0.0);
			std::vector<double> accuY; accuY.resize(rayon+1, 0.0);
			
			double accutime = 0.0;
			// Boucle principale
			for (i=0; i<w; i++){
				// dbi.WriteInfo("Debut de la colonne %d\n",i);
				start = clock();
				accuY = accuYMemory;
				compteurY = compteurYMemory;    
				sommetemp = somme;
				nbpttemp = nbpt;
    
				index = 0;
    
				if (nbpttemp != 0){
					if (imMask(i,0) != 0){
						imOut(i,0) = (Tout)(sommetemp / nbpttemp);
					}
				}
    
				borneInfx = (i-rayon > 0) ? (i-rayon) : 0; // ATTENTION
				borneSupx = (i+rayon+1 < w) ? (i+rayon+1) : w;
    
	
				for (j=1; j<h; j++){
					jj = j+rayon;
					// Si on peut entrer une nouvelle ligne
					if (jj < h){
						for (ii=borneInfx; ii<borneSupx; ii++){
							if (imMask(ii,jj) != 0){
								sommeligneY += double(imIn(ii,jj));
								compteurligneY++;
							}
						}
						sommetemp += sommeligneY;
						nbpttemp += compteurligneY;
						accuY.push_back(sommeligneY);
						compteurY.push_back(compteurligneY);
						sommeligneY = 0.0;
						compteurligneY = 0.0;
					}
	
					// Si on ne touche plus le bord
					if ((j-rayon) > 0){
						sommetemp -= accuY[index];
						nbpttemp -= compteurY[index];
						index++;
					}
	
					if (nbpttemp != 0){
						if (imMask(i,j) != 0){
							imOut(i,j) = (Tout)(sommetemp / nbpttemp);
						}
					}
				}
    
	
				// Reinitialisation
				accuY.resize(rayon+1, 0);
				compteurY.resize(rayon+1, 0);
				index = 0;
    
    
				if (borneSupx >= neighborhoodwidth){
					for (jj=0; jj<(rayon+1); jj++){
						if (imMask(borneInfx,jj) != 0){
							valtemp = double(imIn(borneInfx,jj));
							accuYMemory[jj] -= valtemp;// accuYMemory[jj] - valtemp;
							compteurYMemory[jj]--;// = compteurYMemory[jj] - 1;
							somme -= valtemp;
							nbpt--;
						}
					}
				}
	
				if (borneSupx < w){
					for (jj=0; jj<(rayon+1); jj++){
						if (imMask(borneSupx,jj) != 0){
							valtemp = double(imIn(borneSupx,jj));
							accuYMemory[jj] += valtemp;// accuYMemory[jj] + valtemp;
							compteurYMemory[jj]++;// = compteurYMemory[jj] + 1;
							somme += valtemp;
							nbpt++;
						}
					}
				}
				stop = clock();
				accutime += (double(stop) - double(start))/CLOCKS_PER_SEC;
			}
	
			// stop = clock();
			// time = (double(stop) - double(start)) / CLOCKS_PER_SEC;
			//~ std::cout << "Step 3 = " << accutime / w << std::endl;
	
			dbi.WriteOut("Leaving ImMeanFilter v1\n");
		}
		break;
	
		case v2:{
			if(parameter.size() !=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter.front().type() != typeid(int)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
	
			dbi.WriteEnter("Entering ImMeanFilter V2\n");
	
			ParamT parametertype = parameter.front();
			int rayon = boost::get<int>(parametertype);
			int w = imIn.width();
			int h = imIn.height();

	
			DipImage<double> imMoy(w,h);
			double rapport = 0;
			int i = 0, j = 0, ii = 0, jj = 0;
			double counter = 0;
			double somme = 0;
			int borneInfx = 0, borneMaxx = 0;
			
			
			// DEBUT Nouvelle Boucle 1
			dbi.WriteInfo("Acceleration on\n");
			std::list<int> listecompteur;
			std::list<double> listevaleurs;
			int compteurligne = 0;
			double sommeligne = 0;
			
			for (i=0; i<w; i++){
				j=0;
				somme = 0.0;
				counter = 0;
				listecompteur.clear();
				listevaleurs.clear();
	
				if (i<=rayon){
					borneInfx = 0;
					borneMaxx = i+rayon+1;
				}
				else if ((w-i)<=rayon){
					borneInfx = i-rayon;
					borneMaxx = w;
				}
				else{
					borneInfx = i-rayon;
					borneMaxx = i+rayon+1;
				}
	
				for (jj=0; jj<rayon+1; jj++){
					sommeligne = 0.0;
					compteurligne = 0;
					for (ii= borneInfx; ii< borneMaxx; ii++){
						if (imMask(ii,jj) != 0){
							sommeligne += double(imIn(ii, jj));
							compteurligne += 1;
						}
					}
					somme += sommeligne;
					counter += compteurligne;
					listevaleurs.push_back(sommeligne);
					listecompteur.push_back(compteurligne);
				}
				//std::cout<<"Tu es sorti de la boucle: "<<i<<", "<<j<<std::endl;
				if (imMask(i,j) != 0){
					if (counter != 0){
						rapport = Tout(somme / counter);
						// if (rapport < int(10)) std::cout<<"-----rapport = "<<int(rapport)<<" = "<<i<<" / "<<j<<std::endl;
						// if (rapport > int(255)) std::cout<<"+++++rapport = "<<int(rapport)<<" = "<<i<<" / "<<j<<std::endl;
						imOut(i,j) = Tout(rapport);
						//if (imIn(i,j) > rapport)
						//imOut(i,j) = imIn(i,j);
						//else
						//imOut(i,j) = 0;
					}
					else{
						imOut(i,j) = 0;
						dbi.WriteInfo("Problem : i = %d, j = %d\n", i, j);
					}
					//std::cout<<"Tu es sorti de la boucle 2: "<<somme<<", "<<counter<<std::endl;
				}
	
				for (j=1; j<h; j++){
					if (j<=rayon){
						// Ajout d'une nouvelle ligne
						jj = j+rayon;
						sommeligne = 0.0;
						compteurligne = 0;
						for (ii=borneInfx; ii<borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								sommeligne += double(imIn(ii, jj));
								compteurligne += 1;
							}
						}
						listevaleurs.push_back(sommeligne);
						listecompteur.push_back(compteurligne);
						somme += sommeligne;
						counter += compteurligne;
					}
					else if ((h-j)<=rayon){
						// On retire une ligne
						//std::cout<<"On retire"<<std::endl;
						somme -= listevaleurs.front();
						counter -= listecompteur.front();
						listevaleurs.pop_front();
						listecompteur.pop_front();
					}
					else{
						// On ajoute et on retire une ligne
						jj = j+rayon;
						sommeligne = 0;
						compteurligne = 0;
						for (ii=borneInfx; ii<borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								sommeligne += double(imIn(ii, jj));
								compteurligne += 1;
							}
						}
						listevaleurs.push_back(sommeligne);
						listecompteur.push_back(compteurligne);
						somme = somme + sommeligne - listevaleurs.front();
						counter = counter + compteurligne - listecompteur.front();
						listevaleurs.pop_front();
						listecompteur.pop_front();
					}
					if (imMask(i,j) != 0){
						if (counter != 0){
							rapport = Tout(somme / counter);
							// if (rapport < int(10)) std::cout<<"-----rapport = "<<int(rapport)<<" = "<<i<<" / "<<j<<std::endl;
							// if (rapport > int(255)) std::cout<<"+++++rapport = "<<int(rapport)<<" = "<<i<<" / "<<j<<std::endl;
							imOut(i,j) = Tout(rapport);
						}
						else{
							imOut(i,j) = 0;
							dbi.WriteInfo("Problem : i = %d, j = %d\n", i, j);
						}
					}
				}
			}
	
			dbi.WriteOut("Leaving ImMeanFilter V2\n");
		}

		break;
		case neighborLimits:{
			if(parameter.size() !=2){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(int) && parameter[1].type() != typeid(int)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
	
			dbi.WriteEnter("Entering ImMeanFilter neighborLimits\n");
	
			ParamT parametertype = parameter.front();
			int rayon = boost::get<int>(parametertype);
			parametertype = parameter[1];
			int neighborlimit = boost::get<int>(parametertype);
	
			int w = imIn.width(), h = imIn.height();
			int i = 0, j = 0, ii = 0, jj = 0;
			int borneInfx = 0, borneMaxx = 0;
			
			std::list<double> listecompteur;
			std::list<double> listevaleurs;
			//   double rapport = 0.0;
			double counter = 0.0;
			double somme = 0.0;
			double compteurligne = 0.0;
			double sommeligne = 0.0;
			
			//   exportImage(srcImageRange(imIn), ImageExportInfo("/data/images/test/imagesoriginales/softex/toto.png"));
			//   exportImage(srcImageRange(imMask), ImageExportInfo("/data/images/test/imagesoriginales/softex/totomask.png"));
			dbi.WriteInfo("Valeur du rayon = %d\n", rayon);
			
			// DEBUT Nouvelle Boucle 1
			for (i=0; i<w; i++){
				j=0;
				somme = 0.0;
				counter = 0.0;
				listecompteur.clear();
				listevaleurs.clear();
				if (i<=rayon){
					borneInfx = 0;
					borneMaxx = i+rayon+1;
				}
				else if ((w-i)<=rayon){
					borneInfx = i-rayon;
					borneMaxx = w;
				}
				else{
					borneInfx = i-rayon;
					borneMaxx = i+rayon+1;
				}

				for (jj=0; jj<rayon+1; jj++){
					sommeligne = 0.0;
					compteurligne = 0.0;
					for (ii= borneInfx; ii< borneMaxx; ii++){
						if (imMask(ii,jj) != 0){
							// if (i==50 && j==0)
							// dbi.WriteInfo("Valeur de l'image = %g\n", double(imIn(ii, jj)));
							sommeligne += double(imIn(ii, jj));
							compteurligne += 1.0;
						}
					}
					somme += sommeligne;
					counter += compteurligne;
					listevaleurs.push_back(sommeligne);
					listecompteur.push_back(compteurligne);
				}
				if (imMask(i,j) != 0){
					if (counter > (double)neighborlimit){
						// dbi.WriteInfo("(%d,%d) : somme = %g, counter = %g, rapport = %g\n", i, j, somme, counter, somme/counter);
						imOut(i,j) = Tout(somme / counter);
					}
					else{
						imOut(i,j) = 0;
						// dbi.WriteInfo("Problem : i = %d, j = %d\n", i, j);
					}
				}

				for (j=1; j<h; j++){
					if (j<=rayon){
						// Ajout d'une nouvelle ligne
						jj = j+rayon;
						sommeligne = 0.0;
						compteurligne = 0.0;
						for (ii=borneInfx; ii<borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								sommeligne += double(imIn(ii, jj));
								compteurligne += 1.0;
							}
						}
						listevaleurs.push_back(sommeligne);
						listecompteur.push_back(compteurligne);
						somme += sommeligne;
						counter += compteurligne;
					}
	
					else if ((h-j)<=rayon){
						// On retire une ligne
						somme -= listevaleurs.front();
						counter -= listecompteur.front();
						listevaleurs.pop_front();
						listecompteur.pop_front();
					}
	
					else{
						// On ajoute et on retire une ligne
						jj = j+rayon;
						sommeligne = 0.0;
						compteurligne = 0.0;
						for (ii=borneInfx; ii<borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								sommeligne += double(imIn(ii, jj));
								compteurligne += 1.0;
							}
						}
						listevaleurs.push_back(sommeligne);
						listecompteur.push_back(compteurligne);
						somme = somme + sommeligne - listevaleurs.front();
						counter = counter + compteurligne - listecompteur.front();
						listevaleurs.pop_front();
						listecompteur.pop_front();
					}
	
					if (imMask(i,j) != 0){
						if (counter > (double)neighborlimit){
							imOut(i,j) = Tout(somme / counter);
						}
						else{
							imOut(i,j) = 0;
							// dbi.WriteInfo("Problem : i = %d, j = %d\n", i, j);
						}
					}
				}
			}
			dbi.WriteOut("Leaving ImMeanFilter neighborLimits\n");
		}
		break;
	}
}
/*! 
 * @brief stretching histrgram  
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : clipping, valsetting, regular
 * @param parameter : input parameter
 * @param imOut : output image
 * @todo when new method input, modify brief
 */
template<typename Tin, typename Tmask, typename Tout> void ImStretchHistogram(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, imStretchHistogramMethod method, std::vector<ParamT> parameter, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImStretchHistogram withMask\n");
	const int w = imIn.width(), h = imIn.height();
	
	switch(method){
		case clipping:{
			if(parameter.size() !=4){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(Tin) && parameter[1].type() != typeid(Tin) && parameter[2].type() != typeid(Tout) && parameter[3].type() != typeid(Tout) ){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}			
			ParamT parametertype = parameter[0];
			Tin val_min_in = boost::get<Tin>(parametertype);
			parametertype = parameter[1];
			Tin val_max_in = boost::get<Tin>(parametertype);
			parametertype = parameter[2];
			Tout val_min_out = boost::get<Tout>(parametertype);
			parametertype = parameter[3];
			Tout val_max_out = boost::get<Tout>(parametertype);
			
			std::pair<Tin, Tin> minmax = MeasMinMax(imIn, imMask);
			dbi.WriteInfo("---------> input min = %f\n", float(minmax.first));
			dbi.WriteInfo("---------> input max = %f\n", float(minmax.second));
			Tin minin = val_min_in; // Tin(std::max(int(std::numeric_limits<Tin>::min()),val_min_in));
			Tin maxin = val_max_in; // Tin(std::min(int(std::numeric_limits<Tin>::max()),val_max_in));
			dbi.WriteInfo("---------> val_min_in = %f\n", float(val_min_in));
			dbi.WriteInfo("---------> val_min_out = %f\n", float(val_min_out));
			Tin rangein = maxin - minin;
			dbi.WriteInfo("Range d'entree = %g ( %g - %g)\n", double(rangein), double(maxin), double(minin));
			Tout rangeout = val_max_out - val_min_out;
			dbi.WriteInfo("Range de sortie = %g ( %g - %g)\n", double(rangeout), double(val_max_out), double(val_min_out));

			DipImage<uchar> imTest(w,h); imTest.init(0);

			const Tin * p_in = &imIn(0,0);
			const Tin * p_inend = p_in + w*h;
			const Tmask * p_mask = &imMask(0,0);
			Tout * p_out = &imOut(0,0);
			uchar * p_test = &imTest(0,0);

			for (;p_in != p_inend; p_in++, p_out++, p_mask++, p_test++){
				if (*p_mask != 0){
					if (*p_in <= val_min_in){
						*p_out = val_min_out;
						*p_test = (uchar)100;
					}
					else if (*p_in >= val_max_in){
						*p_out = val_max_out;
						*p_test = (uchar)150;
					}
					else{
						*p_out = static_cast<Tout>(double(val_min_out) + double(*p_in - val_min_in) * double(rangeout) / double(rangein));
						*p_test = (uchar)255;
					}
				}
				else{
					*p_out = (Tout)(0);
				}
			}
			dbi.WriteImage(imTest, "Truc_imTest.png");
			minmax = MeasMinMax(imOut, imMask);
			dbi.WriteInfo("---------> output min = %d\n", minmax.first);
			dbi.WriteInfo("---------> output max = %d\n", minmax.second);			
		}
		break;
		
		case valsetting:{
			if(parameter.size() !=2){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			}
			if(parameter[0].type() != typeid(Tout) && parameter[1].type() != typeid(Tout)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}			
			ParamT parametertype = parameter[0];
			Tout val_min = boost::get<Tout>(parametertype);
			parametertype = parameter[1];
			Tout val_max = boost::get<Tout>(parametertype);
			
			std::pair<Tin, Tin> minmax = MeasMinMax(imIn, imMask);
			double rangein = double(minmax.second - minmax.first);
			Tin minmaxfirst = minmax.first;
			dbi.WriteInfo("Range d'entree = %g ( %g - %g)\n", double(rangein), double(minmax.second), double(minmax.first));
			Tout rangeout = val_max - val_min;
			dbi.WriteInfo("Range de sortie = %g ( %g - %g)\n", double(rangeout), double(val_max), double(val_min));

			const Tin * p_in = &imIn(0,0);
			const Tin * p_inend = p_in + imIn.width() * imIn.height();
			const Tmask * p_mask = &imMask(0,0);
			Tout * p_out = &imOut(0,0);

			for (; p_in != p_inend; p_in++, p_mask++, p_out++){
				if (*p_mask != 0)
					*p_out = val_min + (Tout)((*p_in - minmaxfirst) * double(rangeout) / double(rangein));
				else
					*p_out = 0;
			}
		}
		break;
		case regular:{
			//~ if(parameter.size() !=0){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameter is not correct");
			//~ }
			Tout val_min,val_max;
			
			if (std::numeric_limits<Tout>::min() != 0){
				val_min = -std::numeric_limits<Tout>::max();
			}
			else{
				val_min = std::numeric_limits<Tout>::min();
			}
			
			val_max = std::numeric_limits<Tout>::max();
			std::pair<Tin, Tin> minmax = MeasMinMax(imIn, imMask);
			double rangein = double(minmax.second - minmax.first);
			Tin minmaxfirst = minmax.first;
			dbi.WriteInfo("Range d'entree = %g ( %g - %g)\n", double(rangein), double(minmax.second), double(minmax.first));
			Tout rangeout = val_max - val_min;
			dbi.WriteInfo("Range de sortie = %g ( %g - %g)\n", double(rangeout), double(val_max), double(val_min));

			const Tin * p_in = &imIn(0,0);
			const Tin * p_inend = p_in + imIn.width() * imIn.height();
			const Tmask * p_mask = &imMask(0,0);
			Tout * p_out = &imOut(0,0);

			for (; p_in != p_inend; p_in++, p_mask++, p_out++){
				if (*p_mask != 0)
					*p_out = val_min + (Tout)((*p_in - minmaxfirst) * double(rangeout) / double(rangein));
				else
					*p_out = 0;
			}
			//~ cout << "REGULAR OUT" << endl;
		}
		break;
	}
	
	dbi.WriteOut("Leaving ImStretchHistogram withMask\n");
}

/**
 * 
 * @brief Super contrast enhancement filter(use stretch histogram by setting custom low/high value)
 * @param imIn
 * @param imMask
 * @param rayon
 * @param imLight
 * @param imEnhanced
 */

template<typename Tin, typename Tmask, typename Tlight, typename Tenhanced> void ImSuperContrastEnhancementFilter_WithMask(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, int rayon, DipImage<Tlight> & imLight, DipImage<Tenhanced> & imEnhanced)
{
	dbi.WriteEnter("Entering ImSCE_WithMask\n");
	int w=imIn.width(), h=imIn.height();
	std::vector<ParamT> parameterVector;parameterVector.clear();
	
	DipImage<Tin> imMean(w,h); imMean.init(0);
	DipImage<double> imSub(w,h); imSub.init(0);
	parameterVector.push_back(rayon);
	ImMeanFilter(imIn, imMask, v1,parameterVector,imMean);
	parameterVector.clear();	
	ImSubImage(imIn, imMask , NoClip, imMean , imSub);
	
	std::vector<double> meanstd = MeasStatistics(imSub, imMask, all);
	dbi.WriteInfo("meanSCE = %f\n", float(meanstd[0]));
	dbi.WriteInfo("stdSCE = %f\n", float(meanstd[1]));
	
	parameterVector.push_back((double)(meanstd[0] - 3*meanstd[1]));
	parameterVector.push_back((double)(meanstd[0] + 3*meanstd[1]));
	parameterVector.push_back((Tenhanced)0);
	parameterVector.push_back((Tenhanced)255);
	ImStretchHistogram(imSub, imMask, clipping, parameterVector, imEnhanced);
		
	meanstd.clear();
	meanstd = MeasStatistics(imIn, imMask, all);
	ImAddConst(imSub, imMask, Clip, meanstd[0], imLight);
		
	dbi.WriteOut("Leaving ImSCE_WithMask\n");
}


/**
 * @fn template<typename Tin, typename Tmask, typename Tenhanced> void ImColorSuperContrastEnhancementFilter_WithMask_MT(const vigra::BasicImage<RGBValue<Tin> > & imIn, const vigra::BasicImage<Tmask> & imMask, const int rayon, vigra::BasicImage<RGBValue<Tenhanced> >  & imEnhanced)
 * @brief use threads to accelerate color super contrast enhancement filter
 * @param imIn :input color image
 * @param imMask : mask image
 * @param rayon : radius for mean filter
 * @param imEnhanced : output color image
 */
template<typename Tin, typename Tmask, typename Tlight, typename Tenhanced> void ImEnhancementLighting(const DipColorImage<Tin> & imIn, const DipImage<Tmask> & imMask, const int rayon, DipColorImage<Tlight> & imLight, DipColorImage<Tenhanced> & imEnhanced)
{
	dbi.WriteEnter("Entering ImColorSuperContrastEnhancementFilter_WithMask_MT\n");	
	int w=imIn.width(), h=imIn.height();
	
	DipImage<Tin> inputimages1(w,h);
	DipImage<Tin> inputimages2(w,h);
	DipImage<Tin> inputimages3(w,h);
	
	ImColorBandSeparation(imIn, inputimages1, inputimages2, inputimages3);	
	
	DipImage<Tlight> lightimages1(w,h);
	DipImage<Tlight> lightimages2(w,h);
	DipImage<Tlight> lightimages3(w,h);
	
	DipImage<Tenhanced> outputimages1(w,h);
	DipImage<Tenhanced> outputimages2(w,h);
	DipImage<Tenhanced> outputimages3(w,h);	
	
	ImSuperContrastEnhancementFilter_WithMask(inputimages1, imMask, rayon,lightimages1,outputimages1);
	ImSuperContrastEnhancementFilter_WithMask(inputimages2, imMask, rayon,lightimages2,outputimages2);
	ImSuperContrastEnhancementFilter_WithMask(inputimages3, imMask, rayon,lightimages3,outputimages3);
	
	//~ boost::thread_group tg;

	//~ tg.create_thread(boost::bind(&ImSuperContrastEnhancementFilter_WithMask<Tin,Tmask,Tlight, Tenhanced>,inputimages1, imMask, rayon, boost::ref(lightimages1), boost::ref(outputimages1)) );
	//~ tg.create_thread(boost::bind(&ImSuperContrastEnhancementFilter_WithMask<Tin,Tmask,Tlight, Tenhanced>,inputimages2, imMask, rayon, boost::ref(lightimages2), boost::ref(outputimages2)) );
	//~ tg.create_thread(boost::bind(&ImSuperContrastEnhancementFilter_WithMask<Tin,Tmask,Tlight, Tenhanced>,inputimages3, imMask, rayon, boost::ref(lightimages3), boost::ref(outputimages3)) );	
	//~ tg.join_all();

	ImColorBandMerge(outputimages1, outputimages2, outputimages3, imEnhanced);
	//~ dbi.WriteImage(imEnhanced, "InSide_imEnhanced.png");
	ImColorBandMerge(lightimages1, lightimages2, lightimages3, imLight);
	//~ dbi.WriteImage(imLight, "InSide_imLight.png");
	//~ cin.ignore(1);

	dbi.WriteOut("Leaving ImColorSuperContrastEnhancementFilter_WithMask_MT\n");
}




/** 
 * @brief Computes the mean filter on each channel of the input image
 * @param imIn : input color image
 * @param imMask : mask image
 * @param rayon : neighborhood ray
 * @param imOut : output color image
 */
template<typename Tin, typename Tmask, typename Tout> void ImColorMeanFilter(const DipColorImage<Tin> & imIn, const DipImage<Tmask> & imMask, const int rayon, DipColorImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImColorMeanFilter_WithMask\n");

	int w = imIn.width(), h = imIn.height();

	DipImage<Tin> imRed(w,h); imRed.init(0);
	DipImage<Tin> imGreen(w,h); imGreen.init(0);
	DipImage<Tin> imBlue(w,h); imBlue.init(0);

	DipImage<Tout> imRedOut(w,h); imRedOut.init(0);
	DipImage<Tout> imGreenOut(w,h); imGreenOut.init(0);
	DipImage<Tout> imBlueOut(w,h); imBlueOut.init(0);

	ImColorBandSeparation(imIn, imRed, imGreen, imBlue);
	
	std::vector<ParamT> param;param.clear();
	param.push_back(rayon);

	ImMeanFilter(imRed, imMask, v1, param, imRedOut);
	ImMeanFilter(imGreen, imMask, v1, param, imGreenOut);
	ImMeanFilter(imBlue, imMask, v1, param, imBlueOut);

	ImColorBandMerge(imRedOut, imGreenOut, imBlueOut, imOut);

	dbi.WriteOut("Leaving ImColorMeanFilter_WithMask\n");
}


/*! 
 * 
 * @param imIn Original color image
 * @param imMask Image of ROI
 * @param radius radius of the mean filter
 * @param [out] imOut Reflection mask image
 * @todo independency
 */
template <typename Tin, typename Tmask, typename Tout>
void ImLightingCorrection(
	  const DipColorImage<Tin> & imIn,
	  const DipImage<Tmask> & imMask,
	  const int radius,
	  DipColorImage <Tout> & imOut
	  )
{
	dbi.WriteEnter("Entering ImLightingCorrection\n");
	const int w = imIn.width(), h = imIn.height();
	DipColorImage <uchar> imMean(w,h); imMean.init(0);
	ImColorMeanFilter(imIn, imMask, radius, imMean);

	DipImage<Tin> imRed(w,h); imRed.init(0);
	DipImage<Tin> imGreen(w,h); imGreen.init(0);
	DipImage<Tin> imBlue(w,h); imBlue.init(0);
	ImColorBandSeparation(imIn, imRed, imGreen, imBlue);

	
	std::vector<double> msR = MeasStatistics(imRed, imMask, all);
	std::vector<double> msG = MeasStatistics(imGreen, imMask, all);
	std::vector<double> msB = MeasStatistics(imBlue, imMask, all);

	const DipColor<Tin> * p_in = &imIn(0,0);
	const DipColor<Tin> * p_inend = p_in + w*h;
	DipColor<Tout> * p_out = &imOut(0,0);
	DipColor<Tout> * p_mean = &imMean(0,0);

	const Tmask * p_mask = &imMask(0,0);
	
	for (; p_in != p_inend; p_in++, p_out++, p_mask++, p_mean++)
	{
		if (*p_mask != 0){
			if ((*p_in)[0] + (msR[0] - (*p_mean)[0])<0)
				(*p_out)[0] = 0;
			else if ((*p_in)[0] + (msR[0] - (*p_mean)[0])>std::numeric_limits<Tout>::max())
					(*p_out)[0] = std::numeric_limits<Tout>::max();
			else
				(*p_out)[0] = (*p_in)[0] + (msR[0] - (*p_mean)[0]);
				
			if ((*p_in)[1] + (msG[0] - (*p_mean)[1])<0)
				(*p_out)[1] = 0;
			else if ((*p_in)[1] + (msG[0] - (*p_mean)[1])>std::numeric_limits<Tout>::max())
					(*p_out)[1] = std::numeric_limits<Tout>::max();
			else	
				(*p_out)[1] = (*p_in)[1] + (msG[0] - (*p_mean)[1]);
				
			if ((*p_in)[2] + (msB[0] - (*p_mean)[2])<0)
				(*p_out)[2] = 0;
			else if ((*p_in)[2] + (msB[0] - (*p_mean)[2])>std::numeric_limits<Tout>::max())
					(*p_out)[2] = std::numeric_limits<Tout>::max();
			else	
				(*p_out)[2] = (*p_in)[2] + (msB[0] - (*p_mean)[2]);
		}
	}

	dbi.WriteOut("Leaving ImLightingCorrection\n");
}



template<typename Tmask>void MyFilledEllipse( DipImage<Tmask> & imFaceMask, Rect & roi){
	int thickness = -1;
	int lineType = 8;
	Point center( roi.x + roi.width*0.5, roi.y + roi.height*0.5 );
	ellipse( imFaceMask, center, Size( roi.width*0.5, roi.height*0.5), 0, 0, 360, Scalar( 255, 0, 255 ), thickness, lineType, 0 );
}

template<typename Tmask>void MyFilledRectangle( DipImage<Tmask> & imFaceMask, Rect & roi){
	int thickness = -1;
	int lineType = 8;
	Point pt1( roi.x , roi.y);
	Point pt2( roi.x + roi.width  , roi.y + roi.height);
	rectangle( imFaceMask, pt1,pt2, Scalar( 255, 0, 0 ), thickness, lineType, 0 );
}

template<typename Tin, typename Tmask> cv::Rect detectMouth(DipColorImage<Tin> & frame,DipImage<Tmask> & imFaceMask, std::string path_model)
{
	std::string mouth_cascade_model = path_model + "haarcascade_mcs_mouth.xml";
	CascadeClassifier mouth_cascade;
	if (!mouth_cascade.load(mouth_cascade_model))
	{
		cerr << "Error loading cascade file\n";
		exit(EXIT_FAILURE);
	};
	
	std::vector<Rect> mouth;
	
}
template<typename Tin, typename Tmask> cv::Rect detectFaceAndBody(DipColorImage<Tin> & frame,DipImage<Tmask> & imFaceMask, std::string path_model)
{
	std::string face_cascade_model = path_model + "haarcascade_frontalface_default.xml";
	CascadeClassifier face_cascade;
	if (!face_cascade.load(face_cascade_model))
	{
		cerr << "Error loading cascade file\n";
		exit(EXIT_FAILURE);
	};
	// ---------------------------------------

	std::vector<Rect> faces;
	const int w = frame.width(), h = frame.height();
	DipImage<Tin> frame_gray(w,h); frame_gray.init(0); //grayscale image
	
	// RGB TO GRAYSCALE IMAGE
	cv::cvtColor(frame, frame_gray, COLOR_BGR2GRAY);
	// HISTGRAM EQUALIZATION
	equalizeHist(frame_gray, frame_gray);
	dbi.WriteImage(frame_gray, "gray.png");
	
	// DETECT FACE USING CASCADE CLASSIFIER
	face_cascade.detectMultiScale(frame_gray, faces, 1.1, 2, 0 | CASCADE_SCALE_IMAGE, Size(h/10, w/10));

	
	//SET ROI
	cv::Rect roi_b; //Rect for biggest ROI
	cv::Rect roi_c; // temp Rect for current ROI
	cv::Rect roi_body; // Rect which contains body
	cv::Rect roi_head; // Rect which contains head

	size_t ic = 0; // ic is index of current ROI
	int ac = 0; // ac is area of current ROI

	size_t ib = 0; // ib is index of biggest ROI
	int ab = 0; // ab is area of biggest ROI
	std::cout << faces.size() << endl;
	//~ cin.ignore(1); //pause for debugging

	for (ic = 0; ic < faces.size(); ic++) // Iterate through all current ROIs (detected faces)
	{
		roi_c.x = faces[ic].x;
		roi_c.y = faces[ic].y;
		roi_c.width = (faces[ic].width);
		roi_c.height = (faces[ic].height);

		ac = roi_c.width * roi_c.height; // Get the area of current ROI (detected face)

		roi_b.x = faces[ib].x;
		roi_b.y = faces[ib].y;
		roi_b.width = (faces[ib].width);
		roi_b.height = (faces[ib].height);

		ab = roi_b.width * roi_b.height; // Get the area of biggest ROI, at beginning it is same as "current" ROI

		if (ac > ab)
		{
			ib = ic;
			roi_b.x = faces[ib].x;
			roi_b.y = faces[ib].y;
			roi_b.width = (faces[ib].width);
			roi_b.height = (faces[ib].height);
		}
	}
	roi_head.x = roi_b.x;
	roi_head.y = std::max(roi_b.y - (roi_b.y / 3),0); 
	roi_head.width = roi_b.width;
	roi_head.height = roi_b.height + (roi_b.y / 3);
	MyFilledRectangle(imFaceMask,roi_head);
	
	

	roi_body.x = (roi_b.x - roi_b.width < 0)? 0 : roi_b.x - roi_b.width;
	roi_body.y = roi_b.y + roi_b.height;
	roi_body.height = h - roi_body.y - 1;
	roi_body.width = ((roi_b.x + 0.5 * roi_b.width) - roi_body.x) * 2;
	roi_body.width = (roi_body.x + roi_body.width >= w) ? (w - roi_body.x - 1) : roi_body.width;
	
	MyFilledEllipse(imFaceMask,roi_body);
	
	roi_body.x = (roi_b.x - roi_b.width < 0)? 0 : roi_b.x - roi_b.width;
	roi_body.y = roi_b.y + roi_b.height + (h - roi_b.y - roi_b.height)/2;
	roi_body.height = h - roi_body.y - 1;
	roi_body.width = ((roi_b.x + 0.5 * roi_b.width) - roi_body.x) * 2;
	roi_body.width = (roi_body.x + roi_body.width >= w) ? (w - roi_body.x - 1) : roi_body.width;
	MyFilledRectangle(imFaceMask,roi_body);
	
	
	return roi_head;



}


/*! 
 * @brief Computes the skeleton of a binary image
 * @param imIn : input binary image
 * @param imMask : mask image
 * @param[out] imOut : image with the skeleton
 */
template<typename Tin, typename Tmask, typename Tout> void ImThinning_WithMask(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImThinning\n");

	const int w = imIn.width(), h = imIn.height();
	DipImage<Tin> imTmp(w,h);
	ImCopy(imIn, imTmp);
	dbi.WriteImage(imIn, "imSrcPb.png");
	dbi.WriteImage(imTmp, "imTmpPb.png");

	Tin * p_out = &imTmp(0,0);
	Tin * p_outtmp = &imTmp(0,0);
	Tin * p_outstart = &imTmp(0,0);
	Tin * p_outend = p_out + w * h;

	const Tmask * p_maskstart = &imMask(0,0);
	const Tmask * p_masktmp = &imMask(0,0);
	const Tmask * p_mask = &imMask(0,0);


	std::pair<int,int> minmax = MeasMinMax(imIn, imMask);
	dbi.WriteInfo("minmax = (%d,%d)\n", minmax.first, minmax.second);
	const int taillefah = 2 * minmax.second; //pow(2.0, (double)(8.0 * sizeof(Tin)));
	//dbi.WriteInfo("taille de Tin = %d\n", sizeof(Tin));
	dbi.WriteInfo("taille de la file d'attente = %d\n", int(taillefah + 1));

	std::vector<std::list<int> > fah;
	fah.clear();
	std::list<int> emptylist; emptylist.clear();
	fah.resize(taillefah, emptylist);
	
	DipImage<Tin> imNH(3,3); imNH.init(0);
	DipImage<Tin> imNHStar(3,3); imNHStar.init(0);
	DipImage<Tin> imNHInv(3,3); imNHInv.init(0);
	DipImage<Tin> imNHInvStar(3,3); imNHInvStar.init(0);
	DipImage<uchar> imMaskNH(3,3); imMaskNH.init(255);
	DipImage<int> imLabNH(3,3); imLabNH.init(0);

	Neighborhood<int> nh4, nh8; nh4.clear(); nh8.clear();
	int radius = 1;
	std::vector<ParamT> param; param.clear();
	std::vector<ParamT> paramnb8; paramnb8.clear(); //delete
	param.push_back(radius);
	CreateNeighborhood(Cross, param, nh4);
	CreateNeighborhood(Square, param, nh8);
	
	int nblb = 0;
	int nblbinv = 0;
	int nblbstar = 0;
	int nblbinvstar = 0;

	DipImage<uchar> imFlag(w,h); imFlag.init(0);
	uchar * p_flagstart = &imFlag(0,0);
	uchar * p_flagtmp = &imFlag(0,0);
	uchar * p_flag = &imFlag(0,0);

	int offset = 0;
	std::vector<int> reloffset; reloffset.clear();
	reloffset.push_back(-w-1);
	reloffset.push_back(-w);
	reloffset.push_back(-w+1);
	reloffset.push_back(-1);
	reloffset.push_back(1);
	reloffset.push_back(w-1);
	reloffset.push_back(w);
	reloffset.push_back(w+1);
	std::vector<int>::iterator itro = reloffset.begin(), itroend = reloffset.end();
	
	int ro = 0;
	// uint vol = 0;

	param.clear();
	param.push_back(1);
	param.push_back(nh4);
	
	std::vector<ParamT> param8; param8.clear();
	param8.push_back(1);
	param8.push_back(nh8);
	
	std::vector<ParamT> param4; param4.clear();
	param4.push_back(1);
	param4.push_back(nh4);
				

	// Debut de l'initialisation
	dbi.WriteInfo("Debut de l'initialisation\n");
	for (; p_out != p_outend; p_out++, p_mask++, p_flag++){ // Pour tous les points de l'image
		if (*(p_mask) != 0){ // Si le point est dans le masque
			if ((*p_out) != 0){ // Si le point n'est pas nul ( = appartient a l'objet)
				// Calcul des elements de determination de la simplicite
				offset = p_out - p_outstart;

				CopyNeighborhoodInImage(imTmp, p_out, imNHStar);
				ImBinarisation(imNHStar, imNH);
				// ImCopy(imNH, imNHStar);
				ImInvert(imNH, imNHInv);
				nblb = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
				nblbinv = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
				//nblb = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
				//nblbinv = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
				imNH(1,1) = 0;
				ImInvert(imNH, imNHInv);
				nblbstar = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
				nblbinvstar = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
			
				//nblbstar = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
				//nblbinvstar = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
				// fin de calcul des elements de la simplicite

				if ((nblb == nblbstar) && (nblbinv == nblbinvstar)){ // Si c'est un point simple
					if (*p_flag == 0){ // Si le point n'a pas deja ete traite
						// on l'empile et on flag le point
						*p_flag = 10;
						// if (*p_out > minmax.second)
						//	dbi.WriteInfo("p_out = %d\n", *p_out);
						fah[*p_out].push_back(offset);
					}
				}
			}
		}
	}
	// dbi.WriteImage(imFlag, "imFlagFinInit.png");
	// Fin de l'initialisation
	dbi.WriteInfo("Fin de l'initialisation\n");
	// int nn = 0;
	// for (; nn < fah.size(); nn++){
	// 	std::cout << nn << " -> " << fah[nn].size() << std::endl;
	// }

	int i = 0;
	for (i=0; i<taillefah; i++){ // Pour tous les niveaux de gris de l'image
		// dbi.WriteInfo("Taille de la file %d = %d\n", i, fah[i].size());
		while ( !fah[i].empty() ){ // Tant que la file du niveau de gris courant n'est pas vide
			offset = fah[i].front(); // on recupere l'offset
			// dbi.WriteInfo("offset = %d\n", offset);
			fah[i].pop_front(); // on depile le point
			p_out = p_outstart + offset; // On positionne le pointeur de l'image initiale

			CopyNeighborhoodInImage(imTmp, p_out, imNHStar);
			ImBinarisation(imNHStar, imNH);
			ImInvert(imNH, imNHInv);
			nblb = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
			nblbinv = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
						
			//nblb = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
			//nblbinv = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
			imNH(1,1) = 0;
			ImInvert(imNH, imNHInv);
			
			nblbstar = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
			nblbinvstar = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
			
			//nblbstar = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
			//nblbinvstar = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
			// vol = ImMeasVolume(imNH);
			ImMeasVolume(imNH);
			if ((nblb == nblbstar) && (nblbinv == nblbinvstar)){ // Si c'est toujours un point simple
				if (255 < ImMeasVolume(imNH)){
					*p_out = 0;
					// on empile les voisins qui sont dans le masque, simple et non deja empiles
					itro = reloffset.begin();
					for (;itro != itroend; itro++){
						p_masktmp = p_maskstart + offset + *itro;
						if (*p_masktmp != 0){
							p_flagtmp = p_flagstart + offset + *itro;
							if (*p_flagtmp == 0){
								p_outtmp = p_out + *itro;
								if (*p_outtmp != 0){
									CopyNeighborhoodInImage(imTmp, p_outtmp, imNHStar);
									ImBinarisation(imNHStar, imNH);
									ImInvert(imNH, imNHInv);
									
									nblb = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
									nblbinv = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
			
									
									//nblb = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
									//nblbinv = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
									imNH(1,1) = 0;
									ImInvert(imNH, imNHInv);
									
									nblbstar = ImLabel(imNH, imMaskNH, IterationWithNeighbour, param8, imLabNH);
									nblbinvstar = ImLabel(imNHInv, imMaskNH, IterationWithNeighbour, param4, imLabNH);
			
									//nblbstar = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
									//nblbinvstar = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
									if ((nblb == nblbstar) && (nblbinv == nblbinvstar)){
										ro = p_outtmp - p_outstart;
										fah[*p_outtmp].push_back(ro);
										*p_flagtmp = 10;
										if (*p_outtmp < i)
											i = *p_outtmp;
									}
								}
							}
						}
					}
				}
				else{
					p_flag = p_flagstart + offset;
					*p_flag = 10;
				}

			}
			else{
				p_flag = p_flagstart + offset;
				*p_flag = 0;
			}
		}
	}

	ImBinarisation(imTmp, imOut);
	dbi.WriteImage(imTmp, "imFin.png");
	dbi.WriteOut("Leaving ImThinning\n");
}

/**
 *
 * @param imIn
 * @param imMask
 * @param imOut
 */
template<typename Tin, typename Tmask, typename Tout>
void ImPointsClassification_WithMask(
		const DipImage<Tin> & imIn,
		const DipImage<Tmask> & imMask,
		DipImage<Tout> & imOut
		)
{
	dbi.WriteEnter("Entering ImPointsClassification_WithMask\n");

	const int w = imIn.width(), h = imIn.height();
	DipImage<Tin> imTmp(w,h);
	ImCopy(imIn, imTmp);

	Tin * p_outstart = &imTmp(0,0);
	Tin * p_out = &imTmp(0,0);
	// Tin * p_outtmp = &imTmp(0,0);
	Tin * p_outend = p_out + w * h;

	const Tmask * p_mask = &imMask(0,0);

	DipImage<Tin> imNH(3,3); imNH.init(0);
	DipImage<Tin> imNHStar(3,3); imNHStar.init(0);
	DipImage<Tin> imNHInv(3,3); imNHInv.init(0);
	DipImage<Tin> imNHInvStar(3,3); imNHInvStar.init(0);
	DipImage<uchar> imMaskNH(3,3); imMaskNH.init(255);
	DipImage<TLabel> imLabNH(3,3); imLabNH.init(0);

	Neighborhood<int> nh4, nh8; nh4.clear(); nh8.clear();
	const int radius = 1;
	SquareCross2DNListGenerator<int,int>(radius,nh4);
	Square2DNListGenerator<int,int>(radius, nh8);


	// int nblb = 0;
	int nblbinv = 0;
	int nblbstar = 0;
	// int nblbinvstar = 0;

	DipImage<uchar> imFlag(w,h); imFlag.init(0);
	uchar * p_flag = &imFlag(0,0);

	int offsetref = 132 + 56 *200;
	int offset = 0;

	// Debut de l'initialisation
	dbi.WriteInfo("Debut de l'initialisation\n");
	for (; p_out != p_outend; p_out++, p_mask++, p_flag++){ // Pour tous les points de l'image
		if (*(p_mask) != 0){ // Si le point est dans le masque
			if ((*p_out) != 0){ // Si le point n'est pas nul ( = appartient a l'objet)
				// Calcul des elements de determination de la simplicite
				offset = p_out - p_outstart;

				CopyNeighborhoodInImage(imTmp, p_out, imNHStar);

				ImBinarisation(imNHStar, imNH);
				if (offset == offsetref){
					dbi.WriteImage(imNH, "imNeighborhood_NH.png");
				}
				// ImCopy(imNH, imNHStar);
				ImInvert(imNH, imNHInv);
				if (offset == offsetref){
					dbi.WriteImage(imNHInv, "imNeighborhood_NHInv.png");
				}
				// nblb = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);
				nblbinv = ImLabel_WithMask(imNHInv, imMaskNH, 1, nh4, imLabNH);
				imNH(1,1) = 0;
				if (offset == offsetref){
					dbi.WriteImage(imNH, "imNeighborhood_NHStar.png");
				}
				nblbstar = ImLabel_WithMask(imNH, imMaskNH, 1, nh8, imLabNH);


				switch (nblbinv){
					case 0:// point interieur
						if (p_out == p_outstart){
							std::cout << "Tu passes ici 0" << std::endl;
						}
						*p_out = 130;
						break;
					case 1:
						if (nblbstar == 0){*p_out = 110;
							if (p_out == p_outstart){
								std::cout << "Tu passes ici 1 - 0" << std::endl;
							}
						} // point isole
						if (nblbstar == 1){*p_out = uchar(140);
							if (p_out == p_outstart){
								std::cout << "Tu passes ici ---- p_out = " << int(*p_out) << std::endl;
							}
						} // point terminal 
						break;
					case 2:
						if (p_out == p_outstart){
							std::cout << "Tu passes ici 2" << std::endl;
						}
						*p_out = 120;
						// if (nblbstar == 2){*p_out = 120;} // point de courbe
						// if (nblbstar == 1){*p_out = 170;} // point de surface
						break;
					case 3:
						if (p_out == p_outstart){
							std::cout << "Tu passes ici 3" << std::endl;
						}
						*p_out = 160; // point de bifurcation
						break;
					case 4:
						if (p_out == p_outstart){
							std::cout << "Tu passes ici 4" << std::endl;
						}
						*p_out = 150; // point de croisement
						break;
					default:
						*p_out = 0;
						break;
				}
				//if (p_out == p_outstart){
				//	*p_out = uchar(203);
				//}
			}
		}
	}
	std::vector<ParamT> param; param.clear();
	std::cout << "Valeur du premier point = " << int(imTmp(0,0)) << std::endl;
	DipImage<double> imMean(w,h); imMean.init(0.0);
	param.push_back(1);
	ImMeanFilter(imTmp, imMask, v1,param,imMean);
	dbi.WriteImage(imMean, "imMean.png");
	
	param.clear();
	param.push_back(200.0/3.0);
	param.push_back((Tin)150);
	param.push_back(imTmp);
	ImCompare(imMean, Equal, ssi, param, imTmp);	


	ImCopy(imTmp, imOut);
	dbi.WriteImage(imOut, "imClassInternal.png"); 
	dbi.WriteOut("Leaving ImPointsClassification_WithMask\n");
}

/**
 * ImSkeletonCleaningBranch_WithMask
 * @param imSkeleton
 * @param imMask
 * @param imMask2
 * @param imOut
 */

template<typename Tin, typename Tmask, typename Tout> void ImSkeletonCleaningBranch_WithMask(const DipImage<Tin> & imSkeleton, const DipImage<Tmask> & imMask, const DipImage<Tmask> & imMask2, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImSkeletonCleaningBranch_WithMask\n");
	const int w = imSkeleton.width(), h = imSkeleton.height();
	int i,j,k;

	DipImage<Tin> imSkeletonTmp(w,h); imSkeletonTmp.init(0);
	ImCopy(imSkeleton,imSkeletonTmp);

	DipImage<Tin> imClass(w,h); imClass.init(0);
	ImPointsClassification_WithMask(imSkeleton,imMask,imClass);
	dbi.WriteImage(imClass, "imClass.png");

	DipImage<Tin> imEndPoints(w,h); imEndPoints.init(0);
	ImThreshold(imClass,Tin(139),Tin(141),std::numeric_limits<Tin>::max(), Tin(0),imEndPoints);
	dbi.WriteImage(imEndPoints, "imEndPoints.png");

	int step = 100; // it can be modified with the vessels radius
	int tmpI, tmpJ;
	int sum,sumPatial;
	DipImage<Tin> imClassTmp(w,h); imClassTmp.init(0);

	for(i=0;i<w;i++){
		for(j=0;j<h;j++){
			if(imEndPoints(i,j) !=0 ){
				tmpI = i;
				tmpJ = j;
				for(k=0;k<step;k++){
					imSkeletonTmp(tmpI,tmpJ) = Tin(0);
					imClassTmp.init(0);
					sum = imSkeletonTmp(tmpI-1,tmpJ-1)+imSkeletonTmp(tmpI-1,tmpJ) + imSkeletonTmp(tmpI,tmpJ-1) + imSkeletonTmp(tmpI-1,tmpJ+1) + imSkeletonTmp(tmpI+1,tmpJ-1) + imSkeletonTmp(tmpI+1,tmpJ+1) + imSkeletonTmp(tmpI+1,tmpJ) + imSkeletonTmp(tmpI,tmpJ+1);
					if(sum>255){
						ImPointsClassification_WithMask(imSkeletonTmp,imMask,imClassTmp);
						if(imClassTmp(tmpI-1,tmpJ-1)==Tin(140) || imClassTmp(tmpI-1,tmpJ)==Tin(140) || imClassTmp(tmpI,tmpJ-1)==Tin(140) || imClassTmp(tmpI-1,tmpJ+1)==Tin(140) || imClassTmp(tmpI+1,tmpJ-1)==Tin(140) || imClassTmp(tmpI+1,tmpJ+1)==Tin(140) || imClassTmp(tmpI+1,tmpJ)==Tin(140) || imClassTmp(tmpI,tmpJ+1)==Tin(140)){imSkeletonTmp(tmpI,tmpJ) = Tin(255);}
						//sumPatial =  imSkeletonTmp(tmpI-1,tmpJ-1) + imSkeletonTmp(tmpI-1,tmpJ+1) + imSkeletonTmp(tmpI+1,tmpJ-1) + imSkeletonTmp(tmpI+1,tmpJ+1);
						//if(sumPatial >255){imSkeletonTmp(tmpI,tmpJ) = Tin(255);}
						break;
					}
					else if(sum==0){break;}
					else{
						if(imSkeletonTmp(tmpI-1,tmpJ-1)==Tin(255)){
							tmpI = tmpI-1;
							tmpJ = tmpJ-1;
							//k=k+1;
						}
						else if(imSkeletonTmp(tmpI,tmpJ-1)==Tin(255)){
							tmpI = tmpI;
							tmpJ = tmpJ-1;
						}
						else if(imSkeletonTmp(tmpI+1,tmpJ-1)==Tin(255)){
							tmpI = tmpI+1;
							tmpJ = tmpJ-1;
							//k=k+1;
						}
						else if(imSkeletonTmp(tmpI-1,tmpJ)==Tin(255)){
							tmpI = tmpI-1;
							tmpJ = tmpJ;
						}
						else if(imSkeletonTmp(tmpI+1,tmpJ)==Tin(255)){
							tmpI = tmpI+1;
							tmpJ = tmpJ;
						}
						else if(imSkeletonTmp(tmpI-1,tmpJ+1)==Tin(255)){
							tmpI = tmpI-1;
							tmpJ = tmpJ+1;
							//k=k+1;
						}
						else if(imSkeletonTmp(tmpI,tmpJ+1)==Tin(255)){
							tmpI = tmpI;
							tmpJ = tmpJ+1;
						}
						else if(imSkeletonTmp(tmpI+1,tmpJ+1)==Tin(255)){
							tmpI = tmpI+1;
							tmpJ = tmpJ+1;
							//k=k+1;
						}
					}
				}


			}
		}
	}
ImCopy(imSkeletonTmp,imOut);

	dbi.WriteOut("Leaving ImSkeletonCleaningBranch_WithMask_V2\n");

}



template<typename Tin, typename Tmask> void imageMatting(DipColorImage<Tin> & frame,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting\n");
	const int w = frame.width(), h = frame.height();
	
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	dbi.WriteImage(imInBlue,"imInBlue.png");
	dbi.WriteImage(imInGreen,"imInGreen.png");
	dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	std::vector<ParamT> param; param.clear();
	param.push_back(int(w/50));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");
	
	DipImage<Tmask> imBGTmp(w,h);imBGTmp.init(0);
	DipImage<Tmask> imBGEdge(w,h);imBGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImDilate(imMaskExteral,Fast,param,imBGTmp);
	ImSubImage(imBGTmp,Clip,imMaskExteral,imBGEdge);
	dbi.WriteImage(imBGEdge,"imBGEdge.png");
	
	param.clear();
	param.push_back(int(w/50));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	
	DipImage<Tmask> imFGTmp(w,h);imFGTmp.init(0);
	DipImage<Tmask> imFGEdge(w,h);imFGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImErode(imMaskInteral,Fast,param,imFGTmp);
	ImSubImage(imMaskInteral,Clip,imFGTmp,imFGEdge);
	dbi.WriteImage(imFGEdge,"imFGEdge.png");

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				//~ std::cout << "initial"<<endl;
				//~ std::cout << i << "   " << j <<endl;
				int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				borneInfx = std::max(1,i-w/25);
				borneInfy = std::max(1,j-w/25);
				borneMaxx = std::min(w, i+w/25+1);
				borneMaxy = std::min(h, j+w/25+1);
				std::vector<double> edgeFGStore;edgeFGStore.clear();
				std::vector<double> edgeBGStore;edgeBGStore.clear();
				
				for (int ii=borneInfx; ii<borneMaxx; ii++){
					for (int jj=borneInfy; jj<borneMaxy; jj++){
						if(imFGEdge(ii,jj) != 0){
							edgeFGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));	
						}
						if(imBGEdge(ii,jj) != 0){
							edgeBGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));	
						}
					}
				}
				
				//~ std::vector<double>::iterator it;
				//~ std::cout << "edgeFGStore =" << edgeFGStore.size()<<endl;
				//~ std::cout << "edgeBGStore =" << edgeBGStore.size()<<endl;
				
				//~ cin.get();
				if(edgeFGStore.size() <10){imFGAdd(i,j) = 0;}
				else if(edgeBGStore.size() <10){imFGAdd(i,j) = 255;}
				else{
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());
					//~ edgeFGStore.sort();
					//~ for (it=edgeFGStore.begin(); it!=edgeFGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					//~ edgeBGStore.sort();
					//~ for (it=edgeBGStore.begin(); it!=edgeBGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					while(count !=10){
						//~ it = edgeFGStore.begin();
						double edgeFG = edgeFGStore[0];
						//~ double edgeFG = *it;
						//~ it = edgeBGStore.begin();
						//~ double edgeBG = *it;
						double edgeBG = edgeBGStore[0];
						if(edgeFG <edgeBG){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
							//~ edgeFGStore.pop_front();
						}
						if(edgeBG <edgeFG){
							edgeBGcount  = edgeBGcount +1;
							//~ edgeBGStore.pop_front();
							edgeBGStore.erase(edgeBGStore.begin());
						}
						
						count = count +1;
					}
					if(edgeFGcount > edgeBGcount)
						imFGAdd(i,j) = 255;
				}					
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");
	//~ for (int i= 0; i<w; i++){
		//~ for (int j=0; j<h; j++){
			//~ if(imMaskUnknown(i,j) !=0){
				//~ int bgdistance = 0;
				//~ int fgdistance = 0;
				//~ int bgCount = 0;
				//~ int bgEffectCount = 0;
				//~ int fgCount = 0;
				//~ int fgEffectCount = 0;
				//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				//~ borneInfx = std::max(0,i-w/30);
				//~ borneInfy = std::max(0,j-w/30);
				//~ borneMaxx = std::min(w, i+w/30+1);
				//~ borneMaxy = std::min(h, j+w/30+1);
				//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
					//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
						//~ if(imBGEdge(ii,jj) !=0){
							//~ bgCount = bgCount +1;
							//~ bgdistance = bgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ bgEffectCount = bgEffectCount +1;
							
						//~ }
						//~ if(imFGEdge(ii,jj) !=0){
							//~ fgdistance = fgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ fgCount = fgCount + 1;
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ fgEffectCount = fgEffectCount +1;
							
						//~ }
					//~ }
				//~ }
				//~ if(double(fgdistance)/fgCount < double(bgdistance)/bgCount){
					//~ imFG(i,j) = 255;
				//~ }
				//~ else{
					//~ imBG(i,j) = 255;
				//~ }
			//~ }
		//~ }
	//~ }
	//~ 
	//~ dbi.WriteImage(imFG,"imFG_fii.png");
	
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ DipImage<Tmask> imProcessingPoint(w,h);imProcessingPoint.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp(w,h);imMaskInteralTemp.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp2(w,h);imMaskInteralTemp2.init(0);
	//~ std::vector<ParamT> param; param.clear();
	//~ param.push_back(int(1));
	//~ ImDilate(imMaskInteral,Fast,param,imMaskInteralTemp2);
	//~ ImCopy(imMaskInteral,imMaskInteralTemp);
	//~ ImSubImage(imMaskInteralTemp2,imMaskTotal,Clip,imMaskInteralTemp,imProcessingPoint);
	//~ dbi.WriteImage(imProcessingPoint,"imProcessingPoint.png");
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){			
			//~ if(imProcessingPoint(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	
	
	
	//~ DipImage<Tin> imMaskUnknownInvert(w,h); imMaskUnknownInvert.init(0);
	//~ ImInvert(imMaskUnknown,imMaskUnknownInvert);
	//~ dbi.WriteImage(imMaskUnknownInvert, "imMaskUnknownInvert.png");
	//~ 
	//~ DipImage<Tin> imDistanceVN(w,h); imDistanceVN.init(0);
	//~ ImMorphoDistance(imMaskUnknownInvert,imDistanceVN);
	//~ dbi.WriteImage(imDistanceVN, "imDistanceVN.png");

	//~ DipImage<Tin> imSkeleton(w,h); imSkeleton.init(0);
	//~ ImThinning_WithMask(imDistanceVN, imMaskUnknown, imSkeleton);
	//~ dbi.WriteImage(imSkeleton, "imSkeleton.png");
	//~ 
	//~ DipImage<Tin> imSkeletonCleanSecond(w,h); imSkeletonCleanSecond.init(0);
	//~ ImSkeletonCleaningBranch_WithMask(imSkeleton,imMaskUnknown,imMaskUnknown,imSkeletonCleanSecond);
	//~ dbi.WriteImage(imSkeletonCleanSecond, "imSkeletonCleanSecond.png");
	
	//~ 
	//~ DipImage<uchar> imEdge(w,h);
	//~ MaskedgeDetection(frame,imMaskUnknown, imEdge);
	//~ 
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){
			//~ //if(imSkeleton(i,j) == 255){
			//~ if(imEdge(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	//~ 
	//~ 
	//~ DipImage<int> imPatialBG(w,h); imPatialBG.init(0);
	//~ int count;
	//~ for(int i = 0; i < pointProcessing.size();i++){
		//~ count = count + 1;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/30);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/30);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/30+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/30+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ 
			//~ }
		//~ }
		
		//~ //kmeans algorithm
		//~ Mat src(frame, Rect(borneInfx,borneInfy,borneMaxx - borneInfx,borneMaxy - borneInfy ));
		//~ Mat samples(src.rows * src.cols, 3, CV_32F);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
				//~ for( int z = 0; z < 3; z++)
					//~ samples.at<float>(y + x*src.rows, z) = src.at<Vec3b>(y,x)[z];
//~ 
//~ 
		//~ int clusterCount = 2;
		//~ Mat labels;
		//~ int attempts = 5;
		//~ Mat centers;
		//~ kmeans(samples, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers );
//~ 
//~ 
		//~ DipColorImage<uchar> new_image(src.cols,src.rows);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
			//~ { 
			//~ int cluster_idx = labels.at<int>(y + x*src.rows,0);
			//~ new_image(x,y)[0] = centers.at<float>(cluster_idx, 0);
			//~ new_image(x,y)[1] = centers.at<float>(cluster_idx, 1);
			//~ new_image(x,y)[2] = centers.at<float>(cluster_idx, 2);
		//~ }
		
		//~ dbi.WriteImage(new_image, boost::lexical_cast<string>(count) + "new_image.png");
		

		
		
		
		//~ std::vector<std::vector<int> > data;
		//~ std::vector<std::vector<int> > results;
		//~ std::vector<int> colorInfo;colorInfo.clear();
		//~ std::vector<int> sizeG;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/35);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/35);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/35+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/35+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ colorInfo.clear();
				//~ colorInfo.push_back(imInGreen(ii,jj));
				//~ data.push_back(colorInfo);
			//~ }
		//~ }
		//~ results = KMeans(data,2,sizeG);
		//~ dbi.WriteInfo("results 1 %d",results[0][0]);
		//~ dbi.WriteInfo("results 2 %d",results[1][0]);
		//~ bool isPointBG;
		//~ if(borneInfx < w/2)isPointBG = false;
		//~ else isPointBG = true;
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){				
				//~ if(abs(imInGreen(ii,jj)-results[0][0]) > abs(imInGreen(ii,jj)-results[0][0])
			//~ }
		//~ }
		
		
	//~ }
	
	
	dbi.WriteOut("Leaving imagematting\n");
}





template<typename Tin>
std::vector<std::vector<Tin> > KMeans(const std::vector<std::vector<Tin> > & data, const int & K, std::vector<int> & SizeGroups)
{
	dbi.WriteEnter("Entering Kmeans ... \n");
	uint N = data.size();
	uint D = data[0].size();
	
	std::vector<Tin> MeanInit(D,0);
	std::vector<std::vector<Tin> > Means(K,MeanInit);
	int Classe(0);
	
	// Initialisation aléatoire des centroides
	for (int k=0; k<K; k++)
		Means[k] = data[std::rand()%N];
	
	int it = 0, I = 100;
	bool Changes = true;
	
	// Tant qu'il y a des changements et que le nombre max d'iterations n'est pas atteint
	while ((Changes==true) && (it<I))
	{
		it++;
		std::vector<std::vector<Tin> > NewMeans(K,MeanInit);
		
		std::vector<int> inc(K,0);
		
		// Pour chaque point
		for (uint n=0; n<N; n++)
		{
			// Trouver le centroid le plus proche
			int DistMin = 10000000;
			for (int k = 0; k<K; k++)
			{
				double Dist = 0.0;
				for (uint d=0; d<D; d++)
					Dist+=(data[n][d]-Means[k][d])*(data[n][d]-Means[k][d]);
				Dist = sqrt(Dist);
				if (Dist<DistMin)
				{	
					DistMin = Dist;
					Classe = k;
				}
			}
			
			// Attribuer le point à la classe la plus proche
			for (uint d=0; d<D; d++)
				NewMeans[Classe][d]+=data[n][d];
			inc[Classe]++;
		}
		
		// Pour chaque groupe
		for (int k=0; k<K; k++)
		{
			if (inc[k]==0)
				inc[k]=1;
			// Mettre à jour la moyenne
			for (uint d=0; d<D; d++)
				NewMeans[k][d]=NewMeans[k][d]/inc[k];
		}
		
		// Evaluer les changements	
		int nChanges = K;
		for (int k=0; k<K; k++)
		{
			double Diff = 0.0;
			for (uint d=0; d<D; d++)
				Diff+=(NewMeans[k][d]-Means[k][d])*(NewMeans[k][d]-Means[k][d]);
			Diff = sqrt(Diff);
			if (Diff < 1)
				nChanges--;
		}
		if (nChanges == 0)
			Changes = false;
			
		Means = NewMeans;	
		SizeGroups = inc;
	
	}	 
	
	dbi.WriteOut("Leaving Kmeans. \n");
	
	return(Means);
}


template<typename Tin, typename Tmask,typename Tedge> void MaskedgeDetection(DipColorImage<Tin> & frame,DipImage<Tmask> & imMask, DipImage<Tedge> & imEdge){
	int w  = frame.width();
	int h  = frame.height();
	DipImage<uchar> ImGrayscale(w,h);ImGrayscale.init(0);
	DipImage<uchar> ImGrayscaleMask(w,h);ImGrayscaleMask.init(0);
	cvtColor(frame,ImGrayscale,CV_BGR2GRAY);
	blur(ImGrayscale, ImGrayscale, Size(3,3) );
	DipImage<uchar> edges(w,h);edges.init(0);
	DipImage<uchar> _img(w,h);_img.init(0);
	ImCopy(ImGrayscale,imMask,ImGrayscaleMask);
	
	double otsu_thresh_val = cv::threshold(ImGrayscaleMask, _img, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
	double high_thresh_val  = otsu_thresh_val,lower_thresh_val = 0.5 * otsu_thresh_val;
	Canny(ImGrayscaleMask, imEdge,lower_thresh_val,high_thresh_val);
	dbi.WriteImage(imEdge, "imEdge.png");
}

template<typename Tin, typename Tmask> void maskDetection(DipColorImage<Tin> & frame,std::string & path_model, DipImage<Tmask> & imMaskOut){
	const int w  = frame.width(); const int h  = frame.height();
	std::vector<ParamT> param; param.clear();
	DipImage<uchar> mask(w,h);mask.init(255);
	DipColorImage<uchar> frameLight(w,h);frameLight.init(0);
	DipColorImage<uchar> frameLighttmp(w,h);frameLighttmp.init(0);
	DipColorImage<uchar> frameEnhanced(w,h);frameEnhanced.init(0);
	ImEnhancementLighting(frame, mask, int(w/15),frameLighttmp, frameEnhanced);
	ImLightingCorrection(frame,mask,int(w/15),frameLight);
	dbi.WriteImage(frameLight,"frameLight.png");
	dbi.WriteImage(frameEnhanced,"frameEnhanced.png");

	//~ //YCrcb image
	//~ DipColorImage<uchar> frameYCrCb(w,h);frameYCrCb.init(0);
	//~ cv::cvtColor(frame, frameYCrCb, COLOR_BGR2YCrCb );
	//~ dbi.WriteImage(frameYCrCb,"frameYCrCb.png");
	//~ DipImage<uchar> imY(w,h);imY.init(0);
	//~ DipImage<uchar> imCr(w,h);imCr.init(0);
	//~ DipImage<uchar> imCb(w,h);imCb.init(0);
	//~ ImColorBandSeparation(frameYCrCb,imY,imCr,imCb);
	//~ 
	//~ //Yuv image
	//~ DipColorImage<uchar> frameYUV(w,h);frameYUV.init(0);
	//~ cv::cvtColor(frame, frameYUV, COLOR_BGR2YUV );
	//~ dbi.WriteImage(frameYUV,"frameYUV.png");
	//~ DipImage<uchar> imYUVY(w,h);imYUVY.init(0);
	//~ DipImage<uchar> imU(w,h);imU.init(0);
	//~ DipImage<uchar> imV(w,h);imV.init(0);
	//~ ImColorBandSeparation(frameYUV,imYUVY,imU,imV);
	//~ 
	//~ 
	//~ //face skin classification
	//~ DipImage<uchar> imSkin(w,h);imSkin.init(0);
	//~ for(int i = 0;i < w; i++){
		//~ for( int j = 0;j < h; j++){
			//~ if(imU(i,j) < 130 && imU(i,j) > 80){
				//~ if(imV(i,j) <200 && imV(i,j)>136){
					//~ if(imInRed(i,j) > 80 && imInGreen(i,j) > 30 && imInBlue(i,j) > 15){
						//~ if(abs(imInRed(i,j) - imInGreen(i,j)) > 15){
							//~ imSkin(i,j) = 255;
						//~ }
					//~ }
				//~ }
			//~ }
		//~ }
	//~ }
	//~ dbi.WriteImage(imSkin,"imSkin.png");
			
	DipImage<uchar> imFaceMask(w,h);imFaceMask.init(0);
	cv::Rect roi_head;
	roi_head = detectFaceAndBody(frameLight,imFaceMask,path_model);
	//@todo add detection eye to validate the mask
			
	dbi.WriteImage(imFaceMask,"imFaceMask.png");
	
			
	DipImage<uchar> imFaceMaskDilate(w,h); imFaceMaskDilate.init(0);
	param.clear();
	param.push_back(int(w/40));
	ImDilate(imFaceMask,Fast,param,imFaceMaskDilate);
	dbi.WriteImage(imFaceMaskDilate,"imFaceMaskDilate.png");
	
	
	DipImage<uchar> imGrabMask(w,h);imGrabMask.init(0);
	DipImage<uchar> imFGSure(w,h);imFGSure.init(0);
	DipImage<uchar> imFGSureHalf(w,h);imFGSureHalf.init(0);
	for(int i = 0; i< w;++i){
		for(int j = 0; j< h; ++j){
			if(imFaceMaskDilate(i,j) == 0) {imGrabMask(i,j) = cv::GC_BGD;}
			else{
					if( i> roi_head.x + 1*roi_head.width/4 && i < roi_head.x + 3* roi_head.width/4 && j > roi_head.y + roi_head.height/4 && j < roi_head.y + 3*roi_head.height/4 ){								
						imGrabMask(i,j) = cv::GC_FGD;
						imFGSure(i,j) = 255;
						imFGSureHalf(i,j) = 255;
					}
					else if(i> roi_head.x + 3*roi_head.width/8 && i < roi_head.x + 5* roi_head.width/8 && j >= roi_head.y + 3*roi_head.height/4 && j <= roi_head.y + 5*roi_head.height/4){ 
						imGrabMask(i,j) = cv::GC_FGD;
						imFGSure(i,j) = 255;
						imFGSureHalf(i,j) = 255;
					}
					else if(i> std::max(roi_head.x - roi_head.width/4,0) && i < roi_head.x + 5*roi_head.width/4 && j > roi_head.y + 5*roi_head.height/4){ 
						imGrabMask(i,j) = cv::GC_FGD;
						imFGSure(i,j) = 255;
					}
					
					else imGrabMask(i,j) = cv::GC_PR_FGD;
			}
		}
	}
	int size = roi_head.width /4; 
	dbi.WriteImage(imGrabMask,"imGrabMask.png");
	dbi.WriteImage(imFGSure,"imFGSure.png");
	cv::Mat bgModel,fgModel;
	DipImage<uchar> imGrabMaskTmp(w,h);
	cv::Rect rectangle(1,1,h-1,w-1);
	ImCopy(imGrabMask,imGrabMaskTmp);
	
	//~ // RGB TO GRAYSCALE IMAGE
	//~ DipImage<uchar> frameEnhanced_gray(w,h);
	//~ cv::cvtColor(frameLight, frameEnhanced_gray, COLOR_BGR2GRAY);
	//~ // HISTGRAM EQUALIZATION
	//~ equalizeHist(frameEnhanced_gray, frameEnhanced_gray);
	//~ dbi.WriteImage(frameEnhanced_gray, "frameEnhanced_gray.png");
	//~ 
	//~ DipImage<uchar> imGreen(w,h);
	//~ DipImage<uchar> imRed(w,h);
	//~ DipImage<uchar> imBlue(w,h);
	//~ //ImColorBandSeparation(frameEnhanced,imBlue,imGreen,imRed);
	//~ ImColorBandSeparation(frame,imBlue,imGreen,imRed);
	//~ //equalizeHist(imGreen, imGreen);
	//~ //equalizeHist(imRed, imRed);
	//~ DipColorImage<uchar> frameEnhanced_new(w,h);frameEnhanced_new.init(0);
	//~ ImColorBandMerge(imGreen,frameEnhanced_gray,imRed,frameEnhanced_new);
	//~ dbi.WriteImage(frameEnhanced_new, "frameEnhanced_new.png");
	grabCut(frameEnhanced, imGrabMaskTmp,rectangle,bgModel,fgModel,10, cv::GC_INIT_WITH_MASK); 
	
	DipImage<uchar> imFG(w,h);imFG.init(0);
	for(int i = 0; i< w;++i){
		for(int j = 0; j< h; ++j){
			if(imGrabMaskTmp(i,j) == GC_FGD || imGrabMaskTmp(i,j) == 3 ) {imFG(i,j) = 255;}
	
		}
	}
	dbi.WriteImage(imFG, "imFG.png");
	
	
	//extract the largest one
	DipImage<double> imFGD(w,h); imFGD.init(0);
	DipImage<uchar> imFGLarget(w,h); imFGLarget.init(0);
	DipImage<int> imFGLabel(w,h); imFGLabel.init(0);
	ImLabelFlatZones_WithCriterion(imFG,imFGD,Area,imFGLabel);
	std::pair<double,double> minmax = MeasMinMax(imFGLabel,mask);
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
	
	
	//grab cut again
	DipImage<uchar> imBGAgain(w,h); imBGAgain.init(0);
	param.clear();
	param.push_back(int(w/100));
	ImDilate(imBGFinal,Fast,param,imBGAgain);
	dbi.WriteImage(imBGAgain,"imBGAgain.png");
	DipImage<uchar> imFGAgain(w,h); imFGAgain.init(0);
	DipImage<uchar> imFGAgain2(w,h); imFGAgain2.init(0);
	param.clear();
	param.push_back(int(w/15));
	ImErode(imBGFinal,Fast,param,imFGAgain2);
	ImAddImage(imFGAgain2,Clip,imFGSureHalf,imFGAgain);
	dbi.WriteImage(imFGAgain,"imFGAgain.png");
	DipImage<uchar> imGrabMaskAgain(w,h);imGrabMaskAgain.init(0);
	for(int i = 0; i< w;++i){
		for(int j = 0; j< h; ++j){
			if(imBGAgain(i,j) == 0) {imGrabMaskAgain(i,j) = cv::GC_BGD;}
			else if(imFGAgain(i,j) == 255){imGrabMaskAgain(i,j) = cv::GC_FGD;}
			else imGrabMaskAgain(i,j) = cv::GC_PR_FGD;
		}
	}
	
	
	grabCut(frame, imGrabMaskAgain,rectangle,bgModel,fgModel,10, cv::GC_INIT_WITH_MASK); 
	
	DipImage<uchar> imFGAgainTure(w,h);imFGAgainTure.init(0);
	for(int i = 0; i< w;++i){
		for(int j = 0; j< h; ++j){
			if(imGrabMaskAgain(i,j) == GC_FGD || imGrabMaskAgain(i,j) == 3 ) {imFGAgainTure(i,j) = 255;}
	
		}
	}
	dbi.WriteImage(imFGAgainTure,"imFGAgainTure.png");
	param.clear();
	param.push_back(int(w/25));
	DipImage<uchar> imMaskClose(w,h);imMaskClose.init(0);
	DipImage<uchar> imMaskOpen(w,h);imMaskOpen.init(0);
	ImClose(imFGAgainTure,Fast,param,imMaskClose);
	ImOpen(imMaskClose,Fast,param,imMaskOpen);	
	dbi.WriteImage(imMaskOpen,"imMaskOpen.png");
	//~ ImCopy(imFGAgainTure,imMaskOut);
	ImCopy(imMaskOpen,imMaskOut);
}



template<typename Tin, typename Tmask> void imageMatting2(DipColorImage<Tin> & frame,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting2\n");
	const int w = frame.width(), h = frame.height();
	
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	dbi.WriteImage(imInBlue,"imInBlue.png");
	dbi.WriteImage(imInGreen,"imInGreen.png");
	dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	std::vector<ParamT> param; param.clear();
	param.push_back(int(w/50));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");
	
	DipImage<Tmask> imBGTmp(w,h);imBGTmp.init(0);
	DipImage<Tmask> imBGEdge(w,h);imBGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImDilate(imMaskExteral,Fast,param,imBGTmp);
	ImSubImage(imBGTmp,Clip,imMaskExteral,imBGEdge);
	dbi.WriteImage(imBGEdge,"imBGEdge.png");
	
	param.clear();
	param.push_back(int(w/50));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	
	DipImage<Tmask> imFGTmp(w,h);imFGTmp.init(0);
	DipImage<Tmask> imFGEdge(w,h);imFGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImErode(imMaskInteral,Fast,param,imFGTmp);
	ImSubImage(imMaskInteral,Clip,imFGTmp,imFGEdge);
	dbi.WriteImage(imFGEdge,"imFGEdge.png");

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	//~ std::vector<std::pair<int,int> > edgeFGPosition;edgeFGPosition.clear();
	//~ std::vector<std::pair<int,int> > edgeBGPosition;edgeBGPosition.clear();
	//~ 
	//~ for (int i= 0; i<w; i++){
		//~ for (int j=0; j<h; j++){
			//~ if(imFGEdge(i,j) != 0){
				//~ edgeFGPosition.push_back(std::make_pair(i,j));
			//~ }
			//~ if(imBGEdge(i,j) != 0){
				//~ edgeBGPosition.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				std::vector<double> edgeFGStore;edgeFGStore.clear();
				std::vector<double> edgeBGStore;edgeBGStore.clear();
				DipImage<Tmask> imPoint(w,h);imPoint.init(0);
				DipImage<Tmask> imPointCandidate(w,h);imPointCandidate.init(0);
				DipImage<Tmask> imPoint_2(w,h);imPoint_2.init(0);
				std::cout << i << "   " << j <<endl;
				imPoint(i,j) = 255;
				param.clear();
				param.push_back(int(20));
				while((edgeFGStore.size() < 10) || (edgeBGStore.size() < 10)){
					imPointCandidate.init(0);
					ImDilate(imPoint,Fast,param,imPoint_2);
					ImSubImage(imPoint_2,Clip,imPoint,imPointCandidate);
					for (int ii=0; ii<w; ii++){
						for (int jj=0; jj<h; jj++){
							if(imFGEdge(ii,jj) != 0 && imPointCandidate(ii,jj) != 0){
								edgeFGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));
							}
							if(imBGEdge(ii,jj) != 0 && imPointCandidate(ii,jj) != 0){
								edgeBGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));
							}
						}
					}
				std::cout << edgeFGStore.size() << "   " << edgeBGStore.size() <<endl;
				ImCopy(imPoint_2,imPoint);
				imPoint_2.init(0);
									
				}			
				

				std::sort(edgeFGStore.begin(),edgeFGStore.end());
				std::sort(edgeBGStore.begin(),edgeBGStore.end());

				int count = 0;
				int edgeFGcount = 0;
				int edgeBGcount = 0;
				while(count !=10){
					double edgeFG = edgeFGStore[0];
	
					double edgeBG = edgeBGStore[0];
					if(edgeFG <edgeBG){
						edgeFGcount  = edgeFGcount +1;
						edgeFGStore.erase(edgeFGStore.begin());
					}
					if(edgeBG <edgeFG){
						edgeBGcount  = edgeBGcount +1;	
						edgeBGStore.erase(edgeBGStore.begin());
					}
					
					count = count +1;
				}
				if(edgeFGcount >= edgeBGcount)
					imFGAdd(i,j) = 255;					
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");
	
	
	
	dbi.WriteOut("Leaving imagematting2\n");
}


template<typename Tin, typename Tmask> void imageMatting_initial(DipColorImage<Tin> & frame,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting_initial\n");
	const int w = frame.width(), h = frame.height();
	
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	dbi.WriteImage(imInBlue,"imInBlue.png");
	dbi.WriteImage(imInGreen,"imInGreen.png");
	dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	std::vector<ParamT> param; param.clear();
	param.push_back(int(w/50));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");
	
	DipImage<Tmask> imBGTmp(w,h);imBGTmp.init(0);
	DipImage<Tmask> imBGEdge(w,h);imBGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImDilate(imMaskExteral,Fast,param,imBGTmp);
	ImSubImage(imBGTmp,Clip,imMaskExteral,imBGEdge);
	dbi.WriteImage(imBGEdge,"imBGEdge.png");
	
	param.clear();
	param.push_back(int(w/70));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	
	DipImage<Tmask> imFGTmp(w,h);imFGTmp.init(0);
	DipImage<Tmask> imFGEdge(w,h);imFGEdge.init(0);
	param.clear();
	param.push_back(int(1));
	ImErode(imMaskInteral,Fast,param,imFGTmp);
	ImSubImage(imMaskInteral,Clip,imFGTmp,imFGEdge);
	dbi.WriteImage(imFGEdge,"imFGEdge.png");

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				//~ std::cout << "initial"<<endl;
				//~ std::cout << i << "   " << j <<endl;
				int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				borneInfx = std::max(1,i-w/25);
				borneInfy = std::max(1,j-w/25);
				borneMaxx = std::min(w, i+w/25+1);
				borneMaxy = std::min(h, j+w/25+1);
				std::vector<double> edgeFGStore;edgeFGStore.clear();
				std::vector<double> edgeBGStore;edgeBGStore.clear();
				
				for (int ii=borneInfx; ii<borneMaxx; ii++){
					for (int jj=borneInfy; jj<borneMaxy; jj++){
						if(imFGEdge(ii,jj) != 0){
							edgeFGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));	
						}
						if(imBGEdge(ii,jj) != 0){
							edgeBGStore.push_back(double(abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j)) + abs(imInBlue(ii,jj)- imInBlue(i,j)) * abs(imInBlue(ii,jj)- imInBlue(i,j)) + abs(imInRed(ii,jj)- imInRed(i,j)) *abs(imInRed(ii,jj)- imInRed(i,j))));	
						}
					}
				}
				
				//~ std::vector<double>::iterator it;
				//~ std::cout << "edgeFGStore =" << edgeFGStore.size()<<endl;
				//~ std::cout << "edgeBGStore =" << edgeBGStore.size()<<endl;
				
				//~ cin.get();
				if(edgeFGStore.size() <10){imFGAdd(i,j) = 0;}
				else if(edgeBGStore.size() <10){imFGAdd(i,j) = 255;}
				else{
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());
					//~ edgeFGStore.sort();
					//~ for (it=edgeFGStore.begin(); it!=edgeFGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					//~ edgeBGStore.sort();
					//~ for (it=edgeBGStore.begin(); it!=edgeBGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					while(count !=10){
						//~ it = edgeFGStore.begin();
						double edgeFG = edgeFGStore[0];
						//~ double edgeFG = *it;
						//~ it = edgeBGStore.begin();
						//~ double edgeBG = *it;
						double edgeBG = edgeBGStore[0];
						if(edgeFG <edgeBG){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
							//~ edgeFGStore.pop_front();
						}
						if(edgeBG <edgeFG){
							edgeBGcount  = edgeBGcount +1;
							//~ edgeBGStore.pop_front();
							edgeBGStore.erase(edgeBGStore.begin());
						}
						
						count = count +1;
					}
					if(edgeFGcount > edgeBGcount)
						imFGAdd(i,j) = 255;
				}					
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");	
	
	dbi.WriteOut("Leaving imageMatting_initial\n");
}



template<typename Tin, typename Tmask> void imageMatting_2(DipColorImage<Tin> & frame,DipColorImage<Tin> & frame_previous,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting_2\n");
	const int w = frame.width(), h = frame.height();
	
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	//~ dbi.WriteImage(imInBlue,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen,"imInGreen.png");
	//~ dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<uchar> imInGreen_p(w,h);
	DipImage<uchar> imInRed_p(w,h);
	DipImage<uchar> imInBlue_p(w,h);
	ImColorBandSeparation(frame_previous,imInBlue_p,imInGreen_p,imInRed_p);
	//~ dbi.WriteImage(imInBlue_p,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen_p,"imInGreen.png");
	//~ dbi.WriteImage(imInRed_p,"imInRed.png");
	
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	std::vector<ParamT> param; param.clear();
	param.push_back(int(w/100));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");
	
	DipImage<Tmask> imBGTmp(w,h);imBGTmp.init(0);
	DipImage<Tmask> imBGEdge(w,h);imBGEdge.init(0);
	param.clear();
	param.push_back(int(2));
	ImDilate(imMaskExteral,Fast,param,imBGTmp);
	ImSubImage(imBGTmp,Clip,imMaskExteral,imBGEdge);
	dbi.WriteImage(imBGEdge,"imBGEdge.png");
	
	param.clear();
	param.push_back(int(w/100));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	
	DipImage<Tmask> imFGTmp(w,h);imFGTmp.init(0);
	DipImage<Tmask> imFGEdge(w,h);imFGEdge.init(0);
	param.clear();
	param.push_back(int(2));
	ImErode(imMaskInteral,Fast,param,imFGTmp);
	ImSubImage(imMaskInteral,Clip,imFGTmp,imFGEdge);
	dbi.WriteImage(imFGEdge,"imFGEdge.png");

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				//~ std::cout << "initial"<<endl;
				//~ std::cout << i << "   " << j <<endl;
				int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				borneInfx = std::max(1,i-w/40);
				borneInfy = std::max(1,j-w/40);
				borneMaxx = std::min(w, i+w/40+1);
				borneMaxy = std::min(h, j+w/40+1);
				std::vector<int> edgeFGStore;edgeFGStore.clear();
				std::vector<int> edgeBGStore;edgeBGStore.clear();
				
				for (int ii=borneInfx; ii<borneMaxx; ii++){
					for (int jj=borneInfy; jj<borneMaxy; jj++){
						if(imFGEdge(ii,jj) != 0){
							edgeFGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));	
							//~ edgeFGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) ));	
						}
						if(imBGEdge(ii,jj) != 0){
							edgeBGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));
							//~ edgeBGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) ));
						}
					}
				}
				
				//~ std::vector<double>::iterator it;
				//~ std::cout << "edgeFGStore =" << edgeFGStore.size()<<endl;
				//~ std::cout << "edgeBGStore =" << edgeBGStore.size()<<endl;
				
				//~ cin.get();
				if(edgeFGStore.size() <20){imFGAdd(i,j) = 0;}
				else if(edgeBGStore.size() <20){imFGAdd(i,j) = 255;}
				else{
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());
					//~ edgeFGStore.sort();
					//~ for (it=edgeFGStore.begin(); it!=edgeFGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					//~ edgeBGStore.sort();
					//~ for (it=edgeBGStore.begin(); it!=edgeBGStore.end(); ++it){
						//~ std::cout << ' ' << *it;
					//~ }
					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					while(count !=20){
						//~ it = edgeFGStore.begin();
						double edgeFG = edgeFGStore[0];
						//~ double edgeFG = *it;
						//~ it = edgeBGStore.begin();
						//~ double edgeBG = *it;
						double edgeBG = edgeBGStore[0];
						if(edgeFG <=edgeBG){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
							//~ edgeFGStore.pop_front();
						}
						if(edgeBG <edgeFG){
							edgeBGcount  = edgeBGcount +1;
							//~ edgeBGStore.pop_front();
							edgeBGStore.erase(edgeBGStore.begin());
						}
						
						count = count +1;
					}
					if(edgeFGcount >= edgeBGcount)
						imFGAdd(i,j) = 255;
				}					
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");
	//~ for (int i= 0; i<w; i++){
		//~ for (int j=0; j<h; j++){
			//~ if(imMaskUnknown(i,j) !=0){
				//~ int bgdistance = 0;
				//~ int fgdistance = 0;
				//~ int bgCount = 0;
				//~ int bgEffectCount = 0;
				//~ int fgCount = 0;
				//~ int fgEffectCount = 0;
				//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				//~ borneInfx = std::max(0,i-w/30);
				//~ borneInfy = std::max(0,j-w/30);
				//~ borneMaxx = std::min(w, i+w/30+1);
				//~ borneMaxy = std::min(h, j+w/30+1);
				//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
					//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
						//~ if(imBGEdge(ii,jj) !=0){
							//~ bgCount = bgCount +1;
							//~ bgdistance = bgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ bgEffectCount = bgEffectCount +1;
							
						//~ }
						//~ if(imFGEdge(ii,jj) !=0){
							//~ fgdistance = fgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ fgCount = fgCount + 1;
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ fgEffectCount = fgEffectCount +1;
							
						//~ }
					//~ }
				//~ }
				//~ if(double(fgdistance)/fgCount < double(bgdistance)/bgCount){
					//~ imFG(i,j) = 255;
				//~ }
				//~ else{
					//~ imBG(i,j) = 255;
				//~ }
			//~ }
		//~ }
	//~ }
	//~ 
	//~ dbi.WriteImage(imFG,"imFG_fii.png");
	
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ DipImage<Tmask> imProcessingPoint(w,h);imProcessingPoint.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp(w,h);imMaskInteralTemp.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp2(w,h);imMaskInteralTemp2.init(0);
	//~ std::vector<ParamT> param; param.clear();
	//~ param.push_back(int(1));
	//~ ImDilate(imMaskInteral,Fast,param,imMaskInteralTemp2);
	//~ ImCopy(imMaskInteral,imMaskInteralTemp);
	//~ ImSubImage(imMaskInteralTemp2,imMaskTotal,Clip,imMaskInteralTemp,imProcessingPoint);
	//~ dbi.WriteImage(imProcessingPoint,"imProcessingPoint.png");
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){			
			//~ if(imProcessingPoint(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	
	
	
	//~ DipImage<Tin> imMaskUnknownInvert(w,h); imMaskUnknownInvert.init(0);
	//~ ImInvert(imMaskUnknown,imMaskUnknownInvert);
	//~ dbi.WriteImage(imMaskUnknownInvert, "imMaskUnknownInvert.png");
	//~ 
	//~ DipImage<Tin> imDistanceVN(w,h); imDistanceVN.init(0);
	//~ ImMorphoDistance(imMaskUnknownInvert,imDistanceVN);
	//~ dbi.WriteImage(imDistanceVN, "imDistanceVN.png");

	//~ DipImage<Tin> imSkeleton(w,h); imSkeleton.init(0);
	//~ ImThinning_WithMask(imDistanceVN, imMaskUnknown, imSkeleton);
	//~ dbi.WriteImage(imSkeleton, "imSkeleton.png");
	//~ 
	//~ DipImage<Tin> imSkeletonCleanSecond(w,h); imSkeletonCleanSecond.init(0);
	//~ ImSkeletonCleaningBranch_WithMask(imSkeleton,imMaskUnknown,imMaskUnknown,imSkeletonCleanSecond);
	//~ dbi.WriteImage(imSkeletonCleanSecond, "imSkeletonCleanSecond.png");
	
	//~ 
	//~ DipImage<uchar> imEdge(w,h);
	//~ MaskedgeDetection(frame,imMaskUnknown, imEdge);
	//~ 
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){
			//~ //if(imSkeleton(i,j) == 255){
			//~ if(imEdge(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	//~ 
	//~ 
	//~ DipImage<int> imPatialBG(w,h); imPatialBG.init(0);
	//~ int count;
	//~ for(int i = 0; i < pointProcessing.size();i++){
		//~ count = count + 1;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/30);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/30);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/30+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/30+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ 
			//~ }
		//~ }
		
		//~ //kmeans algorithm
		//~ Mat src(frame, Rect(borneInfx,borneInfy,borneMaxx - borneInfx,borneMaxy - borneInfy ));
		//~ Mat samples(src.rows * src.cols, 3, CV_32F);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
				//~ for( int z = 0; z < 3; z++)
					//~ samples.at<float>(y + x*src.rows, z) = src.at<Vec3b>(y,x)[z];
//~ 
//~ 
		//~ int clusterCount = 2;
		//~ Mat labels;
		//~ int attempts = 5;
		//~ Mat centers;
		//~ kmeans(samples, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers );
//~ 
//~ 
		//~ DipColorImage<uchar> new_image(src.cols,src.rows);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
			//~ { 
			//~ int cluster_idx = labels.at<int>(y + x*src.rows,0);
			//~ new_image(x,y)[0] = centers.at<float>(cluster_idx, 0);
			//~ new_image(x,y)[1] = centers.at<float>(cluster_idx, 1);
			//~ new_image(x,y)[2] = centers.at<float>(cluster_idx, 2);
		//~ }
		
		//~ dbi.WriteImage(new_image, boost::lexical_cast<string>(count) + "new_image.png");
		

		
		
		
		//~ std::vector<std::vector<int> > data;
		//~ std::vector<std::vector<int> > results;
		//~ std::vector<int> colorInfo;colorInfo.clear();
		//~ std::vector<int> sizeG;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/35);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/35);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/35+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/35+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ colorInfo.clear();
				//~ colorInfo.push_back(imInGreen(ii,jj));
				//~ data.push_back(colorInfo);
			//~ }
		//~ }
		//~ results = KMeans(data,2,sizeG);
		//~ dbi.WriteInfo("results 1 %d",results[0][0]);
		//~ dbi.WriteInfo("results 2 %d",results[1][0]);
		//~ bool isPointBG;
		//~ if(borneInfx < w/2)isPointBG = false;
		//~ else isPointBG = true;
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){				
				//~ if(abs(imInGreen(ii,jj)-results[0][0]) > abs(imInGreen(ii,jj)-results[0][0])
			//~ }
		//~ }
		
		
	//~ }
	
	
	dbi.WriteOut("Leaving imagematting\n");
}




template<typename Tin, typename Tmask> void imageMatting_3(DipColorImage<Tin> & frame,DipColorImage<Tin> & frame_previous,DipColorImage<Tin> & frame_0,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskInitial,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting_3\n");
	const int w = frame.width(), h = frame.height();
	std::vector<ParamT> param; param.clear();
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	//~ dbi.WriteImage(imInBlue,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen,"imInGreen.png");
	//~ dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<uchar> imInGreen_p(w,h);
	DipImage<uchar> imInRed_p(w,h);
	DipImage<uchar> imInBlue_p(w,h);
	ImColorBandSeparation(frame_previous,imInBlue_p,imInGreen_p,imInRed_p);
	//~ dbi.WriteImage(imInBlue_p,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen_p,"imInGreen.png");
	//~ dbi.WriteImage(imInRed_p,"imInRed.png");
	
	DipImage<uchar> imInGreen_0(w,h);
	DipImage<uchar> imInRed_0(w,h);
	DipImage<uchar> imInBlue_0(w,h);
	ImColorBandSeparation(frame_0,imInBlue_0,imInGreen_0,imInRed_0);
	
	DipImage<Tmask> imMaskInitialInteral(w,h);imMaskInitialInteral.init(0);
	param.clear();
	param.push_back(int(w/100));
	ImErode(imMaskInitial,Fast,param,imMaskInitialInteral);
	dbi.WriteImage(imMaskInitialInteral,"imMaskInitialInteral.png");
	
	DipImage<Tmask> imFGTmp_0(w,h);imFGTmp_0.init(0);
	DipImage<Tmask> imFGEdge_0(w,h);imFGEdge_0.init(0);
	param.clear();
	param.push_back(int(5));
	ImErode(imMaskInitialInteral,Fast,param,imFGTmp_0);
	ImSubImage(imMaskInitialInteral,Clip,imFGTmp_0,imFGEdge_0);
	dbi.WriteImage(imFGEdge_0,"imFGEdge_0.png");
	
	
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	
	param.push_back(int(w/50));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");
	
	DipImage<Tmask> imBGTmp(w,h);imBGTmp.init(0);
	DipImage<Tmask> imBGEdge(w,h);imBGEdge.init(0);
	param.clear();
	param.push_back(int(2));
	ImDilate(imMaskExteral,Fast,param,imBGTmp);
	ImSubImage(imBGTmp,Clip,imMaskExteral,imBGEdge);
	dbi.WriteImage(imBGEdge,"imBGEdge.png");
	
	param.clear();
	param.push_back(int(w/50));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	
	DipImage<Tmask> imFGTmp(w,h);imFGTmp.init(0);
	DipImage<Tmask> imFGEdge(w,h);imFGEdge.init(0);
	param.clear();
	param.push_back(int(2));
	ImErode(imMaskInteral,Fast,param,imFGTmp);
	ImSubImage(imMaskInteral,Clip,imFGTmp,imFGEdge);
	dbi.WriteImage(imFGEdge,"imFGEdge.png");

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				//~ std::cout << "initial"<<endl;
				//~ std::cout << i << "   " << j <<endl;
				int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				borneInfx = std::max(1,i-w/20);
				borneInfy = std::max(1,j-w/20);
				borneMaxx = std::min(w, i+w/20+1);
				borneMaxy = std::min(h, j+w/20+1);
				std::vector<int> edgeFGStore;edgeFGStore.clear();
				std::vector<int> edgeBGStore;edgeBGStore.clear();
				std::vector<int> edgeFG_0Store;edgeFG_0Store.clear();
				
				for (int ii=borneInfx; ii<borneMaxx; ii++){
					for (int jj=borneInfy; jj<borneMaxy; jj++){
						if(imFGEdge(ii,jj) != 0){
							edgeFGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));	
						}
						if(imBGEdge(ii,jj) != 0){
							edgeBGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));	
						}
						if(imFGEdge_0(ii,jj) != 0){
							edgeFG_0Store.push_back(int(abs(imInGreen_0(ii,jj)- imInGreen(i,j)) * abs(imInGreen_0(ii,jj)- imInGreen(i,j)) + abs(imInGreen_0(ii,jj)- imInBlue(i,j)) * abs(imInGreen_0(ii,jj)- imInBlue(i,j)) + abs(imInGreen_0(ii,jj)- imInRed(i,j)) *abs(imInGreen_0(ii,jj)- imInRed(i,j))));	
						}
					}
				}
				
				//~ std::vector<double>::iterator it;
				//~ std::cout << "edgeFGStore =" << edgeFGStore.size()<<endl;
				//~ std::cout << "edgeBGStore =" << edgeBGStore.size()<<endl;
				//~ std::cout << "edgeFG_0Store =" << edgeFG_0Store.size()<<endl;
				//~ cin.get();
				if(edgeFGStore.size() <10){imFGAdd(i,j) = 0;}
				else if(edgeBGStore.size() <10){imFGAdd(i,j) = 255;}
				else if(edgeFG_0Store.size() <10){
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());
					
					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					while(count !=10){						
						double edgeFG = edgeFGStore[0];
						double edgeBG = edgeBGStore[0];
						if(edgeFG <edgeBG){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
						}
						else if(edgeBG <edgeFG){
							edgeBGcount  = edgeBGcount +1;
							edgeBGStore.erase(edgeBGStore.begin());
						}
						
						count = count +1;
					}
					if(edgeFGcount > edgeBGcount)
						imFGAdd(i,j) = 255;
					
					
				}
				else{
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());
					std::sort(edgeFG_0Store.begin(),edgeFG_0Store.end());

					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					int edgeFG_0count = 0;
					while(count !=10){
						int edgeFG = edgeFGStore[0];
						int edgeBG = edgeBGStore[0];
						int edgeFG_0 = edgeFG_0Store[0];
						//~ std::cout << "edgeFG =" << edgeFG<<endl;
						//~ std::cout << "edgeBG =" << edgeBG<<endl;
						//~ std::cout << "edgeFG_0 =" << edgeFG_0<<endl;
						if((edgeFG <=edgeBG) && (edgeFG <= edgeFG_0) ){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
						}
						else if((edgeBG <edgeFG) && (edgeBG < edgeFG_0) ){
							edgeBGcount  = edgeBGcount +1;
							edgeBGStore.erase(edgeBGStore.begin());
						}
						else if((edgeFG_0 <edgeFG) && (edgeFG_0 <= edgeBG) ){
							edgeFG_0count  = edgeFG_0count +1;
							edgeFG_0Store.erase(edgeFG_0Store.begin());
						}
						else {
							std::cout << "wrong" <<endl;
						}

						
						count = count +1;
					}
					//~ std::cout << "edgeFGcount =" << edgeFGcount<<endl;
					//~ std::cout << "edgeBGcount =" << edgeBGcount<<endl;
					//~ std::cout << "edgeFG_0count =" << edgeFG_0count<<endl;
					//~ std::cout << i << "   " << j <<endl;
					if((edgeFGcount + edgeFG_0count) > edgeBGcount){
						imFGAdd(i,j) = 255;
					}
					
				}					
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");
	//~ for (int i= 0; i<w; i++){
		//~ for (int j=0; j<h; j++){
			//~ if(imMaskUnknown(i,j) !=0){
				//~ int bgdistance = 0;
				//~ int fgdistance = 0;
				//~ int bgCount = 0;
				//~ int bgEffectCount = 0;
				//~ int fgCount = 0;
				//~ int fgEffectCount = 0;
				//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				//~ borneInfx = std::max(0,i-w/30);
				//~ borneInfy = std::max(0,j-w/30);
				//~ borneMaxx = std::min(w, i+w/30+1);
				//~ borneMaxy = std::min(h, j+w/30+1);
				//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
					//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
						//~ if(imBGEdge(ii,jj) !=0){
							//~ bgCount = bgCount +1;
							//~ bgdistance = bgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ bgEffectCount = bgEffectCount +1;
							
						//~ }
						//~ if(imFGEdge(ii,jj) !=0){
							//~ fgdistance = fgdistance + abs(imInGreen(ii,jj)- imInGreen(i,j)) * abs(imInGreen(ii,jj)- imInGreen(i,j));
							//~ fgCount = fgCount + 1;
							//~ if(abs(imInGreen(ii,jj)- imInGreen(i,j)) < 10 && abs(imInRed(ii,jj)- imInRed(i,j)) < 10 && abs(imInBlue(ii,jj)- imInBlue(i,j)) < 10){
								//~ fgEffectCount = fgEffectCount +1;
							
						//~ }
					//~ }
				//~ }
				//~ if(double(fgdistance)/fgCount < double(bgdistance)/bgCount){
					//~ imFG(i,j) = 255;
				//~ }
				//~ else{
					//~ imBG(i,j) = 255;
				//~ }
			//~ }
		//~ }
	//~ }
	//~ 
	//~ dbi.WriteImage(imFG,"imFG_fii.png");
	
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ DipImage<Tmask> imProcessingPoint(w,h);imProcessingPoint.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp(w,h);imMaskInteralTemp.init(0);
	//~ DipImage<Tmask> imMaskInteralTemp2(w,h);imMaskInteralTemp2.init(0);
	//~ std::vector<ParamT> param; param.clear();
	//~ param.push_back(int(1));
	//~ ImDilate(imMaskInteral,Fast,param,imMaskInteralTemp2);
	//~ ImCopy(imMaskInteral,imMaskInteralTemp);
	//~ ImSubImage(imMaskInteralTemp2,imMaskTotal,Clip,imMaskInteralTemp,imProcessingPoint);
	//~ dbi.WriteImage(imProcessingPoint,"imProcessingPoint.png");
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){			
			//~ if(imProcessingPoint(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	
	
	
	//~ DipImage<Tin> imMaskUnknownInvert(w,h); imMaskUnknownInvert.init(0);
	//~ ImInvert(imMaskUnknown,imMaskUnknownInvert);
	//~ dbi.WriteImage(imMaskUnknownInvert, "imMaskUnknownInvert.png");
	//~ 
	//~ DipImage<Tin> imDistanceVN(w,h); imDistanceVN.init(0);
	//~ ImMorphoDistance(imMaskUnknownInvert,imDistanceVN);
	//~ dbi.WriteImage(imDistanceVN, "imDistanceVN.png");

	//~ DipImage<Tin> imSkeleton(w,h); imSkeleton.init(0);
	//~ ImThinning_WithMask(imDistanceVN, imMaskUnknown, imSkeleton);
	//~ dbi.WriteImage(imSkeleton, "imSkeleton.png");
	//~ 
	//~ DipImage<Tin> imSkeletonCleanSecond(w,h); imSkeletonCleanSecond.init(0);
	//~ ImSkeletonCleaningBranch_WithMask(imSkeleton,imMaskUnknown,imMaskUnknown,imSkeletonCleanSecond);
	//~ dbi.WriteImage(imSkeletonCleanSecond, "imSkeletonCleanSecond.png");
	
	//~ 
	//~ DipImage<uchar> imEdge(w,h);
	//~ MaskedgeDetection(frame,imMaskUnknown, imEdge);
	//~ 
	//~ std::vector<std::pair<int,int> > pointProcessing;pointProcessing.clear();
	//~ for(int i = 0;i < w;i++){
		//~ for(int j = 0 ;j< h; j++){
			//~ //if(imSkeleton(i,j) == 255){
			//~ if(imEdge(i,j) == 255){
				//~ pointProcessing.push_back(std::make_pair(i,j));
			//~ }
		//~ }
	//~ }
	//~ 
	//~ 
	//~ DipImage<int> imPatialBG(w,h); imPatialBG.init(0);
	//~ int count;
	//~ for(int i = 0; i < pointProcessing.size();i++){
		//~ count = count + 1;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/30);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/30);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/30+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/30+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ 
			//~ }
		//~ }
		
		//~ //kmeans algorithm
		//~ Mat src(frame, Rect(borneInfx,borneInfy,borneMaxx - borneInfx,borneMaxy - borneInfy ));
		//~ Mat samples(src.rows * src.cols, 3, CV_32F);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
				//~ for( int z = 0; z < 3; z++)
					//~ samples.at<float>(y + x*src.rows, z) = src.at<Vec3b>(y,x)[z];
//~ 
//~ 
		//~ int clusterCount = 2;
		//~ Mat labels;
		//~ int attempts = 5;
		//~ Mat centers;
		//~ kmeans(samples, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers );
//~ 
//~ 
		//~ DipColorImage<uchar> new_image(src.cols,src.rows);
		//~ for( int y = 0; y < src.rows; y++ )
			//~ for( int x = 0; x < src.cols; x++ )
			//~ { 
			//~ int cluster_idx = labels.at<int>(y + x*src.rows,0);
			//~ new_image(x,y)[0] = centers.at<float>(cluster_idx, 0);
			//~ new_image(x,y)[1] = centers.at<float>(cluster_idx, 1);
			//~ new_image(x,y)[2] = centers.at<float>(cluster_idx, 2);
		//~ }
		
		//~ dbi.WriteImage(new_image, boost::lexical_cast<string>(count) + "new_image.png");
		

		
		
		
		//~ std::vector<std::vector<int> > data;
		//~ std::vector<std::vector<int> > results;
		//~ std::vector<int> colorInfo;colorInfo.clear();
		//~ std::vector<int> sizeG;
		//~ int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
		//~ borneInfx = std::max(0,pointProcessing[i].first-w/35);
		//~ borneInfy = std::max(0,pointProcessing[i].second-w/35);
		//~ borneMaxx = std::min(w, pointProcessing[i].first+w/35+1);
		//~ borneMaxy = std::min(h, pointProcessing[i].second+w/35+1);
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){
				//~ colorInfo.clear();
				//~ colorInfo.push_back(imInGreen(ii,jj));
				//~ data.push_back(colorInfo);
			//~ }
		//~ }
		//~ results = KMeans(data,2,sizeG);
		//~ dbi.WriteInfo("results 1 %d",results[0][0]);
		//~ dbi.WriteInfo("results 2 %d",results[1][0]);
		//~ bool isPointBG;
		//~ if(borneInfx < w/2)isPointBG = false;
		//~ else isPointBG = true;
		//~ for (int ii=borneInfx; ii<borneMaxx; ii++){
			//~ for (int jj=borneInfy; jj<borneMaxy; jj++){				
				//~ if(abs(imInGreen(ii,jj)-results[0][0]) > abs(imInGreen(ii,jj)-results[0][0])
			//~ }
		//~ }
		
		
	//~ }
	
	
	dbi.WriteOut("Leaving imagematting\n");
}



template<typename Tin, typename Tmask> void imageMatting_4(DipColorImage<Tin> & frame,DipColorImage<Tin> & frame_previous,DipImage<Tmask> & imMask,DipImage<Tmask> & imMaskOut){
	dbi.WriteEnter(" imageMatting_4\n");
	const int w = frame.width(), h = frame.height();
	
	//rgb images
	DipImage<uchar> imInGreen(w,h);
	DipImage<uchar> imInRed(w,h);
	DipImage<uchar> imInBlue(w,h);
	ImColorBandSeparation(frame,imInBlue,imInGreen,imInRed);
	//~ dbi.WriteImage(imInBlue,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen,"imInGreen.png");
	//~ dbi.WriteImage(imInRed,"imInRed.png");
	
	DipImage<uchar> imInGreen_p(w,h);
	DipImage<uchar> imInRed_p(w,h);
	DipImage<uchar> imInBlue_p(w,h);
	ImColorBandSeparation(frame_previous,imInBlue_p,imInGreen_p,imInRed_p);
	//~ dbi.WriteImage(imInBlue_p,"imInBlue.png");
	//~ dbi.WriteImage(imInGreen_p,"imInGreen.png");
	//~ dbi.WriteImage(imInRed_p,"imInRed.png");
	
	
	DipImage<Tmask> imMaskInteral(w,h);imMaskInteral.init(0);
	DipImage<Tmask> imMaskExteral(w,h);imMaskExteral.init(0);
	DipImage<Tmask> imMaskTotal(w,h);imMaskTotal.init(255);
	std::vector<ParamT> param; param.clear();
	param.push_back(int(w/30));
	ImDilate(imMask,Fast,param,imMaskExteral);
	dbi.WriteImage(imMaskExteral,"imMaskExteral.png");

	param.clear();
	param.push_back(int(w/50));
	ImErode(imMask,Fast,param,imMaskInteral);
	dbi.WriteImage(imMaskInteral,"imMaskInteral.png");
	

	
	DipImage<Tmask> imMaskUnknown(w,h);imMaskUnknown.init(0);
	ImSubImage(imMaskExteral,NoClip,imMaskInteral,imMaskUnknown);
	dbi.WriteImage(imMaskUnknown,"imMaskUnknown.png");
	
	DipImage<Tmask> imBGMask(w,h);imBGMask.init(0);
	ImInvert(imMaskExteral,imBGMask);
	dbi.WriteImage(imBGMask,"imBGMask.png");
	
	
	DipImage<Tmask> imFG(w,h);imFG.init(0);
	DipImage<Tmask> imBG(w,h);imBG.init(0);
	DipImage<Tmask> imFGAdd(w,h);imFGAdd.init(0);
	ImCopy(imMaskInteral,imFG);
	
	for (int i= 0; i<w; i++){
		for (int j=0; j<h; j++){
			if(imMaskUnknown(i,j) != 0){
				//~ std::cout << "initial"<<endl;
				//~ std::cout << i << "   " << j <<endl;
				int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
				borneInfx = std::max(1,i-w/10);
				borneInfy = std::max(1,j-w/10);
				borneMaxx = std::min(w, i+w/10+1);
				borneMaxy = std::min(h, j+w/10+1);
				//~ std::cout << "borneInfx = "<< borneInfx << endl;
				//~ std::cout << "borneInfy = "<< borneInfy << endl;
				//~ std::cout << "borneMaxx = "<< borneMaxx << endl;
				//~ std::cout << "borneMaxy = "<< borneMaxy << endl;
				std::vector<int> edgeFGStore;edgeFGStore.clear();
				std::vector<int> edgeBGStore;edgeBGStore.clear();
				
				int count1 = 0;
				while(edgeFGStore.size()<50 || edgeBGStore.size()<50){
					if(count1>500){break;}
					int ii = borneInfx + (rand()%1000) * (borneMaxx-borneInfx)/1000;
					int jj = borneInfy + (rand()%1000) * (borneMaxy-borneInfy)/1000;
					//~ std::cout << "ii =" << ii<< "jj =" <<jj <<endl;
					if(imMaskInteral(ii,jj) !=0){
						edgeFGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));						
					}
					if(imBGMask(ii,jj) !=0){
						edgeBGStore.push_back(int(abs(imInGreen_p(ii,jj)- imInGreen(i,j)) * abs(imInGreen_p(ii,jj)- imInGreen(i,j)) + abs(imInBlue_p(ii,jj)- imInBlue(i,j)) * abs(imInBlue_p(ii,jj)- imInBlue(i,j)) + abs(imInRed_p(ii,jj)- imInRed(i,j)) *abs(imInRed_p(ii,jj)- imInRed(i,j))));
							
					}
					count1 = count1 + 1;					
				}
						
				//~ std::vector<double>::iterator it;
				//~ std::cout << "edgeFGStore =" << edgeFGStore.size()<<endl;
				//~ std::cout << "edgeBGStore =" << edgeBGStore.size()<<endl;
				

				if(edgeFGStore.size() <40){imFGAdd(i,j) = 0;}
				else if(edgeBGStore.size() <40){imFGAdd(i,j) = 255;}
				else{
					std::sort(edgeFGStore.begin(),edgeFGStore.end());
					std::sort(edgeBGStore.begin(),edgeBGStore.end());

					int count = 0;
					int edgeFGcount = 0;
					int edgeBGcount = 0;
					while(count !=30){
						//~ it = edgeFGStore.begin();
						double edgeFG = edgeFGStore[0];
						//~ double edgeFG = *it;
						//~ it = edgeBGStore.begin();
						//~ double edgeBG = *it;
						double edgeBG = edgeBGStore[0];
						if(edgeFG <edgeBG){
							edgeFGcount  = edgeFGcount +1;
							edgeFGStore.erase(edgeFGStore.begin());
							//~ edgeFGStore.pop_front();
						}
						if(edgeBG <edgeFG){
							edgeBGcount  = edgeBGcount +1;
							//~ edgeBGStore.pop_front();
							edgeBGStore.erase(edgeBGStore.begin());
						}
						
						count = count +1;
					}
					if(edgeFGcount > edgeBGcount)
						imFGAdd(i,j) = 255;
				}
								
			}
		}
	}
	
	ImAddImage(imFG,Clip,imFGAdd,imMaskOut);
	dbi.WriteImage(imMaskOut,"imMaskOut.png");

	
	dbi.WriteOut("Leaving imagematting\n");
}



#endif
