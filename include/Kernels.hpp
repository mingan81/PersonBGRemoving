/*! @file Filters.hpp
 *  @brief Filters toolbox 
 */
#ifndef KERNELS_HPP
#define KERNELS_HPP

using namespace std;
/**
 * @class Kernel
 */
template<typename T> class Kernel: public std::list<std::pair<std::pair<int,int>, T> > {};


/*!
 * Neighborhood: is a list of coordinates that is used as a structuring element for morphological processing 
 */
 template<typename T> class Neighborhood: public std::list<std::pair<T,T> >
 {

 };



/**
 *
 * @param radius
 * @param val
 * @param ker
 */
template<typename Tkernel> void DiscKernelGenerator(const int radius, const Tkernel val, Kernel<Tkernel> & ker)
{
	dbi.WriteEnter("Entering DiscKernelGenerator\n");
	ker.clear();
	int j = 0, i = 0, jj=0;
	const int sqradius = radius * radius;
	for (j = -radius; j <= radius; j++){
		jj = j*j;
		for (i = -radius; i <= radius; i++){
			if ((i*i + jj) <= sqradius){
				ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), val));
			}
		}
	}
	dbi.WriteOut("Leaving DiscKernelGenerator\n");
}


/**
 * @fn template<typename Tkernel> void CircleKernelGenerator(const int radius, const Tkernel val, Kernel<Tkernel> & ker)
 * @brief Generate unfilled circle kernel  
 * @param radius : circle radius
 * @param val : value for circle
 * @param ker : output kernel
 */
template<typename Tkernel> void CircleKernelGenerator(const int radius, const Tkernel val, Kernel<Tkernel> & ker)
{
	dbi.WriteEnter("Entering CricleKernelGenerator\n");
	ker.clear();
	int j = 0, i = 0, jj=0;
	const int sqradius = radius * radius;
	for (j = -radius; j <= radius; j++){
		jj = j*j;
		for (i = -radius; i <= radius; i++){
			if ((double(i*i + jj) >= (double(sqradius-radius)+0.25)) && (double(i*i + jj) <= (double(sqradius+radius)+0.25))){
				ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), val));
			}
			//else{
			//	ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), Tkernel(0)));
			//}
		}
	}
	dbi.WriteInfo("Kernel size = %d", ker.size());
	dbi.WriteOut("Leaving CircleKernelGenerator\n");
}


/**
 * @fn template<typename Tkernel> void TwoCirclesKernelGenerator(const int radius, const Tkernel val, Kernel<Tkernel> & ker)
 * @brief Generate unfilled circle kernel  
 * @param radius : circle radius
 * @param val : value for circle
 * @param ker : output kernel
 */
template<typename Tkernel> void TwoCirclesKernelGenerator(const int radiusInternal,const int radiusExternal, const Tkernel valInternal,const Tkernel valExternal, Kernel<Tkernel> & ker)
{
	dbi.WriteEnter("Entering CricleKernelGenerator\n");
	ker.clear();
	int j = 0, i = 0, jj=0;
	int sqradius = radiusInternal * radiusInternal;
	for (j = -radiusInternal; j <= radiusInternal; j++){
		jj = j*j;
		for (i = -radiusInternal; i <= radiusInternal; i++){
			if ((double(i*i + jj) >= (double(sqradius-radiusInternal)+0.25)) && (double(i*i + jj) <= (double(sqradius+radiusInternal)+0.25))){
				ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), valInternal));
			}
		}
	}
	sqradius = radiusExternal * radiusExternal;
	for (j = -radiusExternal; j <= radiusExternal; j++){
		jj = j*j;
		for (i = -radiusExternal; i <= radiusExternal; i++){
			if ((double(i*i + jj) >= (double(sqradius-radiusExternal)+0.25)) && (double(i*i + jj) <= (double(sqradius+radiusExternal)+0.25))){
				ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), valExternal));
			}
		}
	}
	dbi.WriteInfo("Kernel size = %d", ker.size());
	dbi.WriteOut("Leaving CircleKernelGenerator\n");
}



/**
 * @fn template<typename Tkernel> void CircleKernelGenerator(const int radiusmin, const int radiusmax, const Tkernel val, Kernel<Tkernel> & ker)
 * @brief Generate unfilled circle kernel  
 * @param radius : circle radius
 * @param radius : circle radius
 * @param val : value for circle
 * @param ker : output kernel
 */ 
template<typename Tkernel> void CircleKernelGenerator(const int radiusmin, const int radiusmax, const Tkernel val, Kernel<Tkernel> & ker)
{
	dbi.WriteEnter("Entering CircleKernelGenerator\n");
	ker.clear();
	int j = 0, i = 0, jj=0;
	for (j = -radiusmax; j <= radiusmax; j++){
		jj = j*j;
		for (i = -radiusmax; i <= radiusmax; i++){
			if ((double(i*i + jj) >= double(radiusmin*radiusmin)) && (double(i*i + jj) <= double(radiusmax*radiusmax)))
			{
				ker.push_back(std::make_pair(std::make_pair<int,int>(i,j), val));
			}
		}
	}
	dbi.WriteInfo("Kernel size = %d", ker.size());
	dbi.WriteOut("Leaving CircleKernelGenerator\n");
}


/*! @fn template<typename Tsigma, typename Tkernel> void GaussianKernelGenerator_2D(const int radx, const int rady, const int mux, const int muy, const Tsigma sigmax, const Tsigma sigmay, const double corr, const Tkernel ampl, Kernel<Tkernel> & kernel)
 * @brief 2-Dimention Gaussian Kernel Generation 
 * @param radx
 * @param rady
 * @param mux
 * @param muy
 * @param sigmax
 * @param sigmay
 * @param corr
 * @param ampl
 * @param kernel :output 
 */
template<typename Tsigma, typename Tkernel> void GaussianKernelGenerator_2D(const int radx, const int rady, const int mux, const int muy, const Tsigma sigmax, const Tsigma sigmay, const double corr, const Tkernel ampl, Kernel<Tkernel> & kernel)
{
  dbi.WriteEnter("Entering GaussianKernelGenerator\n");
  int i=0, j=0;
  kernel.clear();
  //const double K = 1.0 / double(sigmax*sigmay * 2 * PI * sqrt(1 - corr*corr));
  double valtemp = 0.0;
  double mx = 0.0;
  double my = 0.0;
  double diffx = 0.0;
  double diffy = 0.0;
  double mcorr = 0.0;
  double ccorr = corr;
  if (ccorr == 1.0)
	  ccorr = 0.9999999999999;
  for (j=-rady; j<=rady; j++)
  {
    for (i=-radx; i<=radx; i++)
    {
	diffx = double(i) - double(mux);
	diffy = double(j) - double(muy);
	mx = diffx * diffx / (double(sigmax) * double(sigmax));
	my = diffy * diffy / (double(sigmay) * double(sigmay));
	mcorr = 2 * ccorr * diffx * diffy / (sigmax * sigmay);
	valtemp = ampl/*K*/ * exp(-(1.0 / double(2*(1-ccorr*ccorr))) * (mx + my - mcorr)  );
	// valtemp = K * exp(-(1.0 / double(2*(1-corr*corr))) * (mx + my - mcorr)  );
	kernel.push_back(std::make_pair(std::make_pair<int,int>(i,j),Tkernel(valtemp)));
    }
  }
  dbi.WriteOut("Leaving GaussianKernelGenerato\n");
}




/*! @fn template<typename Tkernel> void KernelNormalisation(const Kernel<Tkernel> & kin, Kernel<double> & kout)
 * @brief normalization for kernel 
 * @param kin
 * @param kout
 */
template<typename Tkernel> void KernelNormalisation(const Kernel<Tkernel> & kin, Kernel<double> & kout)
{
	dbi.WriteEnter("Entering KernelNormalisation\n");
	typename Kernel<Tkernel>::const_iterator itk = kin.begin(), itkend = kin.end();
	double somme = 0.0;
	for (; itk != itkend; itk++){
		// dbi.WriteInfo("valeur = %f\n", (*itk).second);
		somme += (*itk).second;
	}

	// dbi.WriteInfo("Taille de la liste d'entree = %d\n", kin.size());
	// dbi.WriteInfo("somme = %f\n", somme);
	kout.clear();
	itk = kin.begin();
	for (; itk != itkend; itk++){
		kout.push_back(std::make_pair<std::pair<int,int>, double>((*itk).first,double((*itk).second) / somme));
	}

	// dbi.WriteInfo("Taille de la liste de sortie = %d\n", kout.size());

	dbi.WriteOut("Leaving KernelNormalisation\n");
}


/**
 *
 * @return
 */
template<typename Tkernel> static const Kernel<Tkernel> SobelXKernelGenerator_2D()
{
	Kernel<Tkernel> ktmp; ktmp.clear();
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,-1),-1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,-1),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,-1),1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,0),-2));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,0),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,0),2));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,1),-1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,1),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,1),1));
	return(ktmp);
}

/**
 * SobelXKernel_2D
 */
static const Kernel<int> SobelXKernel_2D = SobelXKernelGenerator_2D<int>();

/**
 *
 * @return
 */
template<typename Tkernel> static const Kernel<Tkernel> SobelYKernelGenerator_2D()
{
	Kernel<Tkernel> ktmp; ktmp.clear();
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,-1),-1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,-1),-2));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,-1),-1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,0),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,0),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,0),0));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,1),1));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,1),2));
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,1),1));
	return(ktmp);
}

/**
 * SobelYKernel_2D
 */
static const Kernel<int> SobelYKernel_2D = SobelYKernelGenerator_2D<int>();



/**
 *
 * @return
 */
template<typename Tkernel> static const Kernel<Tkernel> GradientXKernelGenerator_2D()
{
	Kernel<Tkernel> ktmp; ktmp.clear();	
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(-1,0),-1));	
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(1,0),1));
		
	return(ktmp);
}

/**
 * GradientXKernel_2D
 */
static const Kernel<int> GradientXKernel_2D = GradientXKernelGenerator_2D<int>();

/**
 *
 * @return
 */
template<typename Tkernel> static const Kernel<Tkernel> GradientYKernelGenerator_2D()
{
	Kernel<Tkernel> ktmp; ktmp.clear();
	
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,-1),-1));	
	ktmp.push_back(make_pair<pair<int,int>, Tkernel>(make_pair<int,int>(0,1),1));	
	return(ktmp);
}

/**
 * GradientYKernel_2D
 */
static const Kernel<int> GradientYKernel_2D = GradientYKernelGenerator_2D<int>();




/*! @fn void CreateNeighborhood(const NeighborhoodShape shape, const std::vector<ParamT> param, Neighborhood<T> & liste)
 * @brief Create a structuring element with user-defined shape and size
 * @param shape : shape of the neighborhood (Disc, Segment, Rect, Cross, Square)
 * @param param : vector of parameters (either a single radius for Disc, Square and Cross, or two radius for Rect and Cross, or a radius and an angle for Segment)
 * @param liste : resulting Neighborhood
 */
template<typename T> void CreateNeighborhood(const NeighborhoodShape shape, const std::vector<ParamT> param, Neighborhood <T> & liste)
{
	//~ dbi.WriteEnter("Entering CreateNeighborhood \n");
	int Nparam = param.size();
	liste.clear();
	int i(0), j(0); 
	int radius(0), radiusx(0), radiusy(0);
	double angle(0);
	double radius_segm(0);
	
	switch (shape)
	{
		case Disc:
		{
			if(Nparam!=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if(param[0].type() != typeid(int)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			radius = boost::get<int>(param.front());	
			
			int rayint = int(std::abs(radius) + 0.5);
			int raysq = radius * radius;
			for (i=-rayint; i<=rayint; i++)
			{
				for (j=-rayint; j<=rayint; j++)
				{
					if ((i*i + j*j) <= raysq)
					liste.push_back(std::pair<T,T>(i,j));
				}
			}
			break;
		}
		case DIP_Rect : 
		{
			if(Nparam!=2){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if((param[0].type() != typeid(int)) || (param[1].type()==typeid(int))){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			radiusx = boost::get<int>(param.front());					
			radiusy = boost::get<int>(param.back());
							
			int rayintx = int(std::abs(radiusx) + 0.5);
			int rayinty = int(std::abs(radiusy) + 0.5);
			for (i=-rayintx; i<=rayintx; i++)
			{
				for (j=-rayinty; j<=rayinty; j++)
				{
				  liste.push_back(std::pair<T,T>(i,j));
				}
			}
			break;
		}
		case Cross:
		{
			if((Nparam!=1) && (Nparam!=2)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
						
			if(Nparam==2) 
			{
				if((param[0].type() != typeid(int)) || (param[1].type()==typeid(int))){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				
				radiusx = boost::get<int>(param.front());					
				radiusy = boost::get<int>(param.back());					
				
				const int rayintx = int(std::abs(radiusx) + 0.5);
				const int rayinty = int(std::abs(radiusy) + 0.5);
			  
				liste.push_back(std::pair<T,T>(0,0));
				
				for (i=1; i<=rayinty; i++)
				{
					liste.push_back(std::pair<T,T>(0,i));
					liste.push_back(std::pair<T,T>(0,-i));
				}
				for (i=1; i<=rayintx; i++)
				{
					liste.push_back(std::pair<T,T>(i,0));
					liste.push_back(std::pair<T,T>(-i,0));
				}
			}
			else if(Nparam==1) 
			{
				if((param[0].type() != typeid(int))){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				
				radiusx = boost::get<int>(param.front());					
				radiusy = boost::get<int>(param.front());					
				
				const int rayintx = int(std::abs(radiusx) + 0.5);
				const int rayinty = int(std::abs(radiusy) + 0.5);
			  
				liste.push_back(std::pair<T,T>(0,0));
				
				for (i=1; i<=rayinty; i++)
				{
					liste.push_back(std::pair<T,T>(0,i));
					liste.push_back(std::pair<T,T>(0,-i));
				}
				for (i=1; i<=rayintx; i++)
				{
					liste.push_back(std::pair<T,T>(i,0));
					liste.push_back(std::pair<T,T>(-i,0));
				}
			}
			break;
		}
		case Square:
		{
			if((Nparam!=1)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if((param[0].type() != typeid(int))){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			radius = boost::get<int>(param.front());					
				
			int rayintx = int(std::abs(radius) + 0.5);
			int rayinty = int(std::abs(radius) + 0.5);
			for (i=-rayintx; i<=rayintx; i++)
			{
				for (j=-rayinty; j<=rayinty; j++)
				{
				  liste.push_back(std::pair<T,T>(i,j));
				}
			}
			break;
		}
		case Segment:
		{			
			if((Nparam!=2)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if((param[0].type() != typeid(double)) || (param[1].type()!=typeid(double))){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			radius_segm = boost::get<double>(param.front());
			angle = boost::get<double>(param.back());
			
			double alpha = (angle) * PI / 180;
				  
			if ((abs(angle) <= 45) || (abs(angle) >= 135))
			{   
				double cosangle = cos(alpha);
				int rayintx = int(radius_segm * cosangle + sign(cosangle) * 0.5);
				int bi=std::min(0, rayintx), bf=std::max(0, rayintx);
				for (i=bi; i<=bf; i++)
				{
					j=int(i*tan(alpha)+sign(alpha)*0.5);
					liste.push_back(std::pair<T,T>(i,j));
					if (i != 0)
					liste.push_back(std::pair<T,T>(-i,-j));
				}
			}
			else
			{   
				double sinangle = sin(alpha);
				double tanangle = tan(alpha);
				int rayinty = int(radius_segm * sinangle + sign(sinangle) * 0.5);
				int bi=std::min(0, rayinty), bf=std::max(0, rayinty);
				for (j=bi; j<=bf; j++)
				{
					i = int(j/tanangle + sign(j)*sign(tanangle) * 0.5);
					liste.push_back(std::pair<T,T>(i,j));
					if (j != 0)
					liste.push_back(std::pair<T,T>(-i,-j));
				}
			}			
			break;
		}
		
	}
	//~ dbi.WriteOut("Leaving CreateNeighborhood.\n");
};


/**
 *
 * @param ray
 * @param liste
 */
template<typename Tray, typename Tlist> void Circle2DNListGenerator(Tray ray, Neighborhood<Tlist> & liste)
{
  dbi.WriteEnter("Entering Circular2DNListGenerator\n");
  int i=0, j=0;
  int rayint = int(std::abs(ray) + 0.5);
  int raysq = ray * ray;
  liste.clear();
  for (i=-rayint; i<=rayint; i++)
  {
    for (j=-rayint; j<=rayint; j++)
    {
      if ((i*i + j*j) == raysq)
	liste.push_back(std::pair<Tlist,Tlist>(i,j));
    }
  }
  dbi.WriteOut("Leaving Circular2DNListGenerator\n");
}

/**
 *
 * @param ray
 * @param liste
 */
template<typename Tray, typename Tlist> void Disc2DNListGenerator(Tray ray, Neighborhood<Tlist> & liste)
{
 // dbi.WriteEnter("Entering Disc2DNListGenerator\n");
  int i=0, j=0;
  int rayint = int(std::abs(ray) + 0.5);
  int raysq = ray * ray;
  liste.clear();
  for (i=-rayint; i<=rayint; i++)
  {
    for (j=-rayint; j<=rayint; j++)
    {
      if ((i*i + j*j) <= raysq)
	liste.push_back(std::pair<Tlist,Tlist>(i,j));
    }
  }
  // dbi.WriteOut("Leaving Disc2DNListGenerator\n");
}

/**
 *
 * @param rayx
 * @param rayy
 * @param liste
 */
template<typename Tray, typename Tlist> void Rectangle2DNListGenerator(Tray rayx, Tray rayy, Neighborhood<Tlist> & liste)
{
  int i=0, j=0;
  liste.clear();
  int rayintx = int(std::abs(rayx) + 0.5);
  int rayinty = int(std::abs(rayy) + 0.5);
  for (i=-rayintx; i<=rayintx; i++)
  {
    for (j=-rayinty; j<=rayinty; j++)
    {
      liste.push_back(std::pair<int,int>(i,j));
    }
  }
}

/**
 *
 * @param rayx
 * @param rayy
 * @param liste
 */
template<typename Tray, typename Tlist> void RectangleCross2DNListGenerator(const Tray rayx, const Tray rayy, Neighborhood<Tlist> & liste)
{
  int i=0;
  liste.clear();
  const int rayintx = int(std::abs(rayx) + 0.5);
  const int rayinty = int(std::abs(rayy) + 0.5);
  
  liste.push_back(std::pair<int,int>(0,0));
  
  for (i=1; i<=rayinty; i++)
  {
    liste.push_back(std::pair<int,int>(0,i));
    liste.push_back(std::pair<int,int>(0,-i));
  }
  for (i=1; i<=rayintx; i++)
  {
    liste.push_back(std::pair<int,int>(i,0));
    liste.push_back(std::pair<int,int>(-i,0));
  }
}

/**
 *
 * @param radius
 * @param liste
 */
template<typename Tray, typename Tlist> void SquareCross2DNListGenerator(const Tray radius, Neighborhood<Tlist> & liste)
{
  RectangleCross2DNListGenerator(radius, radius, liste);
}

/**
 *
 * @param ray
 * @param angle
 * @param liste
 */
template<typename Tray, typename Tlist> void Linear2DNListGenerator_old(double ray, double angle, std::list<std::pair<int,int> > & liste)
{
  double alpha = (-angle) * PI / 180;
  int i=0, j=0;
  
  if ((abs(angle) <= 45) || (abs(angle) >= 135))
  {   
    double cosangle = cos(alpha);
    int rayintx = int(ray * cosangle + sign(cosangle) * 0.5);
    int bi=std::min(0, rayintx), bf=std::max(0, rayintx);
    for (i=bi; i<=bf; i++)
    {
	j=int(i*tan(alpha)+sign(alpha)*0.5);
	liste.push_back(std::pair<int,int>(i,j));
	if (i != 0)
	  liste.push_back(std::pair<int,int>(-i,-j));
    }
  }
  else
  {   
    double sinangle = sin(alpha);
    double tanangle = tan(alpha);
    int rayinty = int(ray * sinangle + sign(sinangle) * 0.5);
    int bi=std::min(0, rayinty), bf=std::max(0, rayinty);
    for (j=bi; j<=bf; j++)
    {
	i = int(j/tanangle + sign(j)*sign(tanangle) * 0.5);
	liste.push_back(std::pair<int,int>(i,j));
	if (j != 0)
	  liste.push_back(std::pair<int,int>(-i,-j));
    }
  }
}

/**
 *
 * @param ray
 * @param angle
 * @param liste
 */
template<typename Trad, typename Tlist> void Linear2DNListGenerator(Trad ray, double angle, std::list<std::pair<int,int> > & liste)
{
  double alpha = (angle) * PI / 180;
  int i=0, j=0;
  
  if ((abs(angle) <= 45) || (abs(angle) >= 135))
  {   
    double cosangle = cos(alpha);
    int rayintx = int(ray * cosangle + sign(cosangle) * 0.5);
    int bi=std::min(0, rayintx), bf=std::max(0, rayintx);
    for (i=bi; i<=bf; i++)
    {
	j=int(i*tan(alpha)+sign(alpha)*0.5);
	liste.push_back(std::pair<int,int>(i,j));
	if (i != 0)
	  liste.push_back(std::pair<int,int>(-i,-j));
    }
  }
  else
  {   
    double sinangle = sin(alpha);
    double tanangle = tan(alpha);
    int rayinty = int(ray * sinangle + sign(sinangle) * 0.5);
    int bi=std::min(0, rayinty), bf=std::max(0, rayinty);
    for (j=bi; j<=bf; j++)
    {
	i = int(j/tanangle + sign(j)*sign(tanangle) * 0.5);
	liste.push_back(std::pair<int,int>(i,j));
	if (j != 0)
	  liste.push_back(std::pair<int,int>(-i,-j));
    }
  }
}

/**
 *
 * @param ray
 * @param liste
 */
template<typename Tray, typename Tlist> void Square2DNListGenerator(Tray ray, Neighborhood<Tlist> & liste)
{
  Rectangle2DNListGenerator<Tray, Tlist>(ray, ray, liste);
}

/**
 *
 * @param nh
 * @param imOut
 */
template<typename Tlist, typename Tout> void NeighborhoodToImage(const Neighborhood<Tlist> nh, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering NeighborhoodToImage\n");
	typename Neighborhood<Tlist>::const_iterator it = nh.begin(), itend = nh.end();
	const int imax = std::numeric_limits<int>::max();
	int xmini = imax, ymini = imax;
	int xmaxi = -imax+1, ymaxi = -imax+1;

	for (; it != itend; it++){
		if (it->first > xmaxi) xmaxi = it->first;
		else if (it->first < xmini) xmini = it->first;
		
		if (it->second > ymaxi) ymaxi = it->second;
		else if (it->second < ymini) ymini = it->second;
	}
	const int w = xmaxi-xmini+1;
	const int h = ymaxi-ymini+1;
	dbi.WriteInfo("(xmini,xmaxi) = (%d,%d)\n", xmini,xmaxi);
	dbi.WriteInfo("(ymini,ymaxi) = (%d,%d)\n", ymini,ymaxi);
	dbi.WriteInfo("(w,h) = (%d,%d)\n", w,h);
	imOut.resize(w,h);
	imOut.init(0);

	it = nh.begin();

	Tout outmax = std::numeric_limits<Tout>::max();
	for (; it != itend; it++){
		imOut(it->first - xmini, it->second - ymini) = outmax;
	}
	dbi.WriteOut("Leaving NeighborhoodToImage\n");
}




#endif
