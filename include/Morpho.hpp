/*! @file Morpho.hpp
 *  @brief Morpho toolbox 
 */
#ifndef MORPHO_HPP
#define MORPHO_HPP

#include <omp.h>




/** 
 * @brief Isolate a 3x3 square neighborhood in the input image, into a smaller output image
 * @param imIn : input image
 * @param p_pos : pointer to current pixel (center of the neighborhood)
 * @param imOut : output image of size 3x3
 */
template<typename T> void CopyNeighborhoodInImage(const DipImage<T> & imIn, const T * p_pos, DipImage<T> & imOut)
{
	const int w = imIn.width();
	const int h = imIn.height();
	imOut.resize(3,3);
	const T * p_instart = &imIn(0,0);
	const T * p_intmp = &imIn(0,0);

	const int offset = p_pos - p_instart;

	const int reste = (offset+1) % w;

	T * p_out = &imOut(0,0);
	T * p_outtmp = &imOut(0,0);
	p_out+=4;

	int deltax=0, deltay=0;
	Neighborhood<int> nh; nh.clear();
	std::vector<ParamT> param; param.clear();
	param.push_back(int(1));
	CreateNeighborhood(Square, param, nh);
	Neighborhood<int>::iterator it = nh.begin(), itend = nh.end();

	if ((offset > (w-1)) && (offset < w*(h-1)) && (reste > 0) && (reste < w)){
		for (;it!=itend; it++){
			deltax = (*it).first;
			deltay = (*it).second;
			p_outtmp = p_out + (deltax + 3 * deltay);
			p_intmp = p_pos + (deltax + w * deltay);
			*p_outtmp = *p_intmp;
		}
	}
}

/*! 
 * @brief Computes morphological erosion with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T, typename Tmask> void ImErode(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImErode with a Mask\n");
	int Nparam = param.size();
	int w = imIn.width(), h = imIn.height();
	dbi.WriteInfo("parameter size = %d\n", Nparam);
	switch (method)
	{
		case Basic:
		{
			if(Nparam!=1){
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if(param[0].type() != typeid(std::list<std::pair<int, int> >)){
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			std::list<std::pair<int,int> > nliste = boost::get<std::list<std::pair<int,int> > > (param.front());
			  
			const T * p_center = &(imIn(0,0));
			const T * p_debutLigne = p_center;
			const T * p_finLigne = p_center + (w-1);
			T * p_out = &(imOut(0,0));
			const T * p_begin = p_center;
			const T * p_temp = p_center;
			const T * p_end = &(imIn(w-1,h-1));
			const T * p_neighbor = p_center;
			const T * p_minimum = p_center;
			const Tmask * p_maskbeg = & imMask(0,0);
			const Tmask * p_mask = p_maskbeg;
			const Tmask * p_masktmp = p_maskbeg;
			T maxi = std::numeric_limits<T>::max();
			std::list<std::pair<int,int> >::const_iterator it = nliste.begin(), itend = nliste.end();
			  
			// Debut de l erosion
			for(;p_center != p_end; p_center++, p_out++, p_mask++){
				if (*p_mask != 0){
					for(;it!=itend; it++){
						p_temp = p_center + (*it).first;
						if ((p_temp <= p_finLigne) && (p_temp >= p_debutLigne)){
							p_neighbor = p_temp + w * (*it).second;
							if ((p_neighbor <= p_end) && (p_neighbor >= p_begin)){
								p_masktmp = p_maskbeg + (p_neighbor - p_begin);
								if (*p_masktmp != 0){
									if (*p_minimum > *p_neighbor){
										p_minimum = p_neighbor;
									}
								}
							}
						}
					}
					*p_out = *p_minimum;
				}
				if (p_center > p_finLigne){
					p_debutLigne = p_debutLigne + w;
					p_finLigne = p_finLigne + w;
				}
				it = nliste.begin();
				p_minimum = &maxi;
			}
			
			break;
		}
		case Fast:
		{			
			bool OK = false;
			int rayonx, rayony;
			
			if((Nparam!=1) && (Nparam!=2)){
				dbi.WriteInfo("1\n");
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			
			if (Nparam==1){
				if(param[0].type() != typeid(int)){
					dbi.WriteInfo("2\n");
					//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param.front());
					rayony = rayonx;
					OK = true;
				}
			}
			else if (Nparam==2){
				if(param[0].type() != typeid(int) || param[1].type() != typeid(int)){
					dbi.WriteInfo("3\n");
					//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param[0]);
					rayony = boost::get<int> (param[1]);
					OK = true;
				}
			}	
			if (OK)
			{			
				typedef T TypeTmp;
					
				int i=0, j=0, ii=0, jj=0;
				int borneInfx = 0, borneMaxx = 0;
				
				TypeTmp minligne = std::numeric_limits<TypeTmp>::max();
				TypeTmp minglobal = std::numeric_limits<TypeTmp>::max();
				TypeTmp valtemp = 0;
				
				std::list<TypeTmp> minlist; if(!minlist.empty()) {printf("Attention: la liste n'est pas vide\n"); minlist.clear();}
				
				typename std::list<TypeTmp>::iterator itminlist = minlist.begin();
				typename std::list<TypeTmp>::iterator itminlistend = minlist.end();
				
				for (i=0; i<w; i++){
					/// Premier point de la colonne
					j=0;
					minlist.clear();
					minglobal = std::numeric_limits<TypeTmp>::max();
					minligne = std::numeric_limits<TypeTmp>::max();
					if (i<=rayonx){
						borneInfx = 0;
						borneMaxx = i+rayonx+1;
					}
					else if ((w-i)<=rayonx){
						borneInfx = i-rayonx;
						borneMaxx = w;
					}
					else{
						borneInfx = i-rayonx;
						borneMaxx = i+rayonx+1;
					}

					for (jj=0; jj<rayony+1; jj++){
						minligne = std::numeric_limits<TypeTmp>::max();
						for (ii= borneInfx; ii< borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
							}
						}
						if (minligne < minglobal) minglobal = minligne;
						minlist.push_back(minligne);
					}
					if (imMask(ii,jj) != 0){
						imOut(i,j) = minglobal;
					}
					/// Fin du traitement du premier point de la colonne
					
					/// Iteration sur la colonne
					for (j=1; j<h; j++){
						if (j <= rayony){
							// Ajout d'une nouvelle ligne
							jj = j+rayony;
							minligne = std::numeric_limits<TypeTmp>::max();
							for (ii=borneInfx; ii<borneMaxx; ii++){
								if (imMask(ii,jj) != 0){
									if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
								}
							}
							if (minligne < minglobal) minglobal = minligne;
							minlist.push_back(minligne);
						}
						else if ((h-j)<=rayony)
						{
							// On retire une ligne
							valtemp = minlist.front();
							minlist.pop_front();
							if (valtemp == minglobal){
								minglobal = std::numeric_limits<TypeTmp>::max();
								itminlist = minlist.begin();
								itminlistend = minlist.end();
								for(;itminlist != itminlistend; itminlist++)
									if (*itminlist < minglobal)
										minglobal = *itminlist;
							}
						}
						else
						{
								// On ajoute et on retire une ligne
							jj = j+rayony;
							minligne = std::numeric_limits<TypeTmp>::max();
							for (ii=borneInfx; ii<borneMaxx; ii++){
								if (imMask(ii,jj) != 0){
									if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
								}
							}
							if (minligne < minglobal) minglobal = minligne;
							minlist.push_back(minligne);

							valtemp = minlist.front();
							minlist.pop_front();
							if (valtemp == minglobal){
								minglobal = std::numeric_limits<TypeTmp>::max();
								itminlist = minlist.begin();
								itminlistend = minlist.end();
								for(;itminlist != itminlistend; itminlist++)
									if (*itminlist < minglobal)
										minglobal = *itminlist;
							}
						}
						if (imMask(i,j) != 0){
							imOut(i,j) = minglobal;
						}
					}
				}
			}
			break;
		}
	}
	
	dbi.WriteInfo("Leaving ImErode with a Mask. \n");
	
}

/*! @fn void ImErode_MT(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
 * @brief Computes morphological erosion with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T, typename Tmask> void ImErode_MT(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImErode with a Mask MT\n");
	int Nparam = param.size();
	int w = imIn.width(), h = imIn.height();
	dbi.WriteInfo("parameter size = %d\n", Nparam);
	switch (method)
	{
		case Basic:
		{
			if(Nparam!=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if(param[0].type() != typeid(std::list<std::pair<int, int> >)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}

			std::list<std::pair<int,int> > nliste = boost::get<std::list<std::pair<int,int> > > (param.front());

			const T * p_center = &(imIn(0,0));
			const T * p_debutLigne = p_center;
			const T * p_finLigne = p_center + (w-1);
			T * p_out = &(imOut(0,0));
			const T * p_begin = p_center;
			const T * p_temp = p_center;
			const T * p_end = &(imIn(w-1,h-1));
			const T * p_neighbor = p_center;
			const T * p_minimum = p_center;
			const Tmask * p_maskbeg = & imMask(0,0);
			const Tmask * p_mask = p_maskbeg;
			const Tmask * p_masktmp = p_maskbeg;
			T maxi = std::numeric_limits<T>::max();
			std::list<std::pair<int,int> >::const_iterator it = nliste.begin(), itend = nliste.end();

			// Debut de l erosion
			for(;p_center != p_end; p_center++, p_out++, p_mask++){
				if (*p_mask != 0){
					for(;it!=itend; it++){
						p_temp = p_center + (*it).first;
						if ((p_temp <= p_finLigne) && (p_temp >= p_debutLigne)){
							p_neighbor = p_temp + w * (*it).second;
							if ((p_neighbor <= p_end) && (p_neighbor >= p_begin)){
								p_masktmp = p_maskbeg + (p_neighbor - p_begin);
								if (*p_masktmp != 0){
									if (*p_minimum > *p_neighbor){
										p_minimum = p_neighbor;
									}
								}
							}
						}
					}
					*p_out = *p_minimum;
				}
				if (p_center > p_finLigne){
					p_debutLigne = p_debutLigne + w;
					p_finLigne = p_finLigne + w;
				}
				it = nliste.begin();
				p_minimum = &maxi;
			}

			break;
		}
		case Fast:
		{
			bool OK = false;
			int rayonx, rayony;

			if((Nparam!=1) && (Nparam!=2)){
				dbi.WriteInfo("1\n");
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}

			if (Nparam==1){
				if(param[0].type() != typeid(int)){
					dbi.WriteInfo("2\n");
					//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param.front());
					rayony = rayonx;
					OK = true;
				}
			}
			else if (Nparam==2){
				if(param[0].type() != typeid(int) || param[1].type() != typeid(int)){
					dbi.WriteInfo("3\n");
					//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param[0]);
					rayony = boost::get<int> (param[1]);
					OK = true;
				}
			}
			if (OK)
			{
				typedef T TypeTmp;

				int i=0, j=0, ii=0, jj=0;
				int borneInfx = 0, borneMaxx = 0;

				#pragma omp parallel private( i, j, ii, jj, borneInfx, borneMaxx ) firstprivate( w, h, rayonx, rayony ) shared( imIn, imMask, imOut )
				{
					TypeTmp minligne = std::numeric_limits<TypeTmp>::max();
					TypeTmp minglobal = std::numeric_limits<TypeTmp>::max();
					TypeTmp valtemp = 0;

					std::list<TypeTmp> minlist;

					typename std::list<TypeTmp>::iterator itminlist = minlist.begin();
					typename std::list<TypeTmp>::iterator itminlistend = minlist.end();

					#pragma omp for schedule( guided ) nowait
					for (i=0; i<w; i++){
						/// Premier point de la colonne
						j=0;
						minlist.clear();
						minglobal = std::numeric_limits<TypeTmp>::max();
						minligne = std::numeric_limits<TypeTmp>::max();
						if (i<=rayonx){
							borneInfx = 0;
							borneMaxx = i+rayonx+1;
						}
						else if ((w-i)<=rayonx){
							borneInfx = i-rayonx;
							borneMaxx = w;
						}
						else{
							borneInfx = i-rayonx;
							borneMaxx = i+rayonx+1;
						}

						for (jj=0; jj<rayony+1; jj++){
							minligne = std::numeric_limits<TypeTmp>::max();
							for (ii= borneInfx; ii< borneMaxx; ii++){
								if (imMask(ii,jj) != 0){
									if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
								}
							}
							if (minligne < minglobal) minglobal = minligne;
							minlist.push_back(minligne);
						}
						if (imMask(ii,jj) != 0){
							imOut(i,j) = minglobal;
						}
						/// Fin du traitement du premier point de la colonne

						/// Iteration sur la colonne
						for (j=1; j<h; j++){
							if (j <= rayony){
								// Ajout d'une nouvelle ligne
								jj = j+rayony;
								minligne = std::numeric_limits<TypeTmp>::max();
								for (ii=borneInfx; ii<borneMaxx; ii++){
									if (imMask(ii,jj) != 0){
										if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
									}
								}
								if (minligne < minglobal) minglobal = minligne;
								minlist.push_back(minligne);
							}
							else if ((h-j)<=rayony)
							{
								// On retire une ligne
								valtemp = minlist.front();
								minlist.pop_front();
								if (valtemp == minglobal){
									minglobal = std::numeric_limits<TypeTmp>::max();
									itminlist = minlist.begin();
									itminlistend = minlist.end();
									for(;itminlist != itminlistend; itminlist++)
										if (*itminlist < minglobal)
											minglobal = *itminlist;
								}
							}
							else
							{
									// On ajoute et on retire une ligne
								jj = j+rayony;
								minligne = std::numeric_limits<TypeTmp>::max();
								for (ii=borneInfx; ii<borneMaxx; ii++){
									if (imMask(ii,jj) != 0){
										if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
									}
								}
								if (minligne < minglobal) minglobal = minligne;
								minlist.push_back(minligne);

								valtemp = minlist.front();
								minlist.pop_front();
								if (valtemp == minglobal){
									minglobal = std::numeric_limits<TypeTmp>::max();
									itminlist = minlist.begin();
									itminlistend = minlist.end();
									for(;itminlist != itminlistend; itminlist++)
										if (*itminlist < minglobal)
											minglobal = *itminlist;
								}
							}
							if (imMask(i,j) != 0){
								imOut(i,j) = minglobal;
							}
						}
					}
				}
			}
			break;
		}
	}

	dbi.WriteInfo("Leaving ImErode with a Mask MT \n");

}

#ifdef _GPU_
/*! @fn void ImErode_GPU(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
 * @brief Computes morphological erosion with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T, typename Tmask> void ImErode_GPU( const DipGPUImage<T> & imInBI,
                                                       const DipGPUImage<Tmask> & imMaskBI,
                                                       const MorphoMethod method,
                                                       const std::vector<ParamT> & param,
                                                       DipGPUImage<T> & imOutBI )
{
	dbi.WriteInfo("Entering ImErode_GPU\n");
	int Nparam = param.size();
	dbi.WriteInfo("parameter size = %d\n", Nparam);
	switch (method)
	{
		case Basic:
		{
			throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("Not implemented yet");

			break;
		}
		case Fast:
		{
			bool OK = false;
			int rayonx, rayony;

			if((Nparam!=1) && (Nparam!=2)){
				dbi.WriteInfo("1\n");
				throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}

			if (Nparam==1){
				if(param[0].type() != typeid(int)){
					dbi.WriteInfo("2\n");
					//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param.front());
					rayony = rayonx;
					OK = true;
				}
			}
			else if (Nparam==2){
				if(param[0].type() != typeid(int) || param[1].type() != typeid(int)){
					dbi.WriteInfo("3\n");
					//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param[0]);
					rayony = boost::get<int> (param[1]);
					OK = true;
				}
			}

			if (OK)
			{
				const T MinValueStart = std::numeric_limits<T>::max( );

				const int width = imInBI.Width( );
				const int height = imInBI.Height( );

				// Do some padding to remove coslty tests
				const unsigned int GroupSize = std::max( DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeX_H" ),
														 DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeY_V" ) );
				const int NewHeight = Compute::Utils::GetNextMultipleOf( height, (int)GroupSize ) + rayony * 2 + GroupSize;
				const int NewWidth = Compute::Utils::GetNextMultipleOf( width, (int)GroupSize ) + rayonx * 2 + GroupSize;
				const int NewNbPixel = NewHeight * NewWidth;

				const cl::Context Context = Compute::CCompute::GetInstance( )->GetContext( );
				cl::CommandQueue Queue = Compute::CCompute::GetInstance()->GetCommandQueue( );
				cl_int Error;
				cl_mem_flags InOutMemFlags = CL_MEM_READ_WRITE;

				// Watch out for MaxWorkGroupSize
				unsigned int LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeY_H" );
				unsigned int LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeX_H" );
				Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
				// Global size must be a multiple of local size
				unsigned int GlobalSizeY = (height + LocalSizeY - 1) / LocalSizeY;
				GlobalSizeY *= LocalSizeY;
				unsigned int GlobalSizeX = (width + LocalSizeX - 1) / LocalSizeX;
				GlobalSizeX *= LocalSizeX;

				// We need new buffers to get our padding
				cl::Buffer ImNewInBuffer( Context,
										  InOutMemFlags, NewNbPixel * sizeof( T ),
										  0, &Error );
				cl::Buffer ImNewMaskBuffer( Context,
											InOutMemFlags, NewNbPixel * sizeof( unsigned char ),
											0, &Error );

				cl::Buffer ImMinBuffer( Context,
										InOutMemFlags, NewNbPixel * sizeof( T ),
										0, &Error );
				cl::Buffer ImMinBuffer_Transp( Context,
											   InOutMemFlags, NewNbPixel * sizeof( T ),
											   0, &Error );

				cl::Buffer ImNewOutBuffer( Context,
										   InOutMemFlags, NewNbPixel * sizeof( T ),
										   0, &Error );

				cl::Kernel KH = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastErode_WithMask_H_Pad" );

				cl::Kernel KV = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastErode_WithMask_V_Pad_Transpose" );

				cl::Kernel KT = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastErode_TempTranspose" );

				cl::Kernel KZero = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastErode_ZeroOutPaddingBuffers" );

				// Zero out padded buffers
				int Index = 0;
				KZero.setArg( Index++, ImNewInBuffer );
				KZero.setArg( Index++, ImNewMaskBuffer );
				KZero.setArg( Index++, ImMinBuffer );
				KZero.setArg( Index++, ImNewOutBuffer );
				KZero.setArg( Index++, MinValueStart );
				Queue.enqueueNDRangeKernel( KZero, cl::NullRange, cl::NDRange( NewNbPixel ), cl::NullRange, 0, 0 );

				// Copy the image and mask into their new buffers
				cl::size_t<3> SrcOrigin;
				SrcOrigin.push_back( 0 );
				SrcOrigin.push_back( 0 );
				SrcOrigin.push_back( 0 );
				cl::size_t<3> DstOrigin;
				DstOrigin.push_back( rayonx * sizeof( T ) );
				DstOrigin.push_back( rayony );
				DstOrigin.push_back( 0 );
				cl::size_t<3> Region;
				Region.push_back( width * sizeof( T ) );
				Region.push_back( height );
				Region.push_back( 1 );
				size_t SrcRowPitch = width * sizeof( T );
				size_t SrcSlicePitch = height * SrcRowPitch;
				size_t DstRowPitch = NewWidth * sizeof( T );
				size_t DstSlicePitch = NewHeight * DstRowPitch;
				Queue.enqueueCopyBufferRect( imInBI.Buffer( ), ImNewInBuffer, SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );
				DstOrigin[ 0 ] = rayonx * sizeof( unsigned char );
				Region[ 0 ] = width * sizeof( unsigned char );
				SrcRowPitch = width * sizeof( unsigned char );
				SrcSlicePitch = height * SrcRowPitch;
				DstRowPitch = NewWidth * sizeof( unsigned char );
				DstSlicePitch = NewHeight * DstRowPitch;
				Queue.enqueueCopyBufferRect( imMaskBI.Buffer( ), ImNewMaskBuffer, SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );

				// Horizontal pass
				Index = 0;
				KH.setArg( Index++, ImNewInBuffer );
				KH.setArg( Index++, ImNewMaskBuffer );
				KH.setArg( Index++, ImMinBuffer );
				KH.setArg( Index++, NewHeight );
				KH.setArg( Index++, NewWidth );
				KH.setArg( Index++, rayonx );
				KH.setArg( Index++, rayony );
				KH.setArg( Index++, MinValueStart );
				KH.setArg( Index++, LocalSizeX * 2 * sizeof( T ), 0 );
				KH.setArg( Index++, LocalSizeX * 2 * sizeof( unsigned char ), 0 );

				Queue.enqueueNDRangeKernel( KH, cl::NullRange, cl::NDRange( GlobalSizeX, GlobalSizeY ), cl::NDRange( LocalSizeX, LocalSizeY ), 0, 0 );

				// TESTING
				DipGPUImage<T> TmpMinBI( ImMinBuffer, NewWidth, NewHeight );
				DipImage<T> TmpMin( NewWidth, NewHeight ); TmpMin.init( 0 );
				TmpMinBI.ToImage( TmpMin );
				dbi.WriteImage( TmpMin, "TmpMin.png" );

				// Transpose
				// Watch out for MaxWorkGroupSize
				LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeY_T" );
				LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeX_T" );
				Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
				// Global size must be a multiple of local size
				GlobalSizeY = (NewHeight + LocalSizeY - 1) / LocalSizeY;
				GlobalSizeY *= LocalSizeY;
				GlobalSizeX = (NewWidth + LocalSizeX - 1) / LocalSizeX;
				GlobalSizeX *= LocalSizeX;

				Index = 0;
				KT.setArg( Index++, ImMinBuffer );
				KT.setArg( Index++, ImMinBuffer_Transp );
				KT.setArg( Index++, NewHeight );
				KT.setArg( Index++, NewWidth );
				KT.setArg( Index++, LocalSizeY * LocalSizeX * sizeof( T ), 0 );

				Queue.enqueueNDRangeKernel( KT, cl::NullRange, cl::NDRange( GlobalSizeX, GlobalSizeY ), cl::NDRange( LocalSizeX, LocalSizeY ), 0, 0 );

				// Vertical pass
				// Watch out for MaxWorkGroupSize
				LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeY_V" );
				LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImErode_Fast_WithMask_GPU.LocalSizeX_V" );
				Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
				// Global size must be a multiple of local size
				GlobalSizeY = (height + LocalSizeY - 1) / LocalSizeY;
				GlobalSizeY *= LocalSizeY;
				GlobalSizeX = (width + LocalSizeX - 1) / LocalSizeX;
				GlobalSizeX *= LocalSizeX;

				Index = 0;
				KV.setArg( Index++, ImMinBuffer_Transp );
				KV.setArg( Index++, ImNewMaskBuffer );
				KV.setArg( Index++, ImNewOutBuffer );
				KV.setArg( Index++, NewHeight );
				KV.setArg( Index++, NewWidth );
				KV.setArg( Index++, rayonx );
				KV.setArg( Index++, rayony );
				KV.setArg( Index++, LocalSizeY * 2 * sizeof( T ), 0 );

				Queue.enqueueNDRangeKernel( KV, cl::NullRange, cl::NDRange( GlobalSizeY, GlobalSizeX ), cl::NDRange( LocalSizeY, LocalSizeX ), 0, 0 );

				// Copy the relevant sub rect into the output
				SrcOrigin[ 0 ] = rayonx * sizeof( T );
				SrcOrigin[ 1 ] = rayony;
				SrcOrigin[ 2 ] = 0;
				DstOrigin[ 0 ] = 0;
				DstOrigin[ 1 ] = 0;
				DstOrigin[ 2 ] = 0;
				Region[ 0 ] = width * sizeof( T );
				SrcRowPitch = NewWidth * sizeof( T );
				SrcSlicePitch = NewHeight * SrcRowPitch;
				DstRowPitch = width * sizeof( T );
				DstSlicePitch = height * DstRowPitch;
				Queue.enqueueCopyBufferRect( ImNewOutBuffer, imOutBI.Buffer( ), SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );
			}
			break;
		}
	}

	dbi.WriteInfo("Leaving ImErode_GPU\n");

}

/*! @fn void ImErode_GPU(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
 * @brief Computes morphological erosion
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T> void ImErode_GPU(const DipGPUImage<T> & imInBI, const MorphoMethod method, const std::vector<ParamT> & param, DipGPUImage<T> & imOutBI)
{
	dbi.WriteInfo("Entering ImErode_GPU\n");
	DipGPUImage<unsigned char> imMaskBI(imInBI.Width(), imInBI.Height()); imMaskBI.init(255);

	ImErode_GPU(imInBI, imMaskBI, method, param, imOutBI);

	dbi.WriteInfo("Leaving ImErode_GPU\n");
}

/*! @fn void ImErode_GPU(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
 * @brief Computes morphological erosion
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T> void ImErode_GPU(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImErode_GPU\n");
	DipGPUImage<unsigned char> imMaskBI(imIn.width(), imIn.height()); imMaskBI.init(255);
	DipGPUImage<T> imInBI( imIn, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR );
	DipGPUImage<T> imOutBI( imOut, CL_MEM_WRITE_ONLY );
	ImErode_GPU(imInBI, imMaskBI, method, param, imOutBI);

	imOutBI.ToImage( imOut );

	dbi.WriteInfo("Leaving ImErode_GPU\n");
}
#endif // _GPU_

/*! 
 * @brief Computes morphological erosion 
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T> void ImErode(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImErode \n");
	DipImage<unsigned char> imMask(imIn.width(), imIn.height()); imMask.init(255);
	ImErode(imIn, imMask, method, param, imOut);
	
	dbi.WriteInfo("Leaving ImErode. \n");
}

/*! @fn void ImErode_MT(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
 * @brief Computes morphological erosion
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : eroded image
 */
template<typename T> void ImErode_MT(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImErode_MT\n");
	DipImage<unsigned char> imMask(imIn.width(), imIn.height()); imMask.init(255);
	ImErode_MT(imIn, imMask, method, param, imOut);

	dbi.WriteInfo("Leaving ImErode_MT\n");
}

 
/*! 
 * @brief Computes morphological dilatation with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : dilated image
 */
template<typename T, typename Tmask> void ImDilate(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImDilate with a Mask\n");
	int Nparam = param.size();
	int w = imIn.width(), h = imIn.height();
	
	switch (method)
	{
		case Basic:
		{
			if(Nparam!=1){
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if(param[0].type() != typeid(std::list<std::pair<int, int> >)){
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			std::list<std::pair<int,int> > nliste = boost::get<std::list<std::pair<int,int> > > (param.front());
				  
			const T * p_begin = &(imIn(0,0));
			const T * p_debutLigne = p_begin;
			const T * p_finLigne = p_begin + (w-1);
			const T * p_center = p_begin;
			const T * p_temp = p_begin;
			const T * p_end = p_begin + w*h;
			const T * p_neighbor = p_begin;
			const T * p_maximum = p_begin;
			T * p_out = &(imOut(0,0));
			const Tmask * p_maskbeg = & imMask(0,0);
			const Tmask * p_mask = p_maskbeg;
			const Tmask * p_masktmp = p_maskbeg;
			T mini = std::numeric_limits<T>::min();
			std::list<std::pair<int,int> >::const_iterator it = nliste.begin(), itend = nliste.end();
		  
		  
			// Debut de la diltation
			for(;p_center != p_end; p_center++, p_out++, p_mask++){
				if (*p_mask != 0){
					for(;it!=itend; it++){
						p_temp = p_center + (*it).first;
						if ((p_temp <= p_finLigne) && (p_temp >= p_debutLigne)){
							p_neighbor = p_temp + w * (*it).second;
							if ((p_neighbor <= p_end) && (p_neighbor >= p_begin)){
								p_masktmp = p_maskbeg + (p_neighbor - p_begin);
								if (*p_masktmp != 0)
									if (*p_maximum < *p_neighbor)
										p_maximum = p_neighbor;
							}
						}
					}
					*p_out = *p_maximum;
				}
				if (p_center > p_finLigne){
					p_debutLigne = p_debutLigne + w;
					p_finLigne = p_finLigne + w;
				}
				it = nliste.begin();
				p_maximum = &mini;
			}
			break;
		}
		case Fast:
		{
			bool OK = false;
			int rayonx, rayony;
			
			if((Nparam!=1) && (Nparam!=2)){
				//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			
			if (Nparam==1){
				if(param[0].type() != typeid(int)){
					//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param.front());
					rayony = rayonx;
					OK = true;
				}
			}
			else if (Nparam==2){
				if(param[0].type() != typeid(int) || param[1].type() != typeid(int)){
					//throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
				}
				else
				{
					rayonx = boost::get<int> (param[0]);
					rayony = boost::get<int> (param[1]);
					OK = true;
				}
			}	
			if (OK)
			{					
				typedef T TypeTmp;
					
				int i=0, j=0, ii=0, jj=0;
				int borneInfx = 0, borneMaxx = 0;
				
				TypeTmp maxligne = 0;
				TypeTmp maxglobal = 0;
				TypeTmp valtemp = 0;
				
				std::list<TypeTmp> maxlist; if(!maxlist.empty()) {printf("Attention: la liste n'est pas vide\n"); maxlist.clear();}
				
				typename std::list<TypeTmp>::iterator itmaxlist = maxlist.begin();
				typename std::list<TypeTmp>::iterator itmaxlistend = maxlist.end();
				
				for (i=0; i<w; i++){
					/// Premier point de la colone
					j=0;
					maxlist.clear();
					maxglobal = 0;
					maxligne = 0;
					if (i<=rayonx){
						borneInfx = 0;
						borneMaxx = i+rayonx+1;
					}
					else if ((w-i)<=rayonx){
						borneInfx = i-rayonx;
						borneMaxx = w;
					}
					else{
						borneInfx = i-rayonx;
						borneMaxx = i+rayonx+1;
					}

					for (jj=0; jj<rayony+1; jj++){
						maxligne = 0;
						for (ii= borneInfx; ii< borneMaxx; ii++){
							if (imMask(ii,jj) != 0){
								if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
							}
						}
						if (maxligne > maxglobal) maxglobal = maxligne;
						maxlist.push_back(maxligne);
					}
					if (imMask(i,j) != 0){
						imOut(i,j) = maxglobal;
					}
					/// Fin du traitement du premier point de la colone
					
					/// Iteration sur la colone
					for (j=1; j<h; j++){
						if (j<=rayony){
							// Ajout d'une nouvelle ligne
							jj = j+rayony;
							maxligne = 0;
							for (ii=borneInfx; ii<borneMaxx; ii++){
								if (imMask(ii,jj) != 0){
									if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
								}
							}
							if (maxligne > maxglobal) maxglobal = maxligne;
							maxlist.push_back(maxligne);
						}
						else if ((h-j)<=rayony)
						{
							// On retire une ligne
							valtemp = maxlist.front();
							maxlist.pop_front();
							if (valtemp == maxglobal)
							{
								maxglobal = 0;
								itmaxlist = maxlist.begin();
								itmaxlistend = maxlist.end();
								for(;itmaxlist != itmaxlistend; itmaxlist++)
									if (*itmaxlist > maxglobal)
										maxglobal = *itmaxlist;
							}
						}
						else
						{
							// On ajoute et on retire une ligne
							jj = j+rayony;
							maxligne = 0;
							for (ii=borneInfx; ii<borneMaxx; ii++){
								if (imMask(ii,jj) != 0){
									if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
								}
							}
							if (maxligne > maxglobal) maxglobal = maxligne;
							maxlist.push_back(maxligne);

							valtemp = maxlist.front();
							maxlist.pop_front();
							if (valtemp == maxglobal)
							{
								maxglobal = 0;
								itmaxlist = maxlist.begin();
								itmaxlistend = maxlist.end();
								for(;itmaxlist != itmaxlistend; itmaxlist++)
									if (*itmaxlist > maxglobal)
										maxglobal = *itmaxlist;
							}
						}
						if (imMask(i,j) != 0){
							imOut(i,j) = maxglobal;
						}
					}
				}
			}
			break;
		}
	}
	dbi.WriteInfo("Leaving ImDilate with a Mask. \n");
}

/*! 
 * @brief Computes morphological dilatation 
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : dilated image
 */
template<typename T> void ImDilate(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteInfo("Entering ImDilate \n");
	DipImage<unsigned char> imMask(imIn.width(), imIn.height()); imMask.init(255);
	ImDilate(imIn, imMask, method, param, imOut);
	
	dbi.WriteInfo("Leaving ImDilate. \n");
}

/*! 
 * @brief Computes the morphological opening of the input image with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : opened image
 */
template<typename T, typename Tmask> void ImOpen(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteEnter("Entering ImOpen with a Mask\n");
	DipImage<T> imErode(imIn.width(),imIn.height()); imErode.init(0);
	ImErode(imIn, imMask, method, param, imErode);
	ImDilate(imErode, imMask, method, param, imOut);
	dbi.WriteOut("Leaving ImOpen with a Mask.\n");
}

/*! 
 * @brief Computes the morphological opening of the input image 
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : opened image
 */
template<typename T> void ImOpen(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteEnter("Entering ImOpen\n");
	DipImage<T> imErode(imIn.width(),imIn.height()); imErode.init(0);
	ImErode(imIn, method, param, imErode);
	ImDilate(imErode, method, param, imOut);
	dbi.WriteOut("Leaving ImOpen\n");
}

/*! 
 * @brief Computes the morphological closing of the input image with respect to non null points of imMask
 * @param imIn : input image
 * @param imMask : mask image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : closed image
 */
template<typename T, typename Tmask> void ImClose(const DipImage<T> & imIn, const DipImage<Tmask> & imMask, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteEnter("Entering ImClose with a Mask\n");
	DipImage<T> imDilate(imIn.width(),imIn.height()); imDilate.init(0);
	ImDilate(imIn, imMask, method, param, imDilate);
	ImErode(imDilate, imMask, method, param, imOut);
	dbi.WriteOut("Leaving ImClose with a Mask\n");
}

/*! 
 * @brief Computes the morphological closing of the input image 
 * @param imIn : input image
 * @param method : Basic or Fast (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imOut : closed image
 */
template<typename T> void ImClose(const DipImage<T> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<T> & imOut)
{
	dbi.WriteEnter("Entering ImClose\n");
	DipImage<T> imDilate(imIn.width(),imIn.height()); imDilate.init(0);
	ImDilate(imIn, method, param, imDilate);
	ImErode(imDilate, method, param, imOut);
	dbi.WriteOut("Leaving ImClose.\n");
}

/*! 
 * @brief Computes the morphological distance on the object given by all the non zero points with regard to the mask image
 * @param imObject : image containing the object
 * @param imMask : mask image
 * @param imDistance : distance image
 */
template<typename Tin, typename Tmask, typename Tout> void ImMorphoDistance(const DipImage<Tin> & imObject, const DipImage<Tmask> & imMask, DipImage<Tout> & imDistance)
{
	dbi.WriteEnter("Entering ImMorphoDistance with a Mask\n");
	imDistance.init(0);
	Neighborhood<int> nh; nh.clear();
	std::vector<ParamT> param; param.clear();
	param.push_back(int(1));
	CreateNeighborhood(Square, param, nh);

	const int w = imObject.width(), h = imObject.height();
  
	const Tin * p_objectStart = &imObject(0,0);
	const Tin * p_object = p_objectStart;
	const Tin * p_objectEnd = p_objectStart + w * h -1;
	
	const Tmask * p_maskStart = &imMask(0,0);
	const Tmask * p_mask = p_maskStart;
	
	Tout * p_distanceStart = &imDistance(0,0);
	Tout * p_distance = p_distanceStart;
	
	DipImage<uchar> flag(w,h);
	uchar * p_flagStart = &flag(0,0);
	uchar * p_flag = p_flagStart;
	
	std::list<const Tin *> liste1, liste2;
	liste1.clear(); liste2.clear();
	
	Tout current_distance = Tout(0);
	for(;p_object != p_objectEnd; p_object++, p_mask++, p_flag++, p_distance++)
	{
		if ( (*p_object != 0) && (*p_mask != 0))
		{
		(*p_flag) = 255;
		liste1.push_back(p_object);
		(*p_distance) = current_distance;
		}
		else
		*p_flag = 0;
	}
	current_distance++;
	
	int y = 0;
	int offset = 0;
	Neighborhood<int>::iterator nhit = nh.begin(), nhend = nh.end();
	const Tin * p_debutligne = p_objectStart;
	const Tin * p_finligne = p_objectStart + w-1;
	const Tin * p_neighbor = p_objectStart;
	const Tin * p_temp = p_objectStart;
	
	int wbis = w-1;
	typename std::list<const Tin *>::iterator itliste = liste1.begin(), itlisteend=liste1.end();
	while(itliste != itlisteend)
	{
		p_object = *itliste;
		y = (p_object - p_objectStart) / w ;
		p_debutligne = p_objectStart + y * w;
		p_finligne = p_debutligne + wbis;
	
		nhit = nh.begin();
		for (; nhit != nhend; nhit++)
		{
			p_temp = p_object + (*nhit).first;
			if ( (p_temp >= p_debutligne) && (p_temp <= p_finligne) ){
				p_neighbor = p_temp + w * (*nhit).second;
				offset = p_neighbor - p_objectStart;
				p_flag = p_flagStart + offset;
				p_mask = p_maskStart + offset;
				if ( (p_neighbor >= p_objectStart) && (p_neighbor <= p_objectEnd) && (*p_flag == 0) && (*p_mask != 0)){
					p_distance = p_distanceStart + offset;
					*p_distance = Tout(current_distance);
					liste2.push_back(p_neighbor);
					*p_flag = 255;
				}
			}
		}
		itliste++;
		if (itliste == itlisteend)
		{
			liste1.clear();
			liste1.swap(liste2);
			current_distance++;
			itliste = liste1.begin();
			itlisteend = liste1.end();
		}
	}
	dbi.WriteOut("Leaving ImMorphoDistance with a Mask\n");
}

/*! 
 * @brief Computes the morphological distance on the object image
 * @param imObject : image containing the object
 * @param imDistance : distance image
 */
template<typename Tin, typename Tout> void ImMorphoDistance(const DipImage<Tin> & imObject, DipImage<Tout> & imDistance)
{
	dbi.WriteEnter("Entering ImMorphoDistance \n");
	DipImage<unsigned char> imMask(imObject.width(), imObject.height()); imMask.init(255);
	ImMorphoDistance(imObject, imMask, imDistance);
	dbi.WriteOut("Leaving ImMorphoDistance. \n");
}

/*! 
 * @brief This function performs a segmentation of an image (imIn) using a morphological flooding constrained by a mask image.
 * @param imIn : image to be flooded
 * @param imLabel : seed for flooding
 * @param imOut : flooded image
 */
template<typename Tin, typename Tlabel>
void ImFlooding(const DipImage<Tin> & imIn, const DipImage<Tlabel> & imLabel, DipImage<Tlabel> & imOut)
{
	dbi.WriteEnter("Entering ImFlooding\n");
	const int w = imLabel.width(), h = imLabel.height();  
	imOut.init(0);  
	const int rayon = 1;
	int maxi = std::numeric_limits<Tin>::max();
	dbi.WriteInfo("maxi = %d\n", int(maxi)+1);
  
  	std::vector< std::queue<Node<Tin, Tlabel> > > hq;
	hq.resize(maxi+1);
  
	int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
  
	// deque filling
	int i = 0, j = 0, ii = 0, jj = 0;
	int current_level = std::numeric_limits<int>::max(); // Cette valeur doit rester un int (sinon, il y a des problemes lorsqu'on atteind la valeur max)
	Tlabel current_label = 0;
	Tin current_value = 0;
	Node<Tin, Tlabel> current_node;
  
	int flag = false;
  
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			flag = false;
			if (imLabel(i,j) != 0)
			{
				imOut(i,j) = imLabel(i,j);
				borneInfx = std::max(0,i-rayon);
				borneInfy = std::max(0,j-rayon);
				borneMaxx = std::min(w, i+rayon+1);
				borneMaxy = std::min(h, j+rayon+1);
	
				for (ii=borneInfx; ii<borneMaxx; ii++)
					for (jj=borneInfy; jj<borneMaxy; jj++)
						if (imLabel(ii,jj) == 0)
							flag = true;
	
				if (flag==true){
					current_value = imIn(i,j);
					current_node.SetPos(i,j);
					current_node.SetValue(current_value);
					current_node.SetLabel(imLabel(i,j));
					hq[current_value].push(current_node);
					if (int(current_value) < current_level) 
						current_level = int(current_value);
				}
			}	
		}
	}
  
	for (int k=0; k<256; k++)
  
		// processing
		while (current_level < (int(maxi) + 1)){
			while (!hq[current_level].empty()){
				current_node = hq[current_level].front();
				hq[current_level].pop();
				i = current_node.GetXPos();
				j = current_node.GetYPos();
				current_label = current_node.GetLabel();
			  
				if (current_label == 0)
					std::cout<<" Problem with node label"<< std::endl;
			  
			  
				borneInfx = std::max(0,i-rayon);
				borneInfy = std::max(0,j-rayon);
				borneMaxx = std::min(w, i+rayon+1);
				borneMaxy = std::min(h, j+rayon+1);
			  
				for (ii=borneInfx; ii<borneMaxx; ii++){
					for (jj=borneInfy; jj<borneMaxy; jj++){
						if (imLabel(ii,jj) == 0){
							if(imOut(ii,jj) == 0){
								current_value = imIn(ii,jj);
								hq[current_value].push(Node<Tin,Tlabel>(ii,jj,current_value,current_label));
								imOut(ii,jj) = current_label;
								if (int(current_value) < current_level){
									current_level = int(current_value);
								}
							}
						}
					}
				}
			}
			current_level++;
		}
		
	dbi.WriteOut("Leaving ImFlooding\n");
}

/*! 
 * @brief Performs the constrained dilation of imMark under imIn
 * @param imMark : marker image
 * @param imIn : constraint image
 * @param imOut : result image
 */
template<typename T> void ImUnderBuild(const DipImage<T> & imMark, const DipImage<T> & imIn, DipImage<T> & imOut)
{
	dbi.WriteEnter("Entering ImUnderBuild\n");
	const int w = imMark.width(), h = imMark.height();
	DipImage<T> imFlag(w,h); imFlag.init(0);
	int i=0, j=0, ii=0, jj=0, ilocal=0, jlocal=0;
	int borneInfx = 0, borneInfy = 0, borneMaxx = 0, borneMaxy = 0;
	int valtemp = 0;
	int rayon = 1;
	imOut.init(0);
	unsigned int seuil = 10;

	// Creation de la file d'attente hierarchique
	std::vector<std::list<Node<int, unsigned int> > > fah;
	fah.resize(256);

	int current_level = 0;
	T valeurtemp = 0;
	Node<int, unsigned int> noeudtemp;

	// On empile tous les points du marqueur
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if (imMark(i,j) != 0){
				valeurtemp = std::min(imMark(i,j), imIn(i,j));
				imOut(i,j) = valeurtemp;
				Node<int, unsigned int> noeud(i, j, valeurtemp, 0);
				fah[valeurtemp].push_back(noeud);
				if (int(valeurtemp) > current_level) 
					current_level = int(valeurtemp);
			}
		}
	}
	dbi.WriteInfo("End of initialization.\n");

	while(current_level >= 0){
		while (!fah[current_level].empty()){
			noeudtemp = fah[current_level].front();
			fah[current_level].pop_front();
			ilocal = noeudtemp.GetXPos();
			jlocal = noeudtemp.GetYPos();
			valtemp = noeudtemp.GetValue();

			// On empile les points du voisinage
			borneInfx = std::max(0,ilocal - rayon);
			borneInfy = std::max(0,jlocal - rayon);
			borneMaxx = std::min(w, ilocal + rayon+1);
			borneMaxy = std::min(h, jlocal + rayon+1);

			for (ii=borneInfx; ii<borneMaxx; ii++){
				for (jj=borneInfy; jj<borneMaxy; jj++){
					if (imFlag(ii,jj) != seuil){
						if (valtemp > imOut(ii,jj)){
							valeurtemp = std::min(valtemp, int(imIn(ii,jj)));
							Node<int, unsigned int> noeud(ii, jj, valeurtemp, 0);
							fah[valeurtemp].push_back(noeud);
							imFlag(ii,jj) = seuil;
							imOut(ii,jj) = valeurtemp;
						}
					}
				}
			}
		}
		current_level--;
	}
	dbi.WriteOut("Leaving ImUnderBuild\n");
}

/*! 
 * @brief Performs the constrained dilation of imMask under imIn
 * @param imMark : marker image
 * @param imIn : constraint image
 * @param iter : iteration number
 * @param imOut : result image
 */
template<typename Tmark, typename Tin, typename Tout> void ImUnderBuild_WithIterCount(const DipImage<Tmark> & imMark, const DipImage<Tin> & imIn, const int iter, DipImage<Tout> & imOut)
{
  dbi.WriteEnter("Entering ImUnderBuild_WithIterCount\n");
  imOut.init(0);
  int w = imMark.width(), h = imMark.height();
  DipImage<int> imFlag(w,h); imFlag.init(0);
  int i=0, j=0, ii=0, jj=0, ilocal=0, jlocal=0;
  int borneInfx = 0, borneInfy = 0, borneMaxx = 0, borneMaxy = 0;
  int valtemp = 0;
  
  int rayon = 1;
  int seuil = 1;
  int seuiltemp = 0;
  
  // Creation de la file d'attente hierarchique
  std::vector<std::list<Node<int, int> > > fah;
  fah.resize(256);
  
  int current_level = 0;
  
  Tin valeurtemp = 0;
  
  Node<int, int> noeudtemp;
  // On empile tous les points du marqueur
  for (i=0; i<w; i++){
    for (j=0; j<h; j++){
      if (imMark(i,j) != 0){
	imFlag(i,j) = seuil;
	valeurtemp = std::min(imMark(i,j), imIn(i,j));
	imOut(i,j) = Tout(valeurtemp);
	Node<int, int> noeud(i, j, valeurtemp, seuil);
	fah[valeurtemp].push_back(noeud);
	if (int(valeurtemp) > current_level) current_level = int(valeurtemp);
      }
    }
  }
  dbi.WriteInfo("Fin de l'initialisation\n");
  
  while(current_level >= 0)
  {
    while (!fah[current_level].empty())
    {
      noeudtemp = fah[current_level].front();
      fah[current_level].pop_front();
      ilocal = noeudtemp.GetXPos();
      jlocal = noeudtemp.GetYPos();
      valtemp = noeudtemp.GetValue();
      seuiltemp = noeudtemp.GetLabel();
      
      if (seuiltemp < iter)
      {
	// On empile les points du voisinage
	borneInfx = std::max(0,ilocal - rayon);
	borneInfy = std::max(0,jlocal - rayon);
	borneMaxx = std::min(w, ilocal + rayon+1);
	borneMaxy = std::min(h, jlocal + rayon+1);
	
	for (ii=borneInfx; ii<borneMaxx; ii++)
	{
	  for (jj=borneInfy; jj<borneMaxy; jj++)
	  {
	    if (imFlag(ii,jj) != seuil)
	    {
	      valeurtemp = std::min(valtemp, int(imIn(ii,jj)));
	      Node<int, int> noeud(ii, jj, valeurtemp, seuiltemp+1);
	      fah[valeurtemp].push_back(noeud);
	      imFlag(ii,jj) = seuil;
	      imOut(ii,jj) = Tout(valeurtemp);
	    }
	  }
	}
      }
    }
    current_level--;
  }
  dbi.WriteOut("Leaving ImUnderBuild_WithIterCount\n");
}

/*! 
 * @brief Computes the morphological gradient (neighborhood maximum minus neighborhood minimum).
 * @param imIn : input image
 * @param method : specifies if Basic or Fast implementation method (Basic supports different shape of structuring element and might be longer, Fast is faster but supports only rectangular shaped structuring elements of varying radius)
 * @param param : vector of parameters (either a list as structuring element for the Basic method, or the height and width for the rectangular structuring element for the Fast method)
 * @param imGradient : image containing the gradient information
 */
template<typename Tin, typename Tout> void ImMorphoGradient(const DipImage<Tin> & imIn, const MorphoMethod method, const std::vector<ParamT> & param, DipImage<Tout> & imOut)
{
	dbi.WriteInfo("Entering ImMorphoGradient\n");
	int Nparam = param.size();
	int w = imIn.width(), h = imIn.height();
	
	switch (method)
	{
		case Basic:
		{
			// LS: NOT IN ORGINAL CODE, I ADDED IT ON MAY 26, 2014
			if ((Nparam==1) && (param.front().type()==typeid(std::list<std::pair<int, int> >)))
			{
				dbi.WriteInfo("Method : Basic\n");
				cout << "You called ImMorphoGradient with method Basic" << endl;
				
				DipImage<Tin> imDilate(w,h); imDilate.init(0);
				ImDilate(imIn, Basic, param, imDilate);
				DipImage<Tin> imErode(w,h); imErode.init(0);
				ImErode(imIn, Basic, param, imErode);
				ImSubImage(imDilate, NoClip, imErode, imOut); 		// LS: WATCH OUT FOR SUBSEQUENT CHANGES IN MATH.HPP
								
			}
		}
		case Fast:
		{
			if(Nparam!=1){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("the number of parameters is not correct");
			}
			if(param[0].type() != typeid(int)){
				//~ throw dip_file_exception() << function_name_info(__FUNCTION__) << comment_info("parameter data type is not correct");
			}
			
			int rayon = boost::get<int> (param.front());
				
			int i = 0, j = 0, ii = 0, jj = 0;
			int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
			int minimum = std::numeric_limits<int>::max();
			int maximum = 0;
			
			for (j=0; j<h; j++)
			{
				borneInfy = std::max(0,j-rayon);
				borneMaxy = std::min(h, j+rayon+1);
				for (i=0; i<w; i++)
				{
					borneInfx = std::max(0,i-rayon);
					borneMaxx = std::min(w, i+rayon+1);
			
					minimum = std::numeric_limits<int>::max();
					maximum = 0;
			
					for (ii=borneInfx; ii<borneMaxx; ii++)
					{
						for (jj=borneInfy; jj<borneMaxy; jj++)
						{
							maximum = std::max(maximum, int(imIn(ii,jj)));
							minimum = std::min(minimum, int(imIn(ii,jj)));
						}
					}
					imOut(i,j) = (Tout)(maximum - minimum);
				}
			}
			break;
		}
			
	}
	dbi.WriteInfo("Leaving ImMorphoGradient. \n");
}

/**
 *
 * @param imIn
 * @param h
 * @param imOut
 */
template<typename Tin, typename Tout> void ImHMaxima(const DipImage<Tin> & imIn, const Tin h, DipImage<Tout> & imOut)
{
	dbi.WriteEnter("Entering ImHMaxima\n");
	
	typedef Tin TypeIn;
	//typedef Tout TypeOut;
	
	DipImage<TypeIn> imSub(imIn.width(),imIn.height());
	ImSubConst(imIn, Clip, h, imSub);
	ImUnderBuild(imSub, imIn, imOut);
	
	dbi.WriteOut("Leaving ImHMaxima\n");
}


/*! 
 * \brief This function performs a segmentation of an image (imIn) using a morphological flooding constrained by a mask image.
 * \param imIn : image to be segmented
 * \param imMask : mask information
 * \param rayon : ray of the neighborhood
 * \param imLabel : seed for flooding
 * \param imOut : enhanced image
 */
template<typename Tin, typename Tmask, typename Tlabel>
void ImSegmentation(const DipImage<Tin> & imIn, const DipImage<Tmask> & imMask, const int rayon, const DipImage<Tlabel> & imLabel, DipImage<Tlabel> & imOut)
{
	dbi.WriteEnter("Entering ImSegmentation\n");
	imOut.init(0);
	int w = imLabel.width(), h = imLabel.height();
	
	typedef Tin TypeIn;
	typedef Tlabel TypeLabel;
	//typedef Tmask TypeMask;
	typedef Tlabel TypeOut;
	
	int maxi = std::numeric_limits<Tin>::max();
	dbi.WriteInfo("maxi = %d\n", int(maxi)+1);
	
	std::vector< std::queue<Node<TypeIn, TypeLabel> > > hq;
	hq.resize(int(maxi)+1);
// 	hq.resize(256);
	
	int borneInfx = 0, borneMaxx = 0, borneInfy = 0, borneMaxy = 0;
	
	// deque filling
	int i = 0, j = 0, ii = 0, jj = 0;
	int current_level = 0; // Cette valeur doit rester un int (sinon, il y a des problemes lorsqu'on atteind la valeur max)
	TypeLabel current_label = 0;
	TypeIn current_value = 0;
	Node<TypeIn, TypeLabel> current_node;
	
	for (i=0; i<w; i++){
		for (j=0; j<h; j++){
			if (imMask(i,j) != 0){
				current_label = imLabel(i,j);
				if (current_label != 0)
				{
					current_value = imIn(i,j);
					current_node.SetPos(i,j);
					current_node.SetValue(current_value);
					current_node.SetLabel(current_label);
					hq[current_value].push(current_node);
					if (int(current_value) < current_level) current_level = int(current_value);
				}
// 				imOut(i,j) = 0;
			}
		}
	}
	
	// processing
	while (current_level < int(int(maxi) + 1)){
// 	while (current_level < 256){
		while (!hq[current_level].empty()){
			current_node = hq[current_level].front();
			hq[current_level].pop();
			i = current_node.GetXPos();
			j = current_node.GetYPos();
			
			if (current_node.GetLabel() == 0){
				std::cout<<" Problem with node label"<< std::endl;
			}
			
			if (imOut(i,j) == 0)
			{
				imOut(i,j) = (TypeOut)current_node.GetLabel();
				
				borneInfx = std::max(0,i-rayon);
				borneInfy = std::max(0,j-rayon);
				borneMaxx = std::min(w, i+rayon+1);
				borneMaxy = std::min(h, j+rayon+1);
				
				for (ii=borneInfx; ii<borneMaxx; ii++){
					for (jj=borneInfy; jj<borneMaxy; jj++){
						if (imMask(ii,jj) != 0){
							if(imOut(ii,jj) == 0){
								current_value = imIn(ii,jj);
								hq[current_value].push(Node<TypeIn,TypeLabel>(ii,jj,current_value,current_node.GetLabel()));
								if (int(current_value) < current_level){
									current_level = int(current_value);
								}
							}
						}
					}
				}
			}
		}
// 		dbi.WriteInfo("level = %d -> %d\n", current_level, int(current_level+1));
		current_level++;
		
	}
// 	return imOut;
	dbi.WriteOut("Leaving ImSegmentation\n");
}

/*! 
 * @brief This function performs a segmentation of an image (imIn) using a morphological flooding.
 * @param imIn : image to be segmented
 * @param rayon : ray of the neighborhood
 * @param imLabel : seed for flooding
 * @param imOut : enhanced image
 */
template<typename T> void ImSegmentation(DipImage<T> & imIn, int rayon, DipImage<unsigned int> & imLabel, DipImage<unsigned int> & imOut)
{
	dbi.WriteEnter("Entering ImSegmentation\n");
	int imSize_x = imLabel.width();
	int imSize_y = imLabel.height();
	DipImage<unsigned char> imMask(imSize_x, imSize_y); imMask.init(255);
	
	ImSegmentation_WithMask(imIn, imMask, rayon, imLabel, imOut);
	
	dbi.WriteOut("Leaving ImSegmentation\n");
}


/*! 
 * @brief Computes the morphological dilation of the input image with respect to imMask
 * @param imIn : input image
 * @param imMask : points to be processed
 * @param rayon : neighborhood radius
 * @param imDilate : dilated image
 */
template<typename T> void ImFastMorphoGradient_WithMask(const DipImage<T> & imIn, const DipImage<unsigned char> & imMask, const int rayon, DipImage<T> & imGradient)
{
  dbi.WriteEnter("Entering ImFastMorphoGradient_WithMask\n");
  typedef T TypeTmp;

  int w = imIn.width(), h = imIn.height();

  int i=0, j=0, ii=0, jj=0;
  int borneInfx = 0, borneMaxx = 0;


  TypeTmp maxligne = 0;
  TypeTmp maxglobal = 0;
  TypeTmp minligne = std::numeric_limits<TypeTmp>::min();
  TypeTmp minglobal = std::numeric_limits<TypeTmp>::min();
  TypeTmp valtemp = 0;
  TypeTmp valtempmin = std::numeric_limits<TypeTmp>::min();

  std::list<TypeTmp> maxlist; if(!maxlist.empty()) {printf("Attention: la liste n'est pas vide\n"); maxlist.clear();}
  std::list<TypeTmp> minlist; if(!minlist.empty()) {printf("Attention: la liste n'est pas vide\n"); minlist.clear();}

  typename std::list<TypeTmp>::iterator itmaxlist = maxlist.begin();
  typename std::list<TypeTmp>::iterator itminlist = minlist.begin();
  typename std::list<TypeTmp>::iterator itmaxlistend = maxlist.end();
  typename std::list<TypeTmp>::iterator itminlistend = minlist.end();

  for (i=0; i<w; i++){
    /// Premier point de la colone
    j=0;
    maxlist.clear();
    minlist.clear();
    maxglobal = 0;
    minglobal = std::numeric_limits<TypeTmp>::max();
    maxligne = 0;
    minligne = std::numeric_limits<TypeTmp>::max();
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
      maxligne = 0;
      minligne = std::numeric_limits<TypeTmp>::max();
      for (ii= borneInfx; ii< borneMaxx; ii++){
	if (imMask(ii,jj) != 0){
	  if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
	  if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
	}
      }
      if (maxligne > maxglobal) maxglobal = maxligne;
      if (minligne < minglobal) minglobal = minligne;
      maxlist.push_back(maxligne);
      minlist.push_back(minligne);
    }
    if (imMask(i,j) != 0){
      imGradient(i,j) = (T)(maxglobal - minglobal);
    }
    /// Fin du traitement du premier point de la colone
    
    /// Iteration sur la colone
    for (j=1; j<h; j++){
      if (j<=rayon){
	// Ajout d'une nouvelle ligne
	jj = j+rayon;
	maxligne = 0;
	minligne = std::numeric_limits<TypeTmp>::max();
	for (ii=borneInfx; ii<borneMaxx; ii++){
	  if (imMask(ii,jj) != 0){
	    if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
	    if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
	  }
	}
	if (maxligne > maxglobal) maxglobal = maxligne;
	if (minligne < minglobal) minglobal = minligne;
	maxlist.push_back(maxligne);
	minlist.push_back(minligne);
      }
      else if ((h-j)<=rayon)
      {
	// On retire une ligne
	valtemp = maxlist.front();
	maxlist.pop_front();
	if (valtemp == maxglobal)
	{
	  maxglobal = 0;
	  itmaxlist = maxlist.begin();
	  itmaxlistend = maxlist.end();
	  for(;itmaxlist != itmaxlistend; itmaxlist++)
	    if (*itmaxlist > maxglobal)
	      maxglobal = *itmaxlist;
	}
	
	valtempmin = minlist.front();
	minlist.pop_front();
	if (valtempmin == minglobal)
	{
	  minglobal = std::numeric_limits<TypeTmp>::max();
	  itminlist = minlist.begin();
	  itminlistend = minlist.end();
	  for(;itminlist != itminlistend; itminlist++)
	    if (*itminlist < minglobal)
	      minglobal = *itminlist;
	}
      }
      else
      {
	// On ajoute et on retire une ligne
	jj = j+rayon;
	maxligne = 0;
	minligne = std::numeric_limits<TypeTmp>::max();
	for (ii=borneInfx; ii<borneMaxx; ii++){
	  if (imMask(ii,jj) != 0){
	    if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
	    if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
	  }
	}
	if (maxligne > maxglobal) maxglobal = maxligne;
	if (minligne < minglobal) minglobal = minligne;
	maxlist.push_back(maxligne);
	minlist.push_back(minligne);

	valtemp = maxlist.front();
	maxlist.pop_front();
	if (valtemp == maxglobal)
	{
	  maxglobal = 0;
	  itmaxlist = maxlist.begin();
	  itmaxlistend = maxlist.end();
	  for(;itmaxlist != itmaxlistend; itmaxlist++)
	    if (*itmaxlist > maxglobal)
	      maxglobal = *itmaxlist;
	}
	valtempmin = minlist.front();
	minlist.pop_front();
	if (valtempmin == minglobal)
	{
	  minglobal = std::numeric_limits<TypeTmp>::max();
	  itminlist = minlist.begin();
	  itminlistend = minlist.end();
	  for(;itminlist != itminlistend; itminlist++)
	    if (*itminlist < minglobal)
	      minglobal = *itminlist;
	}
      }
      if (imMask(i,j) != 0){
	      imGradient(i,j) = (T)(maxglobal - minglobal);
      }
    }
  }
  dbi.WriteOut("Leaving ImFastMorphoGradient_WithMask\n");
}


/*! @fn ImFastMorphoGradient_WithMask_MT(const DipImage<T> & imIn, const DipImage<unsigned char> & imMask, const int rayon, DipImage<T> & imGradient)
 * @brief Computes the morphological dilation of the input image with respect to imMask
 * @param imIn : input image
 * @param imMask : points to be processed
 * @param rayon : neighborhood radius
 * @param imDilate : dilated image
 */
template<typename T> void ImFastMorphoGradient_WithMask_MT(const DipImage<T> & imIn, const DipImage<unsigned char> & imMask, const int rayon, DipImage<T> & imGradient)
{
	dbi.WriteEnter("Entering ImFastMorphoGradient_WithMask_MT\n");
	typedef T TypeTmp;

	int w = imIn.width(), h = imIn.height();

	int i=0, j=0, ii=0, jj=0;
	int borneInfx = 0, borneMaxx = 0;

	#pragma omp parallel private( i, j, ii, jj, borneInfx, borneMaxx ) firstprivate( w, h ) shared( imIn, imMask, imGradient )
	{
		TypeTmp maxligne = 0;
		TypeTmp maxglobal = 0;
		TypeTmp minligne = std::numeric_limits<TypeTmp>::min();
		TypeTmp minglobal = std::numeric_limits<TypeTmp>::min();
		TypeTmp valtemp = 0;
		TypeTmp valtempmin = std::numeric_limits<TypeTmp>::min();

		std::list<TypeTmp> maxlist;
		std::list<TypeTmp> minlist;

		typename std::list<TypeTmp>::iterator itmaxlist = maxlist.begin();
		typename std::list<TypeTmp>::iterator itminlist = minlist.begin();
		typename std::list<TypeTmp>::iterator itmaxlistend = maxlist.end();
		typename std::list<TypeTmp>::iterator itminlistend = minlist.end();

		#pragma omp for schedule( guided ) nowait
		for (i=0; i<w; i++){
			/// Premier point de la colone
			j=0;
			maxlist.clear();
			minlist.clear();
			maxglobal = 0;
			minglobal = std::numeric_limits<TypeTmp>::max();
			maxligne = 0;
			minligne = std::numeric_limits<TypeTmp>::max();
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
				maxligne = 0;
				minligne = std::numeric_limits<TypeTmp>::max();
				for (ii= borneInfx; ii< borneMaxx; ii++){
					if (imMask(ii,jj) != 0){
						if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
						if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
					}
				}
				if (maxligne > maxglobal) maxglobal = maxligne;
				if (minligne < minglobal) minglobal = minligne;
				maxlist.push_back(maxligne);
				minlist.push_back(minligne);
			}
			if (imMask(i,j) != 0){
				imGradient(i,j) = (T)(maxglobal - minglobal);
			}
			/// Fin du traitement du premier point de la colone

			/// Iteration sur la colone
			for (j=1; j<h; j++){
				if (j<=rayon){
					// Ajout d'une nouvelle ligne
					jj = j+rayon;
					maxligne = 0;
					minligne = std::numeric_limits<TypeTmp>::max();
					for (ii=borneInfx; ii<borneMaxx; ii++){
						if (imMask(ii,jj) != 0){
							if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
							if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
						}
					}
					if (maxligne > maxglobal) maxglobal = maxligne;
					if (minligne < minglobal) minglobal = minligne;
					maxlist.push_back(maxligne);
					minlist.push_back(minligne);
				}
				else if ((h-j)<=rayon)
				{
					// On retire une ligne
					valtemp = maxlist.front();
					maxlist.pop_front();
					if (valtemp == maxglobal)
					{
						maxglobal = 0;
						itmaxlist = maxlist.begin();
						itmaxlistend = maxlist.end();
						for(;itmaxlist != itmaxlistend; itmaxlist++)
							if (*itmaxlist > maxglobal)
								maxglobal = *itmaxlist;
					}

					valtempmin = minlist.front();
					minlist.pop_front();
					if (valtempmin == minglobal)
					{
						minglobal = std::numeric_limits<TypeTmp>::max();
						itminlist = minlist.begin();
						itminlistend = minlist.end();
						for(;itminlist != itminlistend; itminlist++)
							if (*itminlist < minglobal)
								minglobal = *itminlist;
					}
				}
				else
				{
					// On ajoute et on retire une ligne
					jj = j+rayon;
					maxligne = 0;
					minligne = std::numeric_limits<TypeTmp>::max();
					for (ii=borneInfx; ii<borneMaxx; ii++){
						if (imMask(ii,jj) != 0){
							if (imIn(ii, jj) > maxligne) maxligne = imIn(ii, jj);
							if (imIn(ii, jj) < minligne) minligne = imIn(ii, jj);
						}
					}
					if (maxligne > maxglobal) maxglobal = maxligne;
					if (minligne < minglobal) minglobal = minligne;
					maxlist.push_back(maxligne);
					minlist.push_back(minligne);

					valtemp = maxlist.front();
					maxlist.pop_front();
					if (valtemp == maxglobal)
					{
						maxglobal = 0;
						itmaxlist = maxlist.begin();
						itmaxlistend = maxlist.end();
						for(;itmaxlist != itmaxlistend; itmaxlist++)
							if (*itmaxlist > maxglobal)
								maxglobal = *itmaxlist;
					}
					valtempmin = minlist.front();
					minlist.pop_front();
					if (valtempmin == minglobal)
					{
						minglobal = std::numeric_limits<TypeTmp>::max();
						itminlist = minlist.begin();
						itminlistend = minlist.end();
						for(;itminlist != itminlistend; itminlist++)
							if (*itminlist < minglobal)
								minglobal = *itminlist;
					}
				}
				if (imMask(i,j) != 0){
					imGradient(i,j) = (T)(maxglobal - minglobal);
				}
			}
		}
	}
	dbi.WriteOut("Leaving ImFastMorphoGradient_WithMask_MT\n");
}

#ifdef _GPU_
/*! @fn ImFastMorphoGradient_WithMask_GPU(const DipImage<T> & imIn, const DipImage<unsigned char> & imMask, const int rayon, DipImage<T> & imGradient)
 * @brief Computes the morphological gradient of the input image with respect to imMask, on the GPU
 * @param imIn : input image
 * @param imMask : points to be processed
 * @param rayon : neighborhood radius
 * @param imGradient : dilated image
 */
template<typename T> void ImFastMorphoGradient_WithMask_GPU(const DipGPUImage<T> & imInBI,
                                                            const DipGPUImage<unsigned char> & imMaskBI,
                                                            const int rayon,
                                                            DipGPUImage<T> & imGradientBI)
{
  dbi.WriteEnter("Entering ImFastMorphoGradient_WithMask_GPU\n");

  const T MaxValueStart = std::numeric_limits<T>::min( );
  const T MinValueStart = std::numeric_limits<T>::max( );

  const int width = imInBI.Width( );
  const int height = imInBI.Height( );

  // Do some padding to remove coslty tests
  const unsigned int GroupSize = std::max( DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeX_H" ),
                                           DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeY_V" ) );
  const int NewHeight = Compute::Utils::GetNextMultipleOf( height, (int)GroupSize ) + rayon * 2 + GroupSize;
  const int NewWidth = Compute::Utils::GetNextMultipleOf( width, (int)GroupSize ) + rayon * 2 + GroupSize;
  const int NewNbPixel = NewHeight * NewWidth;

  const cl::Context Context = Compute::CCompute::GetInstance( )->GetContext( );
  cl::CommandQueue Queue = Compute::CCompute::GetInstance()->GetCommandQueue( );
  cl_int Error;
  cl_mem_flags InOutMemFlags = CL_MEM_READ_WRITE;

  // Watch out for MaxWorkGroupSize
  unsigned int LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeY_H" );
  unsigned int LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeX_H" );
  Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
  // Global size must be a multiple of local size
  unsigned int GlobalSizeY = (height + LocalSizeY - 1) / LocalSizeY;
  GlobalSizeY *= LocalSizeY;
  unsigned int GlobalSizeX = (width + LocalSizeX - 1) / LocalSizeX;
  GlobalSizeX *= LocalSizeX;

  // We need new buffers to get our padding
  cl::Buffer ImNewInBuffer( Context,
                            InOutMemFlags, NewNbPixel * sizeof( T ),
                            0, &Error );
  cl::Buffer ImNewMaskBuffer( Context,
                              InOutMemFlags, NewNbPixel * sizeof( unsigned char ),
                              0, &Error );

  cl::Buffer ImMinBuffer( Context,
                          InOutMemFlags, NewNbPixel * sizeof( T ),
                          0, &Error );
  cl::Buffer ImMaxBuffer( Context,
                          InOutMemFlags, NewNbPixel * sizeof( T ),
                          0, &Error );
  cl::Buffer ImMinBuffer_Transp( Context,
                                 InOutMemFlags, NewNbPixel * sizeof( T ),
                                 0, &Error );
  cl::Buffer ImMaxBuffer_Transp( Context,
                                 InOutMemFlags, NewNbPixel * sizeof( T ),
                                 0, &Error );

  cl::Buffer ImNewOutBuffer( Context,
                             InOutMemFlags, NewNbPixel * sizeof( T ),
                             0, &Error );

  cl::Kernel KH = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastGradient_WithMask_H_Pad" );

  cl::Kernel KV = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastGradient_WithMask_V_Pad_Transpose" );

  cl::Kernel KT = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastGradient_TempTranspose" );

  cl::Kernel KZero = DProgramHolder<T, T, unsigned char>::GetInstance( )->GetKernel( "FastGradient_ZeroOutPaddingBuffers" );

  // Zero out padded buffers
  int Index = 0;
  KZero.setArg( Index++, ImNewInBuffer );
  KZero.setArg( Index++, ImNewMaskBuffer );
  KZero.setArg( Index++, ImMinBuffer );
  KZero.setArg( Index++, ImMaxBuffer );
  KZero.setArg( Index++, ImNewOutBuffer );
  KZero.setArg( Index++, MinValueStart );
  KZero.setArg( Index++, MaxValueStart );
  Queue.enqueueNDRangeKernel( KZero, cl::NullRange, cl::NDRange( NewNbPixel ), cl::NullRange, 0, 0 );

  // Copy the image and mask into their new buffers
  cl::size_t<3> SrcOrigin;
  SrcOrigin.push_back( 0 );
  SrcOrigin.push_back( 0 );
  SrcOrigin.push_back( 0 );
  cl::size_t<3> DstOrigin;
  DstOrigin.push_back( rayon * sizeof( T ) );
  DstOrigin.push_back( rayon );
  DstOrigin.push_back( 0 );
  cl::size_t<3> Region;
  Region.push_back( width * sizeof( T ) );
  Region.push_back( height );
  Region.push_back( 1 );
  size_t SrcRowPitch = width * sizeof( T );
  size_t SrcSlicePitch = height * SrcRowPitch;
  size_t DstRowPitch = NewWidth * sizeof( T );
  size_t DstSlicePitch = NewHeight * DstRowPitch;
  Queue.enqueueCopyBufferRect( imInBI.Buffer( ), ImNewInBuffer, SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );
  DstOrigin[ 0 ] = rayon * sizeof( unsigned char );
  Region[ 0 ] = width * sizeof( unsigned char );
  SrcRowPitch = width * sizeof( unsigned char );
  SrcSlicePitch = height * SrcRowPitch;
  DstRowPitch = NewWidth * sizeof( unsigned char );
  DstSlicePitch = NewHeight * DstRowPitch;
  Queue.enqueueCopyBufferRect( imMaskBI.Buffer( ), ImNewMaskBuffer, SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );

  // Horizontal pass
  Index = 0;
  KH.setArg( Index++, ImNewInBuffer );
  KH.setArg( Index++, ImNewMaskBuffer );
  KH.setArg( Index++, ImMinBuffer );
  KH.setArg( Index++, ImMaxBuffer );
  KH.setArg( Index++, NewHeight );
  KH.setArg( Index++, NewWidth );
  KH.setArg( Index++, rayon );
  KH.setArg( Index++, MinValueStart );
  KH.setArg( Index++, MaxValueStart );
  KH.setArg( Index++, LocalSizeX * 2 * sizeof( T ), 0 );
  KH.setArg( Index++, LocalSizeX * 2 * sizeof( unsigned char ), 0 );

  Queue.enqueueNDRangeKernel( KH, cl::NullRange, cl::NDRange( GlobalSizeX, GlobalSizeY ), cl::NDRange( LocalSizeX, LocalSizeY ), 0, 0 );

  // Transpose
  // Watch out for MaxWorkGroupSize
  LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeY_T" );
  LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeX_T" );
  Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
  // Global size must be a multiple of local size
  GlobalSizeY = (NewHeight + LocalSizeY - 1) / LocalSizeY;
  GlobalSizeY *= LocalSizeY;
  GlobalSizeX = (NewWidth + LocalSizeX - 1) / LocalSizeX;
  GlobalSizeX *= LocalSizeX;

  Index = 0;
  KT.setArg( Index++, ImMinBuffer );
  KT.setArg( Index++, ImMaxBuffer );
  KT.setArg( Index++, ImMinBuffer_Transp );
  KT.setArg( Index++, ImMaxBuffer_Transp );
  KT.setArg( Index++, NewHeight );
  KT.setArg( Index++, NewWidth );
  KT.setArg( Index++, LocalSizeY * LocalSizeX * sizeof( T ), 0 );
  KT.setArg( Index++, LocalSizeY * LocalSizeX * sizeof( T ), 0 );

  Queue.enqueueNDRangeKernel( KT, cl::NullRange, cl::NDRange( GlobalSizeX, GlobalSizeY ), cl::NDRange( LocalSizeX, LocalSizeY ), 0, 0 );

  // Vertical pass
  // Watch out for MaxWorkGroupSize
  LocalSizeY = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeY_V" );
  LocalSizeX = DComputeConfig::GetInstance( )->GetConfigValue<unsigned int>( "Morpho.ImFastGradient_WithMask_GPU.LocalSizeX_V" );
  Compute::CCompute::GetInstance( )->TestThrowLocalSize( Queue, LocalSizeX, LocalSizeY );
  // Global size must be a multiple of local size
  GlobalSizeY = (height + LocalSizeY - 1) / LocalSizeY;
  GlobalSizeY *= LocalSizeY;
  GlobalSizeX = (width + LocalSizeX - 1) / LocalSizeX;
  GlobalSizeX *= LocalSizeX;

  Index = 0;
  KV.setArg( Index++, ImMinBuffer_Transp );
  KV.setArg( Index++, ImMaxBuffer_Transp );
  KV.setArg( Index++, ImNewMaskBuffer );
  KV.setArg( Index++, ImNewOutBuffer );
  KV.setArg( Index++, NewHeight );
  KV.setArg( Index++, NewWidth );
  KV.setArg( Index++, rayon );
  KV.setArg( Index++, LocalSizeY * 2 * sizeof( T ), 0 );
  KV.setArg( Index++, LocalSizeY * 2 * sizeof( T ), 0 );

  Queue.enqueueNDRangeKernel( KV, cl::NullRange, cl::NDRange( GlobalSizeY, GlobalSizeX ), cl::NDRange( LocalSizeY, LocalSizeX ), 0, 0 );

  // Copy the relevant sub rect into the output
  SrcOrigin[ 0 ] = rayon * sizeof( T );
  SrcOrigin[ 1 ] = rayon;
  SrcOrigin[ 2 ] = 0;
  DstOrigin[ 0 ] = 0;
  DstOrigin[ 1 ] = 0;
  DstOrigin[ 2 ] = 0;
  Region[ 0 ] = width * sizeof( T );
  SrcRowPitch = NewWidth * sizeof( T );
  SrcSlicePitch = NewHeight * SrcRowPitch;
  DstRowPitch = width * sizeof( T );
  DstSlicePitch = height * DstRowPitch;
  Queue.enqueueCopyBufferRect( ImNewOutBuffer, imGradientBI.Buffer( ), SrcOrigin, DstOrigin, Region, SrcRowPitch, SrcSlicePitch, DstRowPitch, DstSlicePitch, 0, 0 );

  dbi.WriteOut("Leaving ImFastMorphoGradient_WithMask_GPU\n");
}

/*! @fn ImFastMorphoGradient_WithMask_GPU(const DipImage<T> & imIn, const DipImage<unsigned char> & imMask, const int rayon, DipImage<T> & imGradient)
 * @brief Computes the morphological gradient of the input image with respect to imMask, on the GPU
 * @param imIn : input image
 * @param imMask : points to be processed
 * @param rayon : neighborhood radius
 * @param imGradient : dilated image
 */
template<typename T> void ImFastMorphoGradient_WithMask_GPU(const DipImage<T> & imIn,
                                                            const DipImage<unsigned char> & imMask,
                                                            const int rayon,
                                                            DipImage<T> & imGradient)
{
    dbi.WriteEnter("Entering ImFastMorphoGradient_WithMask_GPU\n");

    cl_mem_flags InFlags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
    cl_mem_flags OutFlags = CL_MEM_WRITE_ONLY;

    DipGPUImage<T> imInBI( imIn, InFlags );
    DipGPUImage<unsigned char> imMaskBI( imMask, InFlags );
    DipGPUImage<T> imGradientBI( imGradient, OutFlags );

    // Call the version with buffers
    ImFastMorphoGradient_WithMask_GPU( imInBI, imMaskBI, rayon, imGradientBI );

    // Read the result
    imGradientBI.ToImage( imGradient );

    dbi.WriteOut("Leaving ImFastMorphoGradient_WithMask_GPU\n");
}

#endif // _GPU_

#endif // MORPHO_HPP
