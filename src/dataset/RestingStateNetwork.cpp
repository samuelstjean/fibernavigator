/////////////////////////////////////////////////////////////////////////////
// Name:            RestingState.cpp
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////

#include "RestingStateNetwork.h"

#include "DatasetManager.h"
#include "AnatomyHelper.h"
#include "RTFMRIHelper.h"
#include "../Logger.h"
#include "../gfx/ShaderHelper.h"
#include "../gfx/TheScene.h"
#include "../gui/MyListCtrl.h"
#include "../gui/SceneManager.h"
#include "../misc/nifti/nifti1_io.h"

#include <GL/glew.h>
#include <wx/math.h>
#include <wx/xml/xml.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <vector>
#include <sstream>

#include "../gfx/Image.h"
#include "../gfx/BitmapHandling.h"
#include "../gfx/TextureHandling.h"
#include "../main.h"

//////////////////////////////////////////
//Constructor
//////////////////////////////////////////
RestingStateNetwork::RestingStateNetwork():
m_zMin( 999.0f ),
m_zMax(-999.0f ),
m_alpha( 0.5f),
m_pointSize( 10.0f ),
m_isRealTimeOn( false ),
m_dataType( 16 ),
m_bands( 108 ),
m_corrThreshold( 3.0f ),
m_colorSliderValue( 5.0f ),
m_normalize( true )
{
	m_rowsL = DatasetManager::getInstance()->getRows();
	m_columnsL = DatasetManager::getInstance()->getColumns();
	m_framesL =  DatasetManager::getInstance()->getFrames();

	m_xL = DatasetManager::getInstance()->getVoxelX();
	m_yL = DatasetManager::getInstance()->getVoxelY();
	m_zL =  DatasetManager::getInstance()->getVoxelZ();

	m_datasetSizeL = m_rowsL * m_columnsL * m_framesL;
}

//////////////////////////////////////////
//Destructor
//////////////////////////////////////////
RestingStateNetwork::~RestingStateNetwork()
{
    Logger::getInstance()->print( wxT( "RestingStateNetwork destructor called but nothing to do." ), LOGLEVEL_DEBUG );
}

//////////////////////////////////////////
//Load
//////////////////////////////////////////
bool RestingStateNetwork::load( nifti_image *pHeader, nifti_image *pBody )
{
    m_datasetSize = pHeader->dim[1] * pHeader->dim[2] * pHeader->dim[3];
	m_rows = pHeader->dim[1];
	m_columns = pHeader->dim[2];
	m_frames = pHeader->dim[3];
	m_bands = pHeader->dim[4];

    m_voxelSizeX = pHeader->dx;
    m_voxelSizeY = pHeader->dy;
    m_voxelSizeZ = pHeader->dz;
    
	std::vector<short int> fileFloatData( m_datasetSize * m_bands, 0);

	if(pHeader->datatype == 4)
	{
		short int* pData = (short int*)pBody->data;
		//Prepare the data into a 1D vector, side by side
		for( int i( 0 ); i < m_datasetSize; ++i )
		{
			for( int j( 0 ); j < m_bands; ++j )
			{
				//if(!isnan(pData[j * datasetSize + i]))
					fileFloatData[i * m_bands + j] = pData[j * m_datasetSize + i];
			}
		}
	}
	else
	{
		float* pData = (float*)pBody->data;
		//Prepare the data into a 1D vector, side by side
		for( int i( 0 ); i < m_datasetSize; ++i )
		{
			for( int j( 0 ); j < m_bands; ++j )
			{
				//if(!isnan(pData[j * datasetSize + i]))
					fileFloatData[i * m_bands + j] = pData[j * m_datasetSize + i];
			}
		}
		
	}

	//Assign structure to a 2D vector of timelaps
    createStructure( fileFloatData );

	//Load fMRI sprite texture.
	Image<ColorRGB> TmpImage;
	wxString name = wxT ("fMRI.bmp");

	wxString iconPath = MyApp::iconsPath;
	wxString fullname = iconPath.append(name);
	std::string stlstring = std::string(fullname.mb_str());
	
    //Load the color scheme #1 image and send it to the GPU as a texture.
    LoadBmp(stlstring,TmpImage);
    m_lookupTex = LoadTexture(TmpImage);

	//Logger::getInstance()->print( wxT( "Resting-state network initialized" ), LOGLEVEL_MESSAGE );
    return true;
}

void RestingStateNetwork::loadTxt( wxArrayString fileName )
{
    wxString f = fileName[0];
    ifstream myFile(f);

    if (myFile.is_open())
    {
        string line;
        getline(myFile, line);
        istringstream linestream(line);
        float ID;

        while (linestream >> ID)
        {   
            m_IDs.push_back(ID);
            vector<float> V; 
            m_timeCourseMAP[ID] = V;
        }

        while(getline(myFile, line))
        {
            istringstream linestream(line);
            float value;
            for(size_t i=0; i < m_IDs.size(); i++)
            {
                linestream >> value;
                m_timeCourseMAP[m_IDs[i]].push_back(value);
            }
        }
    }

    m_meansAndSigmas.resize(m_IDs.size());
	//Transpose signal for easy acces of timelaps
    for( size_t s(0); s < m_IDs.size(); ++s )
    {
		calculateMeanAndSigma(m_timeCourseMAP[m_IDs[s]], m_meansAndSigmaMAP[m_IDs[s]]);
    }

    myFile.close();
}

void RestingStateNetwork::setClusters( Anatomy* info )
{ 
	m_pClusterAnatomy = info;
	m_columns = m_pClusterAnatomy->getColumns();
	m_rows = m_pClusterAnatomy->getRows();
	m_frames = m_pClusterAnatomy->getFrames();

	
	
	m_pClusters = m_pClusterAnatomy->getFloatDataset(); 

	for(size_t c=0; c< m_pClusters->size(); c++)
	{
		m_pClusters->at(c) *=  m_pClusterAnatomy->getOldMax();
	}


	for(int x = 0; x < m_rows; x++)
	{
		for(int y = 0; y < m_columns; y++)
		{
			for(int z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				int IDcluster = m_pClusters->at(i);
				if(IDcluster != 0)
					m_voxels.insert(pair<int, Vector>(IDcluster,Vector(x,y,z)));
			}
		}
	}
}


//////////////////////////////////////////
//Create structure
//////////////////////////////////////////
bool RestingStateNetwork::createStructure( std::vector< short int > &i_fileFloatData )
{
	int size = m_rows * m_columns * m_frames;
    m_signal.resize( size );
	m_signalNormalized.resize ( size );
    vector< short int >::iterator it;
    int i = 0;

    //Fetching the directions
    for( it = i_fileFloatData.begin(), i = 0; it != i_fileFloatData.end(); it += m_bands, ++i )
    { 
		m_signal[i].insert( m_signal[i].end(), it, it + m_bands );
    }
	
	//Find min/max for normalization
	vector<float> dataMax, dataMin;
	dataMax.assign(size, -std::numeric_limits<float>::infinity());
	dataMin.assign(size, std::numeric_limits<float>::infinity());
    for( int s(0); s < size; ++s )
    {
		for( int b(0); b < m_bands; ++b )
		{
			if (m_signal[s][b] > dataMax[s])
			{
				dataMax[s] = m_signal[s][b];
			}
			if (m_signal[s][b] < dataMin[s])
			{
				dataMin[s] = m_signal[s][b];
			}
		}
    }

	//Min max Rescale
    for( int s(0); s < size; ++s )
    {
		for( int b(0); b < m_bands; ++b )
		{
			if((m_signal[s][b] == 0 && dataMin[s] == 0) || (m_signal[s][b] == 16767 && dataMin[s] == 16767)) //Ensure that we dont divide by 0.
				m_signalNormalized[s].push_back(0);
			else
				m_signalNormalized[s].push_back ((m_signal[s][b] - dataMin[s]) / (dataMax[s] - dataMin[s]));
		}
    }

	
	m_volumes.resize(m_bands);
	m_meansAndSigmas.resize(size);
	//Transpose signal for easy acces of timelaps
    for( int s(0); s < size; ++s )
    {
		calculateMeanAndSigma(m_signalNormalized[s], m_meansAndSigmas[s]);
    }

	//Create texture made of 1st timelaps
	data.assign(size, 0.0f);

    return true;
}

//////////////////////////////////////////////////////////////////////////
//Get 3D indexes to fill the 1x1x1 texture from a 3x3x3 x,y,z point 
//////////////////////////////////////////////////////////////////////////
vector<int> RestingStateNetwork::get3DIndexes(int x, int y, int z)
{
	std::vector<int> indexes;

	for( int padx = 0; padx < 4; padx++)
	{
		for( int pady = 0; pady < 4; pady++)
		{
			for( int padz = 0; padz < 4; padz++)
			{
				int i = (z * floor(float(m_framesL/m_frames)) + padz) * m_columnsL * m_rowsL + (y *floor(float(m_rowsL/m_rows)) + pady) *m_columnsL + (x * floor(float(m_columnsL/m_columns)) + padx);
				indexes.push_back( i );
			}
		}
	}

	return indexes;
}

//////////////////////////////////////////
//Set raw data texture from sliderValue
//////////////////////////////////////////
void RestingStateNetwork::SetTextureFromSlider(int sliderValue)
{
	std::vector<float> vol(m_datasetSizeL* 3, 0.0f);
	std::vector<int> indexes;
	for( float x = 0; x < m_columns; x++)
	{
		for( float y = 0; y < m_rows; y++)
		{
			for( float z = 0; z < m_frames; z++)
			{
				int i = z * m_columns * m_rows + y *m_columns + x;
				int s = z * floor(float(m_framesL/m_frames)) * m_columnsL * m_rowsL + y *floor(float(m_rowsL/m_rows)) *m_columnsL + x * floor(float(m_columnsL/m_columns)); // O
				
				vol[s*3] = m_signalNormalized[i][sliderValue];
				vol[s*3 + 1] = m_signalNormalized[i][sliderValue];
				vol[s*3 + 2] = m_signalNormalized[i][sliderValue];

				//Patch arround for 1x1x1
				if(m_framesL != m_frames)
				{
					indexes = get3DIndexes(x,y,z);
					for(unsigned int s = 0; s < indexes.size(); s++)
					{
						int id = indexes[s];
						vol[id*3] = m_signalNormalized[i][sliderValue];
						vol[id*3 + 1] = m_signalNormalized[i][sliderValue];
						vol[id*3 + 2] = m_signalNormalized[i][sliderValue];
					}
				}
			}
		}
	}

	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(vol);
	pNewAnatomy->generateTexture();
}

//////////////////////////////////////////////////////////////////////////////////////////
//Set texture from Network fmri: NOTE: doesnt work functionally yet, data should be set
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::SetTextureFromNetwork()
{
	Anatomy* pNewAnatomy = (Anatomy *)DatasetManager::getInstance()->getDataset( m_index );
	pNewAnatomy->setFloatDataset(data);
	pNewAnatomy->generateTexture();
}

//////////////////////////////////////////////////////////////////////////////////////////
//Initiate the seed-based algorithm
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::seedBased()
{
	m_3Dpoints.clear();
	m_zScores.clear();
	m_zMin = 999.0f;
	m_zMax = -999.0f;

    SelectionTree::SelectionObjectVector selObjs = SceneManager::getInstance()->getSelectionTree().getAllObjects();

	for( unsigned int b = 0; b < selObjs.size(); b++ )
	{
		int x = (int)(floor(selObjs[b]->getCenter().x));
		int y = (int)(floor(selObjs[b]->getCenter().y));
		int z = (int)(floor(selObjs[b]->getCenter().z));
		
		int i = z * m_columns * m_rows + y *m_columns + x;
		int ID = m_pClusters->at(i);
		if(ID !=0)
			correlate(ID);
	}

	render3D(true);
	RTFMRIHelper::getInstance()->setRTFMRIDirty(false);
}

//////////////////////////////////////////////////////////////////////////////////////////
//Rendering function, for both 3D sprites and textures options.
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::render3D(bool recalculateTexture)
{
	if( m_zScores.size() > 0 )
    {
		//Apply ColorMap
		for(std::map<int,float>::iterator it = m_zScores.begin(); it != m_zScores.end(); it++)
		{
			float R,G,B;

			float mid = (m_zMin + m_zMax) / 2.0f;
			float v = (it->second - m_zMin) / (m_zMax - m_zMin);
			if(it->second < mid)
			{
				R = (it->second - m_zMin) / (mid - m_zMin);
				G = 1.0f;
				B = 0.0f;
			}
			else
			{
				R = 1.0f;
				G = 1 - (v);
				B = 0.0f;
			}

			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_POINT_SPRITE);
			glPointSize(it->second / m_zMax * m_pointSize + 1.0f);
			glColor4f(R,G,B,it->second / m_zMax *m_alpha);

			//Illuminate voxels for each id.

			
			std::pair <std::multimap<int, Vector>::iterator, std::multimap<int,Vector>::iterator> ret;
			ret = m_voxels.equal_range(it->first);

			for(std::multimap<int, Vector>::iterator itClust = ret.first; itClust != ret.second; itClust++)
			{
				glBegin(GL_POINTS);
					glVertex3f(itClust->second.x, itClust->second.y, itClust->second.z);
				glEnd();
			}

			//glDisable( GL_TEXTURE_2D );
			glDisable(GL_POINT_SPRITE);
			glDisable(GL_BLEND);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//Correlation function given a position, with all other time series
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::correlate(int ID)
{

	//Signal DEPART
	//std::vector<float> meanSignal = m_timeCourseMAP[ID];

	//Mean et sigma de DEPART
	//std::pair<float, float> RefMeanAndSigma = m_meansAndSigmaMAP[ID];

	std::vector<pair<float,int> >corrFactors;
	float corrSum = 0.0f;
	int nb = 0;
	int size;

	//Correlate with rest of the brain, i.e find corr factors
	for(map<int,vector<float> >::iterator it = m_timeCourseMAP.begin() ; it != m_timeCourseMAP.end(); it++)
	{
		float num = 0.0f;
		float denum = 0.0f;
		size = it->second.size();
		//int corri = 0;
					
		for(int t = 0; t < size; t++)
		{
			num += (m_timeCourseMAP[ID][t] - m_meansAndSigmaMAP[ID].first) * ( it->second[t] - m_meansAndSigmaMAP[it->first].first);
		}
		float value = num / ( m_meansAndSigmaMAP[ID].second * m_meansAndSigmaMAP[it->first].second);
		value /= (size);
		corrSum+=value;
				
		//if(value > 0)
		//{
		//	
		//	corrFactors[i] = value;
		//	nb++;
		//}
		//else
		//	corrFactors[i] = -1;
		corrFactors.push_back(make_pair<float, int>(value, it->first));
	}


	//Find mean and sigma of all corr factors.
	float meanCorr = corrSum / size;
	float sigma = 0.0f;
	
	for(int c=0; c<corrFactors.size(); c++)
	{
		sigma += (corrFactors[c].first - meanCorr)*(corrFactors[c].first - meanCorr);	
	}

	//Calculate z-scores, and save them.
	sigma /= size;
	sigma = sqrt(sigma);

	for(int cc=0; cc<corrFactors.size(); cc++)
	{
		float zScore = (corrFactors[cc].first - meanCorr) / sigma;
		if(zScore < m_zMin && zScore > 0.0f)
			m_zMin = zScore;
		if(zScore > m_zMax)
			m_zMax = zScore;
		if(zScore > m_corrThreshold)
		{
			//m_3Dpoints.push_back(std::pair<Vector,float>(Vector(x,y,z),zScore));
			m_zScores[corrFactors[cc].second] = zScore;
		}
	}			
}

//////////////////////////////////////////////////////////////////////////////////////////
//Calculate Mean and Sigma for the signal inside the box
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::calculateMeanAndSigma(std::vector<float> signal, std::pair<float, float>& params)
{
	float mean = 0.0f;
	float sigma = 0.0f;
	
	//mean
	for(unsigned int i=0; i < signal.size(); i++)
	{
		mean+=signal[i];
	}
	mean /= signal.size();

	//sigma
    for(unsigned int i = 0; i < signal.size(); i++)
    {
         sigma += (signal[i] - mean) * (signal[i] - mean) ;
    }
    sigma /= signal.size();

	params.first = mean;
	params.second = sqrt(sigma);
}