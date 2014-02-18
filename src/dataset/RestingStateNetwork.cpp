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
m_zMax( 0.0f ),
m_alpha( 0.5f),
m_pointSize( 10.0f ),
m_dataType( 16 ),
m_bands( 108 ),
m_corrThreshold( 3.0f ),
m_export( false )
{

}

//////////////////////////////////////////
//Destructor
//////////////////////////////////////////
RestingStateNetwork::~RestingStateNetwork()
{
    Logger::getInstance()->print( wxT( "RestingStateNetwork destructor called but nothing to do." ), LOGLEVEL_DEBUG );
}

//////////////////////////////////////////
//Load Timecourses
//////////////////////////////////////////
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
	m_Ztext.resize( m_pClusters->size() );

	for(size_t c=0; c< m_pClusters->size(); c++)
	{
		m_pClusters->at(c) *=  m_pClusterAnatomy->getOldMax();
	}

	for(int x = 0; x < m_columns; x++)
	{
		for(int y = 0; y < m_rows; y++)
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

//////////////////////////////////////////////////////////////////////////////////////////
//Initiate the seed-based algorithm
//////////////////////////////////////////////////////////////////////////////////////////
void RestingStateNetwork::seedBased()
{
	m_zScores.clear();
	m_zMin = 999.0f;
	m_zMax = 0.0f;

    SelectionTree::SelectionObjectVector selObjs = SceneManager::getInstance()->getSelectionTree().getAllObjects();

	for( unsigned int b = 0; b < selObjs.size(); b++ )
	{
		int x = floor(selObjs[b]->getCenter().x);
		int y = floor(selObjs[b]->getCenter().y);
		int z = floor(selObjs[b]->getCenter().z);
		
		int i = z * m_columns * m_rows + y *m_columns + x;
		int ID = m_pClusters->at(i);
		if(ID != 0)
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
			bool render = true;

			float quart = 1.0f*(m_zMin + m_zMax) / 4.0f;
			float mid = (m_zMin + m_zMax) / 2.0f;
			float trois_quart = 3.0f* (m_zMin + m_zMax) / 4.0f;
			float v = (it->second - m_zMin) / (m_zMax - m_zMin);

			if( it->second < quart )
			{
                R = (it->second - m_zMin) / (quart - m_zMin);
                G = 0.0f;
                B = 0.0f;
                render = false;
			}
			else if(it->second >= quart && it->second < trois_quart)
			{
				R = 1.0f;
				G = (it->second - quart) / (trois_quart - quart);
				B = 0.0f;
			}
			else
			{
				R = 1.0f;
				G = 1.0f;
				B = v;
			}

			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_POINT_SPRITE);
			glPointSize(it->second / (m_zMax *1.5f) * m_pointSize + 1.0f);
			glColor4f(R,G,B,it->second / (1.5f*m_zMax) *m_alpha);

			//Illuminate voxels for each id.
			if(render)
			{
				std::pair <std::multimap<int, Vector>::iterator, std::multimap<int,Vector>::iterator> ret;
				ret = m_voxels.equal_range(it->first);

				for(std::multimap<int, Vector>::iterator itClust = ret.first; itClust != ret.second; itClust++)
				{
					int pX = itClust->second.x;
					int pY = itClust->second.y;
					int pZ = itClust->second.z;

					glBegin(GL_POINTS);
						glVertex3f(pX,pY,pZ);
					glEnd();

					if(m_export)
					{
						int i = pZ * m_columns * m_rows + pY *m_columns + pX;
						m_Ztext[i] = it->second;
					}
				}
				render = true;
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
		value /= size;

		if(value > 0)
		{
			corrSum+=value;
			corrFactors.push_back(make_pair<float, int>(value, it->first));
			nb++;
		}

	}

	//Find mean and sigma of all corr factors.
	float meanCorr = corrSum / nb;
	float sigma = 0.0f;
	
	for(size_t c=0; c<corrFactors.size(); c++)
	{
		sigma += (corrFactors[c].first - meanCorr)*(corrFactors[c].first - meanCorr);	
	}

	//Calculate z-scores, and save them.
	sigma /= nb;
	sigma = sqrt(sigma);

	for(size_t cc=0; cc<corrFactors.size(); cc++)
	{
		if(corrFactors[cc].first > 0)
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