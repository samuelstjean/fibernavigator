/////////////////////////////////////////////////////////////////////////////
// Name:            RestingStateNetwork.h
// Author:          Maxime Chamberland
/////////////////////////////////////////////////////////////////////////////
#ifndef RESTINGSTATENETWORK_H_
#define RESTINGSTATENETWORK_H_

#include "Anatomy.h"
#include "DatasetInfo.h"
#include "DatasetIndex.h"
#include "../misc/nifti/nifti1_io.h"
#include <map>

class RestingStateNetwork
{
public:

    // Constructor/Destructor
    RestingStateNetwork();
    virtual ~RestingStateNetwork();
	void setNetworkInfo( DatasetIndex index ) { m_index = index; }
	void SetCorrThreshold( float thresh ) { m_corrThreshold = thresh; }
	void SetSizePSliderValue (float value ) { m_pointSize = value; }
	void SetAlphaSliderValue (float value ) { m_alpha = value; }
	void SetExport (bool value) { m_export = value; }
	void render3D(bool recalculateTexture);
	void seedBased();                           
	void clear3DPoints()                           { m_Ztext.clear(); m_Ztext.resize( m_pClusters->size() ); m_zScores.clear(); }
    void setClusters( Anatomy* info );            
    void loadTxt( wxArrayString fileName );

	std::vector< float >* getZscores() { return &m_Ztext; }

private:

	void correlate(int ID);
	void calculateMeanAndSigma(std::vector<float> signal, std::pair<float, float>& params);
	std::vector<std::pair< float, float > > m_meansAndSigmas; 

	float m_zMin;
	float m_zMax;
	float m_alpha;
	float m_pointSize;

    int m_dataType;
	int m_rows;
	int m_columns;
	int m_frames;
	int m_bands;
	int m_datasetSize;
	float m_voxelSizeX;
    float m_voxelSizeY;
    float m_voxelSizeZ;
	DatasetIndex m_index;
	float m_corrThreshold;
	bool m_export; 

    Anatomy     *m_pClusterAnatomy; //Clusters anatomy
    std::vector<float> *m_pClusters; //Vector where each position is associated to a cluster
    std::map<int, std::vector<float> > m_timeCourseMAP; //Mapping a cluster ID to a timecourse vector
    std::vector<float> m_IDs; //List of clusters ID
    std::map<int, std::pair< float, float > > m_meansAndSigmaMAP; //Mapping a cluster ID to a pair mean/sigma of its timecourse.
	std::multimap<int, Vector> m_voxels; //multimap, mapping a cluster ID to all the pixels it contains (voxels to illuminate following the correlation step)
	std::map<int,float> m_zScores; //Map of cluster ID to zScores, telling us which cluster to illuminate.
	std::vector< float > m_Ztext; //Texture for export

};

#endif /* RESTINGSTATENETWORK_H_ */
