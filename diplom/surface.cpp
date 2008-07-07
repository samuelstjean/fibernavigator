#include "surface.h"
#include "math.h"
#include "Fantom/FMatrix.hh"
#include "Fantom/FBSplineSurface.hh"
#include "surface.h"
#include "point.h"
#include "myListCtrl.h"

#include "GL/glew.h"

Surface::Surface(wxTreeCtrl* treeWidget, wxTreeItemId tPointId)
{
	m_treeWidget = treeWidget;
	m_tPointId = tPointId;

	m_radius = 3.0;
	m_my = 2.0;
	m_numDeBoorRows = 8;
	m_numDeBoorCols = 8;
	m_order = 4;
	m_sampleRateT = m_sampleRateU = 0.2;

	m_show = true;
	m_showFS = true;
	m_useTex = true;

	m_type = Surface_;
	m_threshold = 0.5;
	m_name = wxT("spline surface");
}

FTensor Surface::getCovarianceMatrix(std::vector< std::vector< double > > points)
{
	FTensor result(3,2,true);
	m_xAverage = m_yAverage = m_zAverage = 0;

	std::vector< std::vector< double > >::iterator pointsIt;
	for( pointsIt = points.begin(); pointsIt != points.end(); pointsIt++)
	{
		std::vector< double > dmy = *pointsIt;
		m_xAverage += dmy[0];
		m_yAverage += dmy[1];
		m_zAverage += dmy[2];
	}

	m_xAverage /= points.size();
	m_yAverage /= points.size();
	m_zAverage /= points.size();

	/*
    	 /        \
   		| XX XY XZ |
   		| YX YY YZ |
   		| ZX ZY ZZ |
    	 \        /
	 */
	for( pointsIt = points.begin(); pointsIt != points.end(); pointsIt++)
	{
		std::vector< double > dmy = *pointsIt;

		result(0,0) += (dmy[0] - m_xAverage) * (dmy[0] - m_xAverage); //XX
		result(0,1) += (dmy[0] - m_xAverage) * (dmy[1] - m_yAverage); //XY
		result(0,2) += (dmy[0] - m_xAverage) * (dmy[2] - m_zAverage); //XZ

		result(1,1) += (dmy[1] - m_yAverage) * (dmy[1] - m_yAverage); //YY
		result(1,2) += (dmy[1] - m_yAverage) * (dmy[2] - m_zAverage); //YZ

		result(2,2) += (dmy[2] - m_zAverage) * (dmy[2] - m_zAverage); //ZZ
	}

	result(1,0) = result(0,1);
	result(2,0) = result(0,2);
	result(2,1) = result(1,2);

	return result;
}

void Surface::getSplineSurfaceDeBoorPoints(std::vector< std::vector< double > > &givenPoints, std::vector< std::vector< double > > &deBoorPoints, int numRows, int numCols)
{
	double xMin = givenPoints[0][0];
	double xMax = givenPoints[0][0];
	double zMin = givenPoints[0][2];
	double zMax = givenPoints[0][2];

	std::vector< std::vector< double > >::iterator givenPointsIt;

	for( givenPointsIt = givenPoints.begin(); givenPointsIt != givenPoints.end(); givenPointsIt++)
	{
		std::vector< double > dmy = *givenPointsIt;
		if(dmy[0] < xMin)
			xMin = dmy[0];
		if(dmy[0] > xMax)
			xMax = dmy[0];
		if(dmy[2] < zMin)
			zMin = dmy[2];
		if(dmy[2] > zMax)
			zMax = dmy[2];
	}

	double dX = (xMax - xMin) / (numCols - 1);
	double dZ = (zMax - zMin) / (numRows - 1);


	for( int row = 0; row < numRows; row++)
		for( int col = 0; col < numCols; col++)
		{
			std::vector< double > dmy;
			double x = xMin + dX * col;
			double z = zMin + dZ * row;
			dmy.push_back(x);

			double y = 0;
			double numerator = 0, denominator = 0;

			//<local shepard with franke-little-weights>
			for( givenPointsIt = givenPoints.begin(); givenPointsIt != givenPoints.end(); givenPointsIt++)
			{
				std::vector< double > dmy1 = *givenPointsIt;
				FArray dmyArray(dmy1);
				dmyArray[1] = 0;
				FArray thisPoint(x,0,z);

				double xi; //greek alphabet

				if( thisPoint.distance(dmyArray) < m_radius)
					xi = 1 - thisPoint.distance(dmyArray) / m_radius;
				else
					xi = 0;

				numerator += (pow(xi, m_my) * dmy1[1]);
				denominator += (pow(xi, m_my));
			}
			if( denominator == 0)
				y = 0;
			else
				y = numerator / denominator;
			//</local shepard with franke-little-weights>
			dmy.push_back(y);
			dmy.push_back(z);

			deBoorPoints.push_back(dmy);
		}
  return;
}

void Surface::execute (std::vector< std::vector< double > > givenPoints)
{

	std::vector< std::vector< double > > deBoorPoints;
	std::vector< std::vector< double > > splinePoints;

	FTensor myTensor = getCovarianceMatrix(givenPoints);

	FArray eigenValues(3);
	FArray eigenVectors[3];
	eigenVectors[0] = FArray(3);
	eigenVectors[1] = FArray(3);
	eigenVectors[2] = FArray(3);

	myTensor.getEigenSystem(eigenValues, eigenVectors);

	eigenVectors[0].normalize();
	eigenVectors[1].normalize();
	eigenVectors[2].normalize();

	FTensor::sortEigenvectors(eigenValues, eigenVectors);

	FMatrix transMatrix = FMatrix(3,3);
	transMatrix(0,0) = eigenVectors[1][0];
	transMatrix(0,1) = eigenVectors[1][1];
	transMatrix(0,2) = eigenVectors[1][2];

	transMatrix(1,0) = eigenVectors[2][0];
	transMatrix(1,1) = eigenVectors[2][1];
	transMatrix(1,2) = eigenVectors[2][2];

	transMatrix(2,0) = eigenVectors[0][0];
	transMatrix(2,1) = eigenVectors[0][1];
	transMatrix(2,2) = eigenVectors[0][2];

	std::vector< std::vector< double > >::iterator pointsIt;

	//translate and orientate given points to origin
	for( pointsIt = givenPoints.begin(); pointsIt != givenPoints.end(); pointsIt++)
	{
		(*pointsIt)[0] -= m_xAverage;
		(*pointsIt)[1] -= m_yAverage;
		(*pointsIt)[2] -= m_zAverage;

		FArray dmy(*pointsIt);

		FVector result = transMatrix * dmy;
		(*pointsIt)[0] = result[0];
		(*pointsIt)[1] = result[1];
		(*pointsIt)[2] = result[2];
	}

	//get de Boor points using shepard's method
	getSplineSurfaceDeBoorPoints(givenPoints, deBoorPoints, m_numDeBoorRows, m_numDeBoorCols);

	//translate and orientate de Boor points back
	transMatrix.invert();
	for( pointsIt = deBoorPoints.begin(); pointsIt != deBoorPoints.end(); pointsIt++)
	{
		FArray dmy(*pointsIt);

		FVector result = transMatrix * dmy;
		(*pointsIt)[0] = result[0];
		(*pointsIt)[1] = result[1];
		(*pointsIt)[2] = result[2];

		(*pointsIt)[0] += m_xAverage;
		(*pointsIt)[1] += m_yAverage;
		(*pointsIt)[2] += m_zAverage;
	}

	FBSplineSurface splineSurface(m_order, m_order, deBoorPoints, m_numDeBoorCols, m_numDeBoorRows);

	splineSurface.samplePoints(splinePoints, m_sampleRateT, m_sampleRateU);

	std::vector< double > positions;
	for( std::vector< std::vector< double > >::iterator posIt = splinePoints.begin(); posIt != splinePoints.end(); posIt++)
	{
		positions.push_back((*posIt)[0]);
		positions.push_back((*posIt)[1]);
		positions.push_back((*posIt)[2]);
	}

	//shared_ptr< FPositionSet > positionSet( new FPositionSet3DArbitrary( positions ));
	std::vector< int > vertices;

	int renderpointsPerCol = splineSurface.getNumSamplePointsU();
	int renderpointsPerRow = splineSurface.getNumSamplePointsT();

	for(int z = 0; z < renderpointsPerCol - 1; z++)
	{
		for(int x = 0; x < renderpointsPerRow - 1; x++)
		{
			vertices.push_back(z * renderpointsPerCol + x);
			vertices.push_back(z * renderpointsPerCol + x + 1);
			vertices.push_back((z+1) * renderpointsPerCol + x);

			vertices.push_back((z+1) * renderpointsPerCol + x);
			vertices.push_back(z * renderpointsPerCol + x + 1);
			vertices.push_back((z+1) * renderpointsPerCol + x + 1);
		}
	}

	glBegin (GL_TRIANGLES);
	for (uint i = 0 ; i < vertices.size() ; ++i)
	{
		std::vector< double > p = splinePoints[vertices[i]];
		double x = p[0];
		double y = p[1];
		double z = p[2];
		glVertex3f(x,y,z);
	}
	glEnd();
}

void Surface::draw()
{
	std::vector< std::vector< double > > givenPoints;
	int countPoints = m_treeWidget->GetChildrenCount(m_tPointId, true);

	wxTreeItemId id, childid;
	wxTreeItemIdValue cookie = 0;
	for (int i = 0 ; i < countPoints ; ++i)
	{
		id = m_treeWidget->GetNextChild(m_tPointId, cookie);
		Point *point = (Point*)((MyTreeItemData*)m_treeWidget->GetItemData(id))->getData();

		std::vector< double > p;
		p.push_back(point->getCenter().s.X);
		p.push_back(point->getCenter().s.Y);
		p.push_back(point->getCenter().s.Z);
		givenPoints.push_back(p);
	}

	execute(givenPoints);
}
