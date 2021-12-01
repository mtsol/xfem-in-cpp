/*
 * linearTetrahedron.cpp
 *
 *  Created on: March, 28 2013
 *      Author: Christoph Paulus
 */

// includes
#include "LinearTetrahedron.h"

// function definitions
LinearTetrahedron::LinearTetrahedron(void)
{
}


LinearTetrahedron::~LinearTetrahedron()
{
}


void LinearTetrahedron::init(int tetrahedronId,Ints_t connect,int neighboursOfFaces[][2],Points_t &points,double lambda,double mu,int numberOfIntegrationPoints)
{
	this->m_tetrahedronId = tetrahedronId;
	m_connect = connect;
	m_stiffnessMatrix.resize(3*connect.size(),3*connect.size());
	for(int faceId=0;faceId<4;faceId++)
	{
		for(int i=0;i<2;i++)	// i=0 is the element and i=1 is the local face Id
			m_neighboursOfFaces[faceId][i] = neighboursOfFaces[faceId][i];
		m_triangleInFaceId.push_back(faceId);
	}
	this->m_points = points;
	getShapeFunctionsOf_m_points();
	this->m_numberOfElementNodes = 4;
	this->m_pointsDisplaced.resize(4);

	int nodesOnTetrahedronFaces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
	for(int faceId=0;faceId<4;faceId++)
	{
		ConnectTriangle_t currentConnectTriangle;
		for(int pointId=0;pointId<3;pointId++)
		{
			currentConnectTriangle[pointId] = nodesOnTetrahedronFaces[faceId][pointId];
		}
		this->m_connectTriangle.push_back(currentConnectTriangle);
	}

	getGlobalToLocal();

	// for(int i=0;i<NUMBER_OF_ELEMENT_NODES*DIM;i++)
	for(int i=0;i<connect.size()*3;i++)
	{
		// for(int j=0;j<NUMBER_OF_ELEMENT_NODES*DIM;j++)
		for(int j=0;j<connect.size()*3;j++)
		{
			// this->m_stiffnessMatrix[i][j] = 0;
			this->m_stiffnessMatrix(i,j) = 0;
		}
	}

	// get the gauss integration points in the unit element
	getGaussIntegrationPoints(this->m_GaussIntPoints, this->m_GaussIntWeights);

	getNumberOfEquallyDistributedIntegrationPoints(numberOfIntegrationPoints);
	getEquallyDistributedIntegrationPoints();

	this->m_lambda = lambda;
	this->m_mu = mu;
}


void LinearTetrahedron::getGlobalToLocal()
{
	Matrix_3_3_t inverseOfGlobalToLocal;
	for(int i=0;i<DIM;i++)
    {
		inverseOfGlobalToLocal(i,0) = (this->m_points[0])[i]-(this->m_points[2])[i];
		inverseOfGlobalToLocal(i,1) = (this->m_points[1])[i]-(this->m_points[2])[i];
		inverseOfGlobalToLocal(i,2) = (this->m_points[3])[i]-(this->m_points[2])[i];
    }
	this->m_globalToLocal = inverseOfGlobalToLocal.inverse();
}


void LinearTetrahedron::getShapeFunctionsOf_m_points()
{
	for(Points_t::iterator iter=this->m_points.begin();iter!=this->m_points.end();iter++)
	{
		Point_t currentPoint = *iter;
		ShapeFct_t currentShapeFct;
		getShapeFunction(currentPoint,currentShapeFct);
		this->m_shapeFunctions.push_back(currentShapeFct);
	}
}


void LinearTetrahedron::getShapeFunctionDependentOnLocalCoords(Point_t xi, ShapeFct_t &ShapeFunction)
{
	ShapeFunction.push_back(xi[0]);
	ShapeFunction.push_back(xi[1]);
	ShapeFunction.push_back(1-xi[0]-xi[1]-xi[2]);
	ShapeFunction.push_back(xi[2]);
}


void LinearTetrahedron::transformLocalToGlobalCoordinatesByUsingLinearShapeFunctions(Points_t &TetCoords, Point_t xi, Point_t &x)
{
	ShapeFct_t SF;
	getShapeFunctionDependentOnLocalCoords(xi,SF);


	for(int j=0;j<DIM;j++)
	{
		x[j] = 0;
        for(int i=0;i<NUMBER_OF_ELEMENT_NODES;i++)
            x[j] += SF[i]*(TetCoords[i])[j];
	}
}


void LinearTetrahedron::transformGlobalToLocalCoordinatesByUsingLinearShapeFunctions(Point_t x, Point_t &xi)
{
	for(int i=0;i<DIM;i++)
	{
		xi[i] = 0;
        for(int j=0;j<DIM;j++)
            xi[i] += this->m_globalToLocal(i,j)*(x[j]-(this->m_points[2])[j]);
	}
}


void LinearTetrahedron::getShapeFunction(Point_t x, ShapeFct_t &ShapeFunction)
{
	Point_t xi;
	transformGlobalToLocalCoordinatesByUsingLinearShapeFunctions(x,xi);

	getShapeFunctionDependentOnLocalCoords(xi,ShapeFunction);
}


void LinearTetrahedron::getShapeFunctionLocalDerivatives(Point_t xi, ShapeFctDeriv_t &ShapeFunctionLocalDerivatives)
{
	for(int i=0;i<NUMBER_OF_ELEMENT_NODES;i++)
	{
		Point_t currentLocalDerivativeAtElementNode;
		for(int j=0;j<DIM;j++)
		{
			currentLocalDerivativeAtElementNode[j] = (i==0)*(j==0)+(i==1)*(j==1)-(i==2)+(i==3)*(j==2);
		}
		ShapeFunctionLocalDerivatives.push_back(currentLocalDerivativeAtElementNode);
	}
	/* this information is obtained by XFEMForceField<DataTypes>::InitializeShapeFunctionDerivatives
	ShapeFunctionLocalDerivatives[0][0] = 1;
	ShapeFunctionLocalDerivatives[0][1] = 0;
	ShapeFunctionLocalDerivatives[0][2] = 0;

	ShapeFunctionLocalDerivatives[1][0] = 0;
	ShapeFunctionLocalDerivatives[1][1] = 1;
	ShapeFunctionLocalDerivatives[1][2] = 0;

	ShapeFunctionLocalDerivatives[2][0] = -1;
	ShapeFunctionLocalDerivatives[2][1] = -1;
	ShapeFunctionLocalDerivatives[2][2] = -1;

	ShapeFunctionLocalDerivatives[3][0] = 0;
	ShapeFunctionLocalDerivatives[3][1] = 0;
	ShapeFunctionLocalDerivatives[3][2] = 1;*/
}


void LinearTetrahedron::getJacobian(ShapeFctDeriv_t &ShapeFunctionLocalDerivatives, Matrix_3_3_t &jacobian)	// does not work for XFEM since the asymptotic crack tip functions are nonlinear - I think!
{
	// calculation of the jacobian of the transformation between the unit and element tetrahedron

	for(int m=0;m<DIM;m++)
	{
		for(int n=0;n<DIM;n++)
		{
			for(int l=0;l<4;l++) // we only want to go up to the four points that are used for the linear transformation
			{
				// jacobian[m][n] += (this->m_points[l])[m]*(ShapeFunctionLocalDerivatives[l])[n];
				jacobian(m,n) += (this->m_points[l])[m]*(ShapeFunctionLocalDerivatives[l])[n];
			}
		}
	}
}


double LinearTetrahedron::getJacobianDeterminant(Point_t xi)
{
	ShapeFctDeriv_t SFLDs;
	LinearTetrahedron::getShapeFunctionLocalDerivatives(xi,SFLDs);

	Matrix_3_3_t jacobian;
	LinearTetrahedron::getJacobian(SFLDs,jacobian);
	//double jacobianDeterminant = sofa::defaulttype::determinant(jacobian);
    double jacobianDeterminant = jacobian.determinant();

	return jacobianDeterminant;
}


void LinearTetrahedron::getShapeFunctionDerivatives(Point_t xi, ShapeFctDeriv_t &ShapeFunctionDerivative)
{
	// get the local shape function derivative
	ShapeFctDeriv_t SFLDs;
	getShapeFunctionLocalDerivatives(xi,SFLDs);

	/* declare, calculate and invert the jacobian matrix at the natural coordinates r,s,t
	(in case of linear shape functions, the jacobian is constant and therefore independent of r,s,t)*/
	Matrix_3_3_t jacobian;
	getJacobian(SFLDs,jacobian);
	Matrix_3_3_t inverseJacobian = jacobian.inverse();
	// sofa::defaulttype::invertMatrix(inverseJacobian,jacobian);

	// calculate the derivative of the shape function by multiplying the inverse jacobian with the local derivatives of the shape function
	for(int i=0;i<NUMBER_OF_ELEMENT_NODES;i++)
	{
		Point_t currentDerivativeAtElementNode;
		for(int j=0;j<DIM;j++)
		{
			currentDerivativeAtElementNode[j] = 0;
			for(int l=0;l<DIM;l++)
			{
				// currentDerivativeAtElementNode[j] += (SFLDs[i])[l]*inverseJacobian[l][j];
				currentDerivativeAtElementNode[j] += (SFLDs[i])[l]*inverseJacobian(l,j);
			}
		}
		ShapeFunctionDerivative.push_back(currentDerivativeAtElementNode);
	}
}


void LinearTetrahedron::getShapeFunctionDerivativesDependentOnIntegrationPointNumber
	(int integrationPointNumber, ShapeFctDeriv_t &ShapeFunctionDerivative)
{
	// get the current integration point
	Point_t xi;
	for(int j=0;j<DIM;j++)
	{
		xi[j] = this->m_localIntegrationPoints[integrationPointNumber][j];
	}

	// get the shape function derivatives
	getShapeFunctionDerivatives(xi,ShapeFunctionDerivative);
}


void LinearTetrahedron::getNumberOfEquallyDistributedIntegrationPoints(int numberOfIntegrationPointsWanted)
{
	/* since the integration points have a certain distribution over the
	 * tetrahedron they can't be chosen arbitrarily; the functions determines
	 * the smallest possible number of integration points more than
	 * nIntPointsWanted, returns the actual number of integration points
	 * nIntPoints and the number of segments in each direction for further
	 * information I am refering to my Diploma thesis:
	 * "Simulation beliebiger Schnitte in Weichgewebe mit der Extended
	 * Finite Element Method" page 68; follows notes on page 102,103,211-215*/

    // initialize the variables
    this->m_numberOfIntegrationPoints = 0;
    this->m_numberOfSegments = 0;

	// increase number of cubes until we receive more integration points than needed
    while(numberOfIntegrationPointsWanted>=this->m_numberOfIntegrationPoints)
    {
		this->m_numberOfSegments++;
        this->m_numberOfIntegrationPoints = ((((this->m_numberOfSegments-2) + 3)*(this->m_numberOfSegments-2) + 2)*(this->m_numberOfSegments-2)
									+ this->m_numberOfSegments*(this->m_numberOfSegments-1)*5/2
									+ (this->m_numberOfSegments+1)*this->m_numberOfSegments/2)
                                    *NUMBER_OF_INTEGRATION_POINTS_TET;
	}
}


void LinearTetrahedron::getGaussIntegrationPoints(Points_t &integrationPoints, Doubles_t &integrationWeights)
{
	// follows XFEMForceField<DataTypes>::InitializeIntegrationPointsAndWeights()
	double alpha1 = 0.58541020;
	double beta1 = 0.13819660;

	double alpha2 = 0.785714285714286;
	double beta2 = 0.071428571428571;

	double alpha3 = 0.399403576166799;
	double  beta3 = 0.100596423833201;

	double intPoints[NUMBER_OF_INTEGRATION_POINTS_TET*DIM];
	double intWeights[NUMBER_OF_INTEGRATION_POINTS_TET];

	switch(NUMBER_OF_INTEGRATION_POINTS_TET)
	{
		case 1:
		  intPoints[0] = 0.25;
		  intPoints[1] = 0.25;
		  intPoints[2] = 0.25;
		  intWeights[0] = (double)1.0/(double)6.0;
		  break;
		case 4:
		  intPoints[0] = alpha1;
		  intPoints[1] = beta1;
		  intPoints[2] = beta1;
		  intPoints[3] = beta1;
		  intPoints[4] = alpha1;
		  intPoints[5] = beta1;
		  intPoints[6] = beta1;
		  intPoints[7] = beta1;
		  intPoints[8] = alpha1;
		  intPoints[9] = beta1;
		  intPoints[10] = beta1;
		  intPoints[11] = beta1;
		  intWeights[0] = 0.25;
		  intWeights[1] = 0.25;
		  intWeights[2] = 0.25;
		  intWeights[3] = 0.25;
		  break;
		case 5:
		  intPoints[0] = 0.25;
		  intPoints[1] = 0.25;
		  intPoints[2] = 0.25;
		  intPoints[3] = 0.5;
		  intPoints[4] = 1.0/6.0;
		  intPoints[5] = 1.0/6.0;
		  intPoints[6] = 1.0/6.0;
		  intPoints[7] = 0.5;
		  intPoints[8] = 1.0/6.0;
		  intPoints[9] = 1.0/6.0;
		  intPoints[10] = 1.0/6.0;
		  intPoints[11] = 0.5;
		  intPoints[12] = 1.0/6.0;
		  intPoints[13] = 1.0/6.0;
		  intPoints[14] = 1.0/6.0;
		  intWeights[0] = -0.8;
		  intWeights[1] = 9.0/20.0;
		  intWeights[2] = 9.0/20.0;
		  intWeights[3] = 9.0/20.0;
		  intWeights[4] = 9.0/20.0;
		  break;
		case 11:
		  intPoints[0] = 0.25;
		  intPoints[1] = 0.25;
		  intPoints[2] = 0.25;
		  intPoints[3] = alpha2;
		  intPoints[4] = beta2;
		  intPoints[5] = beta2;
		  intPoints[6] = beta2;
		  intPoints[7] = alpha2;
		  intPoints[8] = beta2;
		  intPoints[9] = beta2;
		  intPoints[10] = beta2;
		  intPoints[11] = alpha2;
		  intPoints[12] = beta2;
		  intPoints[13] = beta2;
		  intPoints[14] = beta2;
		  intPoints[15] = alpha3;
		  intPoints[16] = alpha3;
		  intPoints[17] = beta3;
		  intPoints[18] = alpha3;
		  intPoints[19] = beta3;
		  intPoints[20] = alpha3;
		  intPoints[21] = alpha3;
		  intPoints[22] = beta3;
		  intPoints[23] = beta3;
		  intPoints[24] = beta3;
		  intPoints[25] = alpha3;
		  intPoints[26] = alpha3;
		  intPoints[27] = beta3;
		  intPoints[28] = alpha3;
		  intPoints[29] = beta3;
		  intPoints[30] = beta3;
		  intPoints[31] = beta3;
		  intPoints[32] = alpha3;
		  intWeights[0] = -0.013133333333336;
		  intWeights[1] = 0.007622222222222;
		  intWeights[2] = 0.007622222222222;
		  intWeights[3] = 0.007622222222222;
		  intWeights[4] = 0.007622222222222;
		  intWeights[5] = 0.24888888888889;
		  intWeights[6] = 0.24888888888889;
		  intWeights[7] = 0.24888888888889;
		  intWeights[8] = 0.24888888888889;
		  intWeights[9] = 0.24888888888889;
		  intWeights[10] = 0.24888888888889;
		  break;
	}

	for(int intPointId=0;intPointId<NUMBER_OF_INTEGRATION_POINTS_TET;intPointId++)
	{
		Point_t currentintPointId;
		for(int i=0;i<3;i++)
			currentintPointId[i] = intPoints[intPointId*3+i];
		integrationPoints.push_back(currentintPointId);
		integrationWeights.push_back(intWeights[intPointId]);
	}
}


void LinearTetrahedron::getGlobalIntegrationPointsIn_Cube_CubeWithoutTet_Tet(
		int walkX, int walkY, int walkZ,/*sofa::helper::fixed_array<Point_t,8> &CoordsReferenceCube,*/
		Eigen::Matrix<double,8,3> &CoordsReferenceCube,
		int tetrahedronConnectOfCube8[][4], int numberOfTetrahedra)//)
{
	Point_t currentCoords;
	Point_t currentGaussIntegrationPoint;
	Point_t currentLocalIntegrationPoint;
	Point_t currentGlobalIntegrationPoint;
	currentCoords[0] = (double)walkX/(double)this->m_numberOfSegments;
	currentCoords[1] = (double)walkY/(double)this->m_numberOfSegments;
	currentCoords[2] = (double)walkZ/(double)this->m_numberOfSegments;
	for(int tetrahedronNumberInCube=0;tetrahedronNumberInCube<numberOfTetrahedra;tetrahedronNumberInCube++)
	{
		Points_t coordsCurrentTetrahedronInCube;
		for(int j=0;j<NUMBER_OF_ELEMENT_NODES;j++)
		{
			Point_t currentTetrahedronPoint;
			for(int i=0;i<DIM;i++)
			{
				currentTetrahedronPoint[i] = currentCoords[i] +
						// CoordsReferenceCube[tetrahedronConnectOfCube8[tetrahedronNumberInCube][j]][i]/this->m_numberOfSegments;
						CoordsReferenceCube(tetrahedronConnectOfCube8[tetrahedronNumberInCube][j],i)/this->m_numberOfSegments;
			}
			coordsCurrentTetrahedronInCube.push_back(currentTetrahedronPoint);
		}
		for(int integrationPointNumberTet=0;integrationPointNumberTet<NUMBER_OF_INTEGRATION_POINTS_TET;integrationPointNumberTet++)
		{
			currentGaussIntegrationPoint = this->m_GaussIntPoints[integrationPointNumberTet];
			transformLocalToGlobalCoordinatesByUsingLinearShapeFunctions(coordsCurrentTetrahedronInCube, currentGaussIntegrationPoint, currentLocalIntegrationPoint);
			this->m_localIntegrationPoints.push_back(currentLocalIntegrationPoint);
			transformLocalToGlobalCoordinatesByUsingLinearShapeFunctions(this->m_points, currentLocalIntegrationPoint, currentGlobalIntegrationPoint);
			this->m_globalIntegrationPoints.push_back(currentGlobalIntegrationPoint);
			this->m_integrationWeights.push_back(this->m_GaussIntWeights[integrationPointNumberTet]/this->m_numberOfIntegrationPoints);
		}
	}
}


void LinearTetrahedron::getEquallyDistributedIntegrationPoints()
{
	// get the unit cube - [0,1]^3 - and the separation of the unit cube in tetrahedra
	int tetrahedronConnectOfCube8[6][4] =
        {{0,1,3,4} , {1,2,3,4} , {3,2,7,4} , {1,4,5,2} , {5,4,7,2} , {2,6,7,5}};
	Eigen::Matrix<double,8,3> CoordsReferenceCube;

	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
			{
			for(int k=0;k<2;k++)
                {
					CoordsReferenceCube(i*4+j*2+k,0) = (double)(k!=j);
					CoordsReferenceCube(i*4+j*2+k,1) = (double)(j==1);
					CoordsReferenceCube(i*4+j*2+k,2) = (double)(i==1);
				}
			}
	}

    // calculation of the integration points ...
    int numberOfTetrahedra;

	for(int walkZ=0;walkZ<this->m_numberOfSegments;walkZ++)
	{
        for(int walkY=0;walkY<this->m_numberOfSegments-walkZ;walkY++)
	    {
	        // ... the points that come from the cubes
			for(int walkX=0;walkX<this->m_numberOfSegments-walkY-walkZ;walkX++)
			{
				if(walkX<this->m_numberOfSegments-2-walkY-walkZ)
				{
					numberOfTetrahedra = 6;
				}
				else if (walkX==this->m_numberOfSegments-2-walkY-walkZ)
				{
					numberOfTetrahedra = 5;
				}
				else if (walkX==this->m_numberOfSegments-1-walkY-walkZ)
				{
					numberOfTetrahedra = 1;
                }

                getGlobalIntegrationPointsIn_Cube_CubeWithoutTet_Tet(walkX,walkY,walkZ,CoordsReferenceCube, tetrahedronConnectOfCube8,numberOfTetrahedra);
			}
	    }
	}
}


// void LinearTetrahedron::addTetrahedronVertexesToIntegrationPoints()
// {
//	// adjust the integration weights
//	double help=0;
//	double help2=0;
//	for(int integrationPointNumber=0;integrationPointNumber<m_numberOfIntegrationPoints;integrationPointNumber++)
//	{
//		help+=m_integrationWeights[integrationPointNumber];
//		m_integrationWeights[integrationPointNumber] *= (double)(m_numberOfIntegrationPoints)/(double)(m_numberOfIntegrationPoints+4);
//		help2+=m_integrationWeights[integrationPointNumber];
//	}
////	std::cout << help << " " << help2 << " " << m_numberOfIntegrationPoints << " ";
//	for(int additionalIntegrationPointNumber=0;additionalIntegrationPointNumber<4;additionalIntegrationPointNumber++)
//	{
//		m_integrationWeights.push_back(1/(m_numberOfIntegrationPoints+4)/6);
//	}
//
//	// add new points
//	double localTetCoords[4][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,0.0},{0.0,0.0,1.0}};
//	for(int additionalIntegrationPoint=0;additionalIntegrationPoint<4;additionalIntegrationPoint++)
//	{
//		Point_t currentTetrahedronVertex;
//		Point_t currentLocalCoordsTetrahedronVertex;
//
//		for(int i=0;i<DIM;i++)
//		{
//			currentLocalCoordsTetrahedronVertex[i] = localTetCoords[additionalIntegrationPoint][i];
//			currentTetrahedronVertex[i] = (m_points[additionalIntegrationPoint])[i];
//		}
//
//		m_localIntegrationPoints.push_back(currentLocalCoordsTetrahedronVertex);
//		m_globalIntegrationPoints.push_back(currentTetrahedronVertex);
//	}
//
//	// adjust the number of integration points
//	m_numberOfIntegrationPoints += 4;
//
//	help=0;
//	for(int integrationPointNumber=0;integrationPointNumber<m_numberOfIntegrationPoints;integrationPointNumber++)
//		help+=m_integrationWeights[integrationPointNumber];
//
////	std::cout << help << " " << m_numberOfIntegrationPoints<<endl << endl;
// }


void LinearTetrahedron::stiffnessBlock()
{
	double jacobianDeterminant, integrationWeight;
	Point_t xi;
	StiffnessMatrix_t stiffnessBlockTemp;

	for(int integrationPointNumber=0;integrationPointNumber<this->m_numberOfIntegrationPoints;integrationPointNumber++)
	{
		for(int j=0;j<DIM;j++)
		{
			xi[j] = this->m_localIntegrationPoints[integrationPointNumber][j];
		}
		ShapeFctDeriv_t SFDs;
		getShapeFunctionDerivatives(xi,SFDs);

		stiffnessBlockDependentOnShapeFunctionDerivatives(SFDs,stiffnessBlockTemp);
		jacobianDeterminant = getJacobianDeterminant(xi);
		integrationWeight = this->m_integrationWeights[integrationPointNumber];
		for(int i=0;i<this->m_numberOfElementNodes*DIM;i++)
		{
			for(int j=0;j<this->m_numberOfElementNodes*DIM;j++)
			{
				// this->m_stiffnessMatrix[i][j] += stiffnessBlockTemp[i][j]*jacobianDeterminant*integrationWeight;
				this->m_stiffnessMatrix(i,j) += stiffnessBlockTemp(i,j)*jacobianDeterminant*integrationWeight;
			}
		}
	}
}


void LinearTetrahedron::stiffnessBlockDependentOnShapeFunctionDerivatives(ShapeFctDeriv_t &SFDs,StiffnessMatrix_t &stiffnessBlock)
{
	stiffnessBlock.resize(m_numberOfElementNodes*3,m_numberOfElementNodes*3);
	int mn, ij;
	for(int m=0;m<this->m_numberOfElementNodes;m++)
	{
		for(int n=0;n<DIM;n++)
        {
			for(int i=0;i<this->m_numberOfElementNodes;i++)
			{
				for(int j=0;j<DIM;j++)
				{
					mn = 3*m+n;
					ij = 3*i+j;

//					stiffnessBlock[ij][mn] = 0;
					// stiffnessBlock[ij][mn] = this->m_lambda*(SFDs[m])[n]*(SFDs[i])[j]+this->m_mu*(SFDs[m])[j]*(SFDs[i])[n];
					stiffnessBlock(ij,mn) = this->m_lambda*(SFDs[m])[n]*(SFDs[i])[j]+this->m_mu*(SFDs[m])[j]*(SFDs[i])[n];
					if(j==n)
					{
						for(int k=0;k<DIM;k++)
						{
							// stiffnessBlock[ij][mn] += this->m_mu*(SFDs[m])[k]*(SFDs[i])[k];
							stiffnessBlock(ij,mn) += this->m_mu*(SFDs[m])[k]*(SFDs[i])[k];
						}
					}
				}
			}
		}
	}
}


void LinearTetrahedron::stiffnessBlock2()
{
	double jacobianDeterminant, integrationWeight;
    //Point_t xi;
	StiffnessMatrix_t stiffnessBlockTemp;

	for(int integrationPointNumber=0;integrationPointNumber<this->m_numberOfIntegrationPoints;integrationPointNumber++)
	{
        /*for(int j=0;j<DIM;j++)
            xi[j] = this->m_localIntegrationPoints[integrationPointNumber][j];*/
		ShapeFctDeriv_t SFDs;
		getShapeFunctionDerivativesDependentOnIntegrationPointNumber(integrationPointNumber,SFDs);

		stiffnessBlockDependentOnShapeFunctionDerivatives(SFDs,stiffnessBlockTemp);
        jacobianDeterminant = getJacobianDeterminant();     // since we work with linear tetrahedra, the jacobian determinant is independant on its
        integrationWeight = this->m_integrationWeights[integrationPointNumber];
		for(int i=0;i<this->m_numberOfElementNodes*DIM;i++)
		{
			for(int j=0;j<this->m_numberOfElementNodes*DIM;j++)
			{
                // this->m_stiffnessMatrix[i][j] += stiffnessBlockTemp[i][j]*jacobianDeterminant*integrationWeight;
                this->m_stiffnessMatrix(i,j) += stiffnessBlockTemp(i,j)*jacobianDeterminant*integrationWeight;
                assert(!isnan(this->m_stiffnessMatrix(i,j)));
			}
		}
	}
}


void LinearTetrahedron::stiffnessBlock3()
{
	double jacobianDeterminant, integrationWeight;
	Point_t xi;

	for(int integrationPointNumber=0;integrationPointNumber<this->m_numberOfIntegrationPoints;integrationPointNumber++)
	{
		for(int j=0;j<DIM;j++)
		{
			xi[j] = this->m_localIntegrationPoints[integrationPointNumber][j];
		}
		ShapeFctDeriv_t SFDs;
		getShapeFunctionDerivativesDependentOnIntegrationPointNumber(integrationPointNumber,SFDs);
		jacobianDeterminant = getJacobianDeterminant(xi);

		integrationWeight = this->m_integrationWeights[integrationPointNumber];
		int mn, ij;
		for(int m=0;m<this->m_numberOfElementNodes;m++)
		{
			for(int n=0;n<DIM;n++)
			{
				for(int i=0;i<this->m_numberOfElementNodes;i++)
				{
					for(int j=0;j<DIM;j++)
					{
						mn = 3*m+n;
						ij = 3*i+j;

						// this->m_stiffnessMatrix[ij][mn] += (this->m_lambda*(SFDs[m])[n]*(SFDs[i])[j]+this->m_mu*(SFDs[m])[j]*(SFDs[i])[n])*jacobianDeterminant*integrationWeight;
						this->m_stiffnessMatrix(ij,mn) += (this->m_lambda*(SFDs[m])[n]*(SFDs[i])[j]+this->m_mu*(SFDs[m])[j]*(SFDs[i])[n])*jacobianDeterminant*integrationWeight;
						if(j==n)
						{
							for(int k=0;k<DIM;k++)
							{
								// this->m_stiffnessMatrix[ij][mn] += this->m_mu*(SFDs[m])[k]*(SFDs[i])[k]*jacobianDeterminant*integrationWeight;
								this->m_stiffnessMatrix(ij,mn) += this->m_mu*(SFDs[m])[k]*(SFDs[i])[k]*jacobianDeterminant*integrationWeight;
							}
						}
					}
				}
			}
		}
	}
}


void LinearTetrahedron::writeDisplacementInTetrahedronData(Points_t &displacement)
{
	this->m_displacement = displacement;

	int numberOfElementNodes = std::distance(displacement.begin(),displacement.end());
	if(numberOfElementNodes!=this->m_numberOfElementNodes)
		std::cout << "The displacement vector has another size then the number of nodes on the element " << this->m_tetrahedronId << endl;
}


LinearTetrahedron::Point_t LinearTetrahedron::displacePoint(ShapeFct_t &SF, Point_t &initialPosition)
{
    Point_t x_displaced = initialPosition;
    for(int j=0;j<DIM;j++)
        for(int i=0;i<this->m_numberOfElementNodes;i++)
            x_displaced[j] += SF[i]*(this->m_displacement[i])[j];
    return x_displaced;
}


void LinearTetrahedron::displaceTetrahedronPoints()
{
    for(unsigned i = 0 ; i < m_shapeFunctions.size() ; ++i)
        this->m_pointsDisplaced[i] = displacePoint(m_shapeFunctions[i],m_points[i]);
}


void LinearTetrahedron::displaceReferencePoints()
{
    this->m_displacedReferencePointsInTet.resize(m_referencePointsInTet.size());
    for(unsigned i = 0 ; i < m_referencePointShapeFunctions.size() ; ++i)
        this->m_displacedReferencePointsInTet[i] = displacePoint(m_referencePointShapeFunctions[i],m_referencePointsInTet[i]);
}


void LinearTetrahedron::PointsToVTKPointContainer(Points_t &points, vtkPoints* pointContainerVTK, int& numberOfPoints)
{
	numberOfPoints = std::distance(points.begin(),points.end());
	double currentPointDouble[3];

	for(Points_t::iterator pointIterator=points.begin(); pointIterator!=points.end();++pointIterator)
	{
		Point_t currentPoint = *pointIterator;
		for(int i=0;i<3;i++)
		{
			currentPointDouble[i] = currentPoint[i];
		}
		pointContainerVTK->InsertNextPoint(currentPointDouble);
	}
}


void LinearTetrahedron::cellsToVTKCellContainer(Ints_t &cellIds, int numberOfAlreadyExistingPoints,
		vtkCellArray* cellContainerVTK)
{
	int numberOfCellNodes = 3; // this can also be provided as an input and throught that
	vtkIdType* currentCell;
	currentCell = new vtkIdType[numberOfCellNodes];
	for(Ints_t::iterator cellIter=cellIds.begin();cellIter!=cellIds.end();cellIter++)
	{
		int cellId = *cellIter;
		if((m_neighboursOfFaces[this->m_triangleInFaceId[cellId]][0]==-1) || (this->m_triangleInFaceId[cellId] == -1))
		{
			ConnectTriangle_t currentConnectCell = this->m_connectTriangle[cellId];
			for(int nodeId=0; nodeId<numberOfCellNodes; nodeId++)
			{
				currentCell[nodeId] = currentConnectCell[nodeId]+numberOfAlreadyExistingPoints;
			}

			cellContainerVTK->InsertNextCell(numberOfCellNodes, currentCell);
		}
	}
}


void LinearTetrahedron::tetrahedronToVTKContainer(int plotDisplaced, Ints_t &cellIds, int& numberOfAlreadyExistingPoints, // this is in- and output
		vtkPoints* pointContainerVTK,vtkCellArray* cellContainerVTK)
{
	Points_t points;
	if(plotDisplaced)
	{
		points = this->m_pointsDisplaced;
	}
	else
	{
		points = this->m_points;
	}

	int currentNumberOfPoints=0;
	PointsToVTKPointContainer(points,pointContainerVTK, currentNumberOfPoints);
	cellsToVTKCellContainer(cellIds,numberOfAlreadyExistingPoints,cellContainerVTK);

	numberOfAlreadyExistingPoints+=currentNumberOfPoints;
}


void LinearTetrahedron::getCellIdsForOutput(Ints_t &cellIds)
{
	int numberOfElements = std::distance(m_connectTriangle.begin(),m_connectTriangle.end());
	for(int cellId=0;cellId<numberOfElements;cellId++)
	{
		cellIds.push_back(cellId);
	}
}


// the following should be DELETED - I just use it for now!
void LinearTetrahedron::outputFromVTKContainersToVTKFile(int plotDisplaced)
{
	vtkPoints* pointContainerVTK = vtkPoints::New();
	vtkCellArray* triangleContainerVTK = vtkCellArray::New();

	Ints_t cellIds;
	getCellIdsForOutput(cellIds);
	int numberOfAlreadyExistingPoints = 0;
	tetrahedronToVTKContainer(plotDisplaced,cellIds,numberOfAlreadyExistingPoints,
			pointContainerVTK,triangleContainerVTK);

	// define input and output filename
	std::string OutputFolder = "/org/share/home/paulus/SOFA_Results/LinearTetrahedra/";

	std::string filename = OutputFolder;
	std::stringstream ss;
	ss.str(std::string());
	ss << plotDisplaced;
	filename += ss.str();
	filename += "LinearTet_";
	ss.str(std::string());
	ss << this->m_tetrahedronId;
	filename += ss.str();
	filename += ".vtk";

	vtkUnstructuredGrid* myGridTopoDebug = vtkUnstructuredGrid::New();
	myGridTopoDebug->SetPoints(pointContainerVTK);
	myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK);

	vtkUnstructuredGridWriter* writerTopoDebug = vtkUnstructuredGridWriter::New();

	writerTopoDebug->SetInput(myGridTopoDebug);
	writerTopoDebug->SetFileName(filename.c_str());
	writerTopoDebug->Write();

	pointContainerVTK->Delete();
	triangleContainerVTK->Delete();
	myGridTopoDebug->Delete();
	writerTopoDebug->Delete();
}


template<class Element>
void LinearTetrahedron::writeVectorInRow(std::vector<Element> & vector, std::ofstream & myfile) {
    for(unsigned i=0;i<vector.size();++i)
        myfile << vector[i] << " ";
    myfile << "\n";
}


template<class Element>
void LinearTetrahedron::writeMatrixInFile(std::vector<  std::vector<Element>   > & matrix, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::trunc);
    for(unsigned i=0;i<matrix.size();++i)
        writeVectorInRow(matrix[i], myfile);
    myfile.close();
}


void LinearTetrahedron::writeShapeFunctionDerivativeInFile(ShapeFctDeriv_t & dSF, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::trunc);
    for ( unsigned j = 0 ; j < DIM ; ++j )
    {
        for(unsigned i=0;i<dSF.size();++i)
            myfile << (dSF[i])[j] << " ";
        myfile << "\n";
    }
    myfile.close();
}

void LinearTetrahedron::writeMatrixInFile(LinearTetrahedron::StiffnessMatrix_t & matrix, int sizeOfMatrix, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::trunc);
    for(unsigned i=0;i<sizeOfMatrix;++i)
    {
        for(unsigned j=0;j<sizeOfMatrix;++j)
            myfile << matrix(i,j) << " ";
        myfile << "\n";
    }
    myfile.close();
}


template<class Element>
void LinearTetrahedron::writeVectorInRowOfFile(std::vector<Element> & vector, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::app);
    writeVectorInRow(vector, myfile);
    myfile.close();
}
