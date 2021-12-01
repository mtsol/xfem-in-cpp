/*
 * XFemLinearTetrahedron.cpp
 *
 *  Created on: July 12, 2013
 *      Author: Christoph Paulus
 */

// includes
#include "XFemLinearTetrahedron.h"
#include <limits>

// function definitions
XFemLinearTetrahedron::XFemLinearTetrahedron(void)
{
}


XFemLinearTetrahedron::~XFemLinearTetrahedron()
{
}


void XFemLinearTetrahedron::init(int tetrahedronId, Ints_t connect, int neighboursOfFaces[][2], Points_t &points, Ints_t pointIsAbove, Ints_t &signEnrichedVertexIds, Ints_t &branchEnrichedVertexIds,
		CutModelTopology_CGAL *cutModel,ConnectTriangles_t &connectTriangle, Ints_t &triangleInFaceId, Ints_t &triangleAbove,
		double lambda,double mu,int numberOfIntegrationPoints)
{
	this->m_connect = connect;
	m_stiffnessMatrix.resize(connect.size()*3,connect.size()*3);
	this->m_tetrahedronId = tetrahedronId;
	for(int faceId=0;faceId<4;faceId++)
			for(int i=0;i<2;i++)	// i=0 is the element and i=1 is the local face Id
				m_neighboursOfFaces[faceId][i] = neighboursOfFaces[faceId][i];
	// initialize information about the points
	this->m_points = points;

	this->m_pointIsAbove = pointIsAbove;
	for(int pointId=0;pointId<4;pointId++)
	{
		if(this->m_pointIsAbove[pointId]==-1)
		  std::cout << "WATCH OUT, in Element " << m_tetrahedronId << 
		    " the vertex with Id " << pointId << " lies ON the cut!" << endl;
	}
	int numberOfPointsIncludingTheCutPoints = std::distance(points.begin(),points.end());
	this->m_points_r.resize(numberOfPointsIncludingTheCutPoints);
	this->m_points_phi.resize(numberOfPointsIncludingTheCutPoints);
	this->m_shapeFunctions.resize(numberOfPointsIncludingTheCutPoints);
	this->m_pointsDisplaced.resize(numberOfPointsIncludingTheCutPoints);

	getGlobalToLocal();

	for(int i=0;i<connect.size()*3;i++)
	{
		for(int j=0;j<connect.size()*3;j++)
		{
			this->m_stiffnessMatrix(i,j) = 0.0;
		}
	}

	this->m_cutModel = cutModel;
	getInformationAboutPartialAndCompleteCut();

	if(this->m_tetrahedronIsSeparated)
		this->m_signEnrichedVertexIds = signEnrichedVertexIds;
	this->m_numberOfElementNodes =  4+this->m_tetrahedronIsSeparated*std::distance(signEnrichedVertexIds.begin(),signEnrichedVertexIds.end())
		+ 4*std::distance(branchEnrichedVertexIds.begin(),branchEnrichedVertexIds.end());
	this->m_branchEnrichedVertexIds = branchEnrichedVertexIds;

	getShapeFunctionsOf_m_points();

	this->m_connectTriangle = connectTriangle;
	this->m_triangleInFaceId = triangleInFaceId;
	this->m_triangleAbove = triangleAbove;

	// get the gauss integration points in the unit element
	getGaussIntegrationPoints(this->m_GaussIntPoints, this->m_GaussIntWeights);

	// get the integration points
    getNumberOfEquallyDistributedIntegrationPoints(numberOfIntegrationPoints);
    getEquallyDistributedIntegrationPoints();

	//addTetrahedronVertexesToIntegrationPoints();

    // get cgal representation of integration points
	CutModelTopology_CGAL::Points_t integrationPoints_CGAL;
	for(Points_t::iterator iter=this->m_globalIntegrationPoints.begin();
			iter!=this->m_globalIntegrationPoints.end();iter++)
	{
		Point_t currentPoint = *iter;
		CutModelTopology_CGAL::Point_t currentPoint_CGAL(currentPoint[0],currentPoint[1],currentPoint[2]);
		integrationPoints_CGAL.push_back(currentPoint_CGAL);
	}

	// check whether integration points are above or below
	if (signEnrichedVertexIds.size())
	{
		this->m_cutModel->pointsAbove(integrationPoints_CGAL, this->m_integrationPointIsAbove);
		for(int integrationPointId=0;integrationPointId<m_numberOfIntegrationPoints;integrationPointId++)
		{
			if(this->m_integrationPointIsAbove[integrationPointId]<-1)
			  std::cout << "WATCH OUT, in Element " << m_tetrahedronId << 
			    " the integration point with Id " << integrationPointId << " lies ON the cut!" << endl;
		}
	}

	// get the coordinates of the integration points in the (r,phi) coordinate system, 
	if (m_branchEnrichedVertexIds.size())
	{
		this->m_integrationPoint_r.resize(this->m_numberOfIntegrationPoints);
		this->m_integrationPoint_phi.resize(this->m_numberOfIntegrationPoints);
		this->m_integrationPoint_x_dependentOnCutFront.resize(this->m_numberOfIntegrationPoints);
		this->m_integrationPoint_y_dependentOnCutFront.resize(this->m_numberOfIntegrationPoints);
		this->m_integrationPointClosestToCutBoundaryElementId.resize(this->m_numberOfIntegrationPoints);

		this->m_cutModel->getR_Phi_XAndYOfPointsDependentOnCutFront(integrationPoints_CGAL,
				this->m_integrationPoint_r,this->m_integrationPoint_phi,
				this->m_integrationPoint_x_dependentOnCutFront,this->m_integrationPoint_y_dependentOnCutFront,
				this->m_integrationPointClosestToCutBoundaryElementId);
	}

	this->m_lambda = lambda;
	this->m_mu = mu;
}


void XFemLinearTetrahedron::getInformationAboutPartialAndCompleteCut()
{
	int edgesOnTetrahedronFaces[4][3] = {{1,5,3},{4,5,2},{0,3,4},{2,1,0}};	// I want this to be a m_ variable - how do I do this?

	for(int edgeId=0;edgeId<6;edgeId++)
	{
		int currentGlobalEdgeId = m_tetrahedronId*6+edgeId;
		int Id;
		m_cutModel->getVectorIdOfIntEntry(currentGlobalEdgeId,m_cutModel->m_objectSplitSegmentIds,Id);
		m_edgeIsSplit[edgeId] = (Id!=-1);
	}
	int helpEdge=0;
	for(int edgeId=0;edgeId<6;edgeId++)
	{
		helpEdge += m_edgeIsSplit[edgeId];
	}
	if(helpEdge)
	{
		int helpElement=0;
		for(int faceId=0;faceId<4;faceId++)
		{
			int helpFace=0;
			for(int edgeIdFace=0;edgeIdFace<3;edgeIdFace++)
			{
				helpFace += m_edgeIsSplit[edgesOnTetrahedronFaces[faceId][edgeIdFace]];//
			}
			m_faceIsSeparated[faceId] = helpFace>1;
			m_faceIsPartiallySeparated[faceId] = helpFace==1;
			helpElement += m_faceIsSeparated[faceId];
		}
		m_tetrahedronIsSeparated = helpElement>2;
		m_tetrahedronIsPartiallySeparated = helpElement<3;
	}
	else
	{
		m_tetrahedronIsSeparated = 0;
		m_tetrahedronIsPartiallySeparated = 0;
		for(int faceId=0;faceId<4;faceId++)
		{
			m_faceIsSeparated[faceId] = 0;
			m_faceIsPartiallySeparated[faceId] = 0;
		}
	}
}


void XFemLinearTetrahedron::getShapeFunctionsOf_m_points()
{
	for(Points_t::iterator iter=this->m_points.begin();iter!=this->m_points.end();iter++)
	{
		Point_t currentPoint = *iter;
		int currentPointId = std::distance(this->m_points.begin(),iter);

		ShapeFct_t currentShapeFct;
		LinearTetrahedron::getShapeFunction(currentPoint,currentShapeFct);

		ShapeFct_t helper;
		LinearTetrahedron::getShapeFunction(currentPoint,helper);
		getShapeFunctionSignEnrichedNode2(this->m_pointIsAbove[currentPointId],helper,currentShapeFct);

		CutModelTopology_CGAL::Point_t currentPointCGAL(currentPoint[0],currentPoint[1],currentPoint[2]);
		if(currentPointId<4 && m_branchEnrichedVertexIds.size())
		{
			this->m_cutModel->getRAndPhiOfPointDependentOnCutFront(currentPointCGAL,m_points_r[currentPointId],m_points_phi[currentPointId]);
		}
		else if (m_branchEnrichedVertexIds.size())
		{
			m_points_r[currentPointId] = this->m_cutModel->getROfPoint(currentPointCGAL);
			m_points_phi[currentPointId] = atan(1.0)*4*(this->m_pointIsAbove[currentPointId]*2-1);	//	atan(1.0)*4 = pi = 3.1415...
		}
		getShapeFunctionBranchEnrichedNode(m_points_r[currentPointId],m_points_phi[currentPointId],helper,currentShapeFct);
		m_shapeFunctions[currentPointId] = currentShapeFct;
	}
}


bool XFemLinearTetrahedron::isEnriched()
{
	return m_numberOfElementNodes>4;
}


void XFemLinearTetrahedron::getShapeFunctionSignEnrichedNode(Point_t x, ShapeFct_t &ShapeFunction)
{
	ShapeFct_t helpShapeFunction;
	LinearTetrahedron::getShapeFunction(x,helpShapeFunction);
	CutModelTopology_CGAL::Point_t x_CGAL(x[0],x[1],x[2]);
	int pointIsAbove = m_cutModel->pointAbove(x_CGAL);
	if(pointIsAbove == -1)
	  std::cout << "WATCH OUT, in Element " << m_tetrahedronId << " a point lies ON the cut!" << endl;

	getShapeFunctionSignEnrichedNode2(pointIsAbove, helpShapeFunction, ShapeFunction);
}


void XFemLinearTetrahedron::getShapeFunctionSignEnrichedNode2(int pointIsAbove, ShapeFct_t ShapeFunctionIn, ShapeFct_t &ShapeFunctionOut)
{
	for(Ints_t::iterator iter=m_signEnrichedVertexIds.begin();iter!=m_signEnrichedVertexIds.end();iter++)
	{
		int currentNodeId = *iter;
		double currentShapeFctValue = (double)(pointIsAbove-m_pointIsAbove[currentNodeId])*ShapeFunctionIn[currentNodeId];
		ShapeFunctionOut.push_back(currentShapeFctValue);
	}
}


void XFemLinearTetrahedron::getShapeFunctionBranchEnrichedNode(double r, double phi, ShapeFct_t ShapeFunctionIn, ShapeFct_t &ShapeFunctionOut)
{
	Doubles_t currentACTFs;
	currentACTFs.resize(4);
	getAsymptoticCrackTipFunctions(r,phi,currentACTFs);

	for(Ints_t::iterator iter=m_branchEnrichedVertexIds.begin();iter!=m_branchEnrichedVertexIds.end();iter++)
	{
		int currentNodeId = *iter;
		for(int ACTFnumber=0;ACTFnumber<4;ACTFnumber++)
		{
			double currentShapeFctValue = currentACTFs[ACTFnumber]*ShapeFunctionIn[currentNodeId];
			ShapeFunctionOut.push_back(currentShapeFctValue);
		}
	}
}


void XFemLinearTetrahedron::getAsymptoticCrackTipFunctions(double r, double phi, Doubles_t &ACTFs)
{
	ACTFs.resize(4);

	ACTFs[0] = sqrt(r)*sin(phi/2);
	ACTFs[1] = sqrt(r)*cos(phi/2);
	ACTFs[2] = sqrt(r)*sin(phi)*sin(phi/2);
	ACTFs[3] = sqrt(r)*sin(phi)*cos(phi/2);
}


void XFemLinearTetrahedron::getShapeFunction(Point_t x, int pointIsAbove, double r, double phi, ShapeFct_t &ShapeFunction)
{
	ShapeFct_t helper;
	LinearTetrahedron::getShapeFunction(x,helper);

	// for the unenriched nodes
	LinearTetrahedron::getShapeFunction(x,ShapeFunction);



	// for the sign enriched nodes
	getShapeFunctionSignEnrichedNode2(pointIsAbove, helper, ShapeFunction);

	// for the branch enriched nodes
	getShapeFunctionBranchEnrichedNode(r, phi, helper, ShapeFunction);
}


void XFemLinearTetrahedron::derivativeACTFsWithRespectToRAndPhi(double r,double phi,std::vector<Doubles_t> &dACTFs_drdphi)
{
    // initialization of the derivative
    dACTFs_drdphi.resize(4);
    for(int ACTFnumber=0;ACTFnumber<4;ACTFnumber++)
        dACTFs_drdphi[ACTFnumber].resize(2);

    // fill the entries of the derivative
    Doubles_t currentDACTF;
    currentDACTF.resize(2);

    currentDACTF[0] = sin(phi/2)/(2*sqrt(r));
    currentDACTF[1] = 0.5*sqrt(r)*cos(phi/2);
    dACTFs_drdphi[0] = currentDACTF;

    currentDACTF[0] = cos(phi/2)/(2*sqrt(r));
    currentDACTF[1] = -0.5*sqrt(r)*sin(phi/2);
    dACTFs_drdphi[1] = currentDACTF;

    currentDACTF[0] = sin(phi/2)*sin(phi)/(2*sqrt(r));
    currentDACTF[1] = sqrt(r)*(0.5*cos(phi/2)*sin(phi)+sin(phi/2)*cos(phi));
    dACTFs_drdphi[2] = currentDACTF;

    currentDACTF[0] = cos(phi/2)*sin(phi)/(2*sqrt(r));
    currentDACTF[1] = sqrt(r)*(-0.5*sin(phi/2)*sin(phi)+cos(phi/2)*cos(phi));
    dACTFs_drdphi[3] = currentDACTF;
}


void XFemLinearTetrahedron::getDerivativeACTFs(int integrationPointNumber, ShapeFctDeriv_t &dACTFs)
{
	std::vector<Doubles_t> dACTFs_drdphi;
    double r = this->m_integrationPoint_r[integrationPointNumber];
    /*// avoid dividing by zero
    if(abs(r) < 1e-300)
        r = 1e-20 * ( (r>=0) - (r<0) );*/
	double phi = this->m_integrationPoint_phi[integrationPointNumber];
	double x = this->m_integrationPoint_x_dependentOnCutFront[integrationPointNumber];
    double y = this->m_integrationPoint_y_dependentOnCutFront[integrationPointNumber];

	derivativeACTFsWithRespectToRAndPhi(r,phi,dACTFs_drdphi);
	
	Matrix_3_3_t currentCutBoundaryTransformationMatrix = this->m_cutModel->m_cutBoundaryTransformationMatrices[this->m_integrationPointClosestToCutBoundaryElementId[integrationPointNumber]];
	Matrix_3_3_t cICBTM;	// cICBTM = currentInverseCutBoundaryTransformationMatrix

	cICBTM = currentCutBoundaryTransformationMatrix.inverse();

	dACTFs.resize(4);
	for(int ACTFNumber=0;ACTFNumber<4;ACTFNumber++)
	{
        Point_t currentDACTFs;
        //std::cout << std::numeric_limits<double>::epsilon() << std::endl;
        if ( r > std::numeric_limits<double>::epsilon() )   // ensures that we do not devide by 0 !
        {
            double dACTFs_dr = (dACTFs_drdphi[ACTFNumber])[0];
            double dACTFs_dphi = (dACTFs_drdphi[ACTFNumber])[1];
            for(int j=0;j<DIM;j++)
            {
                double dr_dxj = (cICBTM(0,j)*x+cICBTM(1,j)*y)/r;
                double dphi_dxj = (-cICBTM(0,j)*y+cICBTM(1,j)*x)/(r*r);

                currentDACTFs[j] = dACTFs_dr*dr_dxj+dACTFs_dphi*dphi_dxj;
            }
        }
		dACTFs[ACTFNumber] = currentDACTFs;
	}
}


void XFemLinearTetrahedron::getShapeFunctionDerivativesDependentOnIntegrationPointNumber(int integrationPointNumber, ShapeFctDeriv_t &ShapeFunctionDerivative)
{
	// get the standard derivative
	LinearTetrahedron::getShapeFunctionDerivativesDependentOnIntegrationPointNumber(integrationPointNumber,ShapeFunctionDerivative);

	// add the derivatives at the sign enriched element nodes
	Point_t derivativeOfCurrentShapeFunction;
	for(Ints_t::iterator iter=this->m_signEnrichedVertexIds.begin();
			iter!=this->m_signEnrichedVertexIds.end(); iter++)
	{
		int currentElementNodeId = *iter;
		derivativeOfCurrentShapeFunction = ShapeFunctionDerivative[currentElementNodeId];
		for(int i=0;i<DIM;i++)
		{
			derivativeOfCurrentShapeFunction[i] *= (double)(this->m_integrationPointIsAbove[integrationPointNumber]-this->m_pointIsAbove[currentElementNodeId]);
		}
		ShapeFunctionDerivative.push_back(derivativeOfCurrentShapeFunction);
	}

	// add the derivatives at the sign enriched element nodes
	ShapeFct_t ACTFs;
	ShapeFctDeriv_t dACTFs;
	ShapeFct_t currentShapeFunction;
	if(this->m_branchEnrichedVertexIds.begin()!=this->m_branchEnrichedVertexIds.end())
	{
		getAsymptoticCrackTipFunctions(this->m_integrationPoint_r[integrationPointNumber],
            this->m_integrationPoint_phi[integrationPointNumber],ACTFs);
        getDerivativeACTFs(integrationPointNumber,dACTFs);
		LinearTetrahedron::getShapeFunction(this->m_globalIntegrationPoints[integrationPointNumber], currentShapeFunction);
	}
	for(Ints_t::iterator iter=this->m_branchEnrichedVertexIds.begin();
			iter!=this->m_branchEnrichedVertexIds.end(); iter++)
	{
		int currentElementNodeId = *iter;
		for(int ACTFNumber=0;ACTFNumber<4;ACTFNumber++)
		{
			Point_t helper = ShapeFunctionDerivative[currentElementNodeId];
			
			Point_t sd(	ACTFs[ACTFNumber]*helper[0] + currentShapeFunction[currentElementNodeId]*((dACTFs[ACTFNumber])[0]),
                        ACTFs[ACTFNumber]*helper[1] + currentShapeFunction[currentElementNodeId]*((dACTFs[ACTFNumber])[1]),
                        ACTFs[ACTFNumber]*helper[2] + currentShapeFunction[currentElementNodeId]*((dACTFs[ACTFNumber])[2])  );
			ShapeFunctionDerivative.push_back(sd);
        }
	}

    //WARNING("currentDirectory is hard coded");
    std::string currentDirectory = "/home/cpaulus/dev/xfem";
    std::string filenameCurrentDerivative = currentDirectory + "/output/ShapeFunctionDerivatives/" + boost::lexical_cast<std::string>(this->m_tetrahedronId) + "/" + boost::lexical_cast<std::string>(integrationPointNumber);
    writeShapeFunctionDerivativeInFile(ShapeFunctionDerivative,filenameCurrentDerivative);
}


void XFemLinearTetrahedron::getCellIdsForOutput(int above, Ints_t &cellIds)
{
	if(above == -1)
		LinearTetrahedron::getCellIdsForOutput(cellIds);
	else
		m_cutModel->getVectorIdsOfIntEntry(above,m_triangleAbove,cellIds);
}
