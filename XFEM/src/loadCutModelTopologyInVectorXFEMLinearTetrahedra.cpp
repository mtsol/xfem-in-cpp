/* loadCutModelTopologyInVectorXFEMLinearTetrahedra.cpp
 *
 * 		Created On: 19.07.2013
 * 		Author: 	Christoph Joachim Paulus
 *
 * */

// includes
#include "loadCutModelTopologyInVectorXFEMLinearTetrahedra.h"


// function definitions
loadCutModelTopologyInVectorXFEMLinearTetrahedra::loadCutModelTopologyInVectorXFEMLinearTetrahedra()
{
}


loadCutModelTopologyInVectorXFEMLinearTetrahedra::~loadCutModelTopologyInVectorXFEMLinearTetrahedra()
{
}


void loadCutModelTopologyInVectorXFEMLinearTetrahedra::init(std::string ObjectFilename, std::string CutFilename, std::string RefGridFilename, 
		double lambda,double mu,int numberOfIntegrationPoints, Ints_t &pointsReferenceSolutionInObjectTetrahedronId, Ints_t &pointsReferenceSolutionAbove, 
		Doubles_t &pointsReferenceSolution_r, Doubles_t &pointsReferenceSolution_phi, 
        std::vector<XFemLinearTetrahedron> &VectorXFemLinearTetrahedra) // make the last few inputs optional! Maybe use struct for output...!
{
//    DEBUG;
	// initialization of the cut model topology
    m_cutModel.init(ObjectFilename,CutFilename);

	// UNCOMMENT THE FOLLOWING FOR THE REFERENCE SOLUTION
//    m_cutModel.pointsInObjectTetrahedronNumber(RefGridFilename, pointsReferenceSolutionInObjectTetrahedronId, pointsReferenceSolutionAbove, pointsReferenceSolution_r, pointsReferenceSolution_phi);

	int numberOfObjectTetrahedra = std::distance(m_cutModel.m_objectConnectTetrahedra.begin(),m_cutModel.m_objectConnectTetrahedra.end());
	int numberOfObjectPoints = std::distance(m_cutModel.m_intersectionPoints.begin(),m_cutModel.m_intersectionPoints.end());

	// vtkPoints* pointContainerVTK = vtkPoints::New();
	// int numberOfAlreadyExistingPoints = 0;
	// vtkCellArray* cellContainerVTK = vtkCellArray::New();

	// fill the vector of xfem linear tetrahedra with the information of the cut model topology
	for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
	{
		XFemLinearTetrahedron::Points_t pointsInCurrentTetrahedron;
		Ints_t pointIds;
		Ints_t pointIsAbove;
		Ints_t signEnrichedVertexIds;
		Ints_t branchEnrichedVertexIds;

        informationAboutPointsOnTetrahedron(objectTetrahedronId,pointIsAbove,signEnrichedVertexIds,branchEnrichedVertexIds,pointIds,pointsInCurrentTetrahedron);

		XFemLinearTetrahedron::ConnectTriangles_t TetrahedronConnectTriangle;
		Ints_t triangleInFaceId;
		Ints_t triangleAbove;
		getConnectTriangle(objectTetrahedronId,pointIds,pointIsAbove,TetrahedronConnectTriangle,triangleInFaceId,triangleAbove);

		Ints_t connect = m_cutModel.m_objectConnectTetrahedraWithEnrichedNodes[objectTetrahedronId];

		int neighboursOfFaces[4][2];
		for(int faceId=0;faceId<4;faceId++)
		{
			for(int i=0;i<2;i++)
				neighboursOfFaces[faceId][i]=-1;
			int currentNeighbourGlobalFaceId = m_cutModel.m_globalFaceIdHasNeighbourWithglobalFaceId[objectTetrahedronId*4+faceId];
			if(currentNeighbourGlobalFaceId!=-1)
			{
				neighboursOfFaces[faceId][0] = currentNeighbourGlobalFaceId/4;
				neighboursOfFaces[faceId][1] = currentNeighbourGlobalFaceId%4;
			}
		}

		XFemLinearTetrahedron currentXFemLinearTetrahedron;

        currentXFemLinearTetrahedron.init(objectTetrahedronId,connect,neighboursOfFaces,pointsInCurrentTetrahedron,pointIsAbove,signEnrichedVertexIds,branchEnrichedVertexIds,&m_cutModel,TetrahedronConnectTriangle,triangleInFaceId,triangleAbove,lambda,mu,numberOfIntegrationPoints);

		VectorXFemLinearTetrahedra.push_back(currentXFemLinearTetrahedron);
	}
}


void loadCutModelTopologyInVectorXFEMLinearTetrahedra::informationAboutPointsOnTetrahedron(int tetrahedronId,
		Ints_t &pointIsAbove, Ints_t &signEnrichedVertexIds, Ints_t &branchEnrichedVertexIds, Ints_t &pointIds, XFemLinearTetrahedron::Points_t &points)
{
	for(int localPointId=0;localPointId<4;localPointId++)
	{
		int currentPointId = (m_cutModel.m_objectConnectTetrahedraWithEnrichedNodes[tetrahedronId])[localPointId];
		pointIds.push_back(currentPointId);
		int Id;
		m_cutModel.getVectorIdOfIntEntry(currentPointId,m_cutModel.m_signEnrichedObjectPoints,Id);
		if(Id!=-1)
			signEnrichedVertexIds.push_back(localPointId);

		m_cutModel.getVectorIdOfIntEntry(currentPointId,m_cutModel.m_branchEnrichedObjectPoints,Id);
		if(Id!=-1)
			branchEnrichedVertexIds.push_back(localPointId);
	}
	m_cutModel.getVectorIdsOfIntEntry(tetrahedronId,m_cutModel.m_intersectionPointInObjectTetrahedronId,pointIds);
	for(Ints_t::iterator iter=pointIds.begin();iter!=pointIds.end();iter++)
	{
		int currentId = *iter;
		CutModelTopology_CGAL::Point_t currentPointCGAL = m_cutModel.m_intersectionPoints[currentId];
		XFemLinearTetrahedron::Point_t currentPoint;
		for(int i=0;i<3;i++)
			currentPoint[i] = currentPointCGAL[i];
		points.push_back(currentPoint);
		pointIsAbove.push_back(m_cutModel.m_intersectionPointAbove[currentId]);
	}
}


void loadCutModelTopologyInVectorXFEMLinearTetrahedra::getConnectTriangle(int tetrahedronId,Ints_t &pointIds,Ints_t &pointIsAbove,
		XFemLinearTetrahedron::ConnectTriangles_t &connectTriangle, Ints_t &triangleInFaceId, Ints_t &triangleAbove)
{
	int numberOfPointsInCurrentTetrahedron = std::distance(pointIds.begin(),pointIds.end());
	int nodesOnTetrahedronFaces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}}; // I want this to be a m_ variable - how do I do this?


	XFemLinearTetrahedron::ConnectTriangle_t currentConnectTriangle;
	for(int faceId=0;faceId<4;faceId++)
	{
		if(!faceIntersectsWithCut(tetrahedronId*4+faceId))
		{
			for(int pointId=0;pointId<3;pointId++)
			{
				currentConnectTriangle[pointId] = nodesOnTetrahedronFaces[faceId][pointId];
			}
			connectTriangle.push_back(currentConnectTriangle);
			triangleInFaceId.push_back(faceId);
			triangleAbove.push_back(pointIsAbove[nodesOnTetrahedronFaces[faceId][0]]);//vertexIsAbove[nodesOnTetrahedronFaces[faceId][0]]
		}
	}

	Ints_t triangleIds;
	m_cutModel.getVectorIdsOfIntEntry(tetrahedronId,m_cutModel.m_TriangleOfTriangulationInObjectTetrahedronId,triangleIds);
	for(Ints_t::iterator iter=triangleIds.begin();iter!=triangleIds.end();iter++)
	{
		int triangleId = *iter;
		int IdIn_m_cutModel_m_intersectionPoints = m_cutModel.m_triangulationOfCutObjectTetrahedronFacesConnectTriangle[triangleId].get<0>();
		int IdInPointsInCurrentTetrahedron;
		m_cutModel.getVectorIdOfIntEntry(IdIn_m_cutModel_m_intersectionPoints,pointIds,IdInPointsInCurrentTetrahedron);
		currentConnectTriangle[0] = IdInPointsInCurrentTetrahedron;
		IdIn_m_cutModel_m_intersectionPoints = m_cutModel.m_triangulationOfCutObjectTetrahedronFacesConnectTriangle[triangleId].get<1>();
		m_cutModel.getVectorIdOfIntEntry(IdIn_m_cutModel_m_intersectionPoints,pointIds,IdInPointsInCurrentTetrahedron);
		currentConnectTriangle[1] = IdInPointsInCurrentTetrahedron;
		IdIn_m_cutModel_m_intersectionPoints = m_cutModel.m_triangulationOfCutObjectTetrahedronFacesConnectTriangle[triangleId].get<2>();
		m_cutModel.getVectorIdOfIntEntry(IdIn_m_cutModel_m_intersectionPoints,pointIds,IdInPointsInCurrentTetrahedron);
		currentConnectTriangle[2] = IdInPointsInCurrentTetrahedron;

		connectTriangle.push_back(currentConnectTriangle);
		triangleInFaceId.push_back(m_cutModel.m_TriangleOfTriangulationInObjectTriangleId[triangleId]%4);
		triangleAbove.push_back(m_cutModel.m_TriangleOfTriangulationAbove[triangleId]);
	}
}


bool loadCutModelTopologyInVectorXFEMLinearTetrahedra::faceIntersectsWithCut(int globalFaceId)
{
	int edgesOnTetrahedronFaces[4][3] = {{1,5,3},{4,5,2},{0,3,4},{2,1,0}};	// I want this to be a m_ variable - how do I do this?

	for(int edgeId=0;edgeId<3;edgeId++)
	{
		int currentGlobalEdgeId = (globalFaceId/4)*6+edgesOnTetrahedronFaces[globalFaceId%4][edgeId];
		int Id;
		m_cutModel.getVectorIdOfIntEntry(currentGlobalEdgeId,m_cutModel.m_objectSplitSegmentIds,Id);
		if(Id!=-1)
		{
			return true;
		}
	}
	return false;
}


int loadCutModelTopologyInVectorXFEMLinearTetrahedra::getNumberOfSignEnrichedPoints()
{
	int numberOfSignEnrichedVertexes = std::distance(m_cutModel.m_signEnrichedObjectPoints.begin(),m_cutModel.m_signEnrichedObjectPoints.end());

	return numberOfSignEnrichedVertexes;
}


int loadCutModelTopologyInVectorXFEMLinearTetrahedra::getNumberOfBranchEnrichedPoints()
{
	int numberOfBranchEnrichedVertexes = std::distance(m_cutModel.m_branchEnrichedObjectPoints.begin(),m_cutModel.m_branchEnrichedObjectPoints.end());

	return 4*numberOfBranchEnrichedVertexes;
}


int loadCutModelTopologyInVectorXFEMLinearTetrahedra::getNumberOfPoints()
{
	int numberOfObjectPoints = std::distance(m_cutModel.m_objectPoints.begin(),m_cutModel.m_objectPoints.end());
    int numberOfSignEnrichedVertexes = std::distance(m_cutModel.m_signEnrichedObjectPoints.begin(),m_cutModel.m_signEnrichedObjectPoints.end());
    int numberOfBranchEnrichedVertexes = std::distance(m_cutModel.m_branchEnrichedObjectPoints.begin(),m_cutModel.m_branchEnrichedObjectPoints.end());

    return (numberOfObjectPoints+numberOfSignEnrichedVertexes+4*numberOfBranchEnrichedVertexes);
}


//void loadCutModelTopologyInVectorXFEMLinearTetrahedra::outputFromVTKContainersToVTKFile(vtkPoints* pointContainer, vtkCellArray* cellContainer, int numberOfCellNodes, const char* filename)
//{
//	vtkUnstructuredGrid* myGridTopoDebug = vtkUnstructuredGrid::New();
//	myGridTopoDebug->SetPoints(pointContainer);
//	if(numberOfCellNodes==3)
//	{
//		myGridTopoDebug->SetCells(VTK_TRIANGLE, cellContainer);
//	}
//	else if (numberOfCellNodes==2)
//	{
//		myGridTopoDebug->SetCells(VTK_LINE, cellContainer);
//	}
//
//	vtkUnstructuredGridWriter* writerTopoDebug = vtkUnstructuredGridWriter::New();
//
//	writerTopoDebug->SetInput(myGridTopoDebug);
//	writerTopoDebug->SetFileName(filename);
//	writerTopoDebug->Write();
//
//	pointContainer->Delete();
//	cellContainer->Delete();
//	myGridTopoDebug->Delete();
//	writerTopoDebug->Delete();
//}
