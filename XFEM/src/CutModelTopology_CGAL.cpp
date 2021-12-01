/* CutModelTopology_CGAL.cpp
 *
 * 		Created On: 14.07.2013
 * 		Author: 	Christoph Joachim Paulus
 *
 * */

// includes
#include "CutModelTopology_CGAL.h"


// function definitions
CutModelTopology_CGAL::CutModelTopology_CGAL(void)
{
}


CutModelTopology_CGAL::~CutModelTopology_CGAL()
{
}


void CutModelTopology_CGAL::init(std::string ObjectFilename, std::string CutFilename)
{
	/*INTERNAL DATA*/
	// define properties of the tetrahedra and triangles
	int numberOfNodesSegment = 2;
	int numberOfNodesTriangle = 3;
	int numberOfNodesTetrahedron = 4;

	int nodesOnTriangleFaces[1][3] = {0,1,2};
	int nodesOnTriangleEdges[3][2] = {{0,1},{1,2},{2,0}};	// this helps that the normal always shows in the same direction, no matter whether we turn the triangle in a direction without interchanging the vertexes
	int nodeNotOnTriangleEdges[3] = {2,0,1};	// the node number that is not in edge number...
	int edgesOnTriangleNodes[3][2] = {{0,2},{0,1},{1,2}};

	int nodesOnTetrahedronFaces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
	int nodesOnTetrahedronEdges[6][2] = {{0,1},{1,2},{0,2},{1,3},{0,3},{2,3}};
	int edgesOnTetrahedronFaces[4][3] = {{1,5,3},{4,5,2},{0,3,4},{2,1,0}};
//	int nodesOnTetrahedronFaces[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
//	int nodesOnTetrahedronEdges[6][2] = {{0,1},{1,2},{0,2},{1,3},{0,3},{2,3}};
//	int edgesOnTetrahedronFaces[4][3] = {{1,5,3},{2,5,4},{0,3,4},{0,1,2}};





	/*GET OBJECT*/

    // define the reader, the grid, the points and the array, to receive information read from the vtk file
    vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
    reader->SetFileName(ObjectFilename.c_str());
    reader->Update();
    vtkUnstructuredGrid* myGrid = reader->GetOutput();
    vtkPoints* myPoints = myGrid->GetPoints();
    vtkCellArray* myArray = myGrid->GetCells();

    // get the connect of the tetrahedra grid
    vtkIdType numberOfCells = myGrid->GetNumberOfCells();
    int objectConnectTetrahedra[int(numberOfCells)][4];
    int numberOfObjectPoints;
    getTetConnectFromVTKCellArray(myArray, numberOfNodesTetrahedron, objectConnectTetrahedra, numberOfCells, numberOfObjectPoints);
    int numberOfObjectTetrahedra = numberOfCells;

    // receive the points of the tetrahedra grid
    double objectPointsDouble[numberOfObjectPoints][3];
    getPointsFromVTKCellArray(numberOfObjectPoints, myPoints, 1, objectPointsDouble, this->m_objectPoints);

    // delete the reader that it can be used for the triangle grid (the cut)
    reader->Delete();

    // write the connect of the object in two different ways - as triangle connect and segment connect
    int objectConnectTriangle[numberOfObjectTetrahedra*4][3];
    int objectConnectSegment[numberOfObjectTetrahedra*6][2];
    getConnectTriangleOfConnectTet(numberOfObjectTetrahedra, objectConnectTetrahedra, nodesOnTetrahedronFaces, objectConnectTriangle);
    getConnectSegmentOfConnectTet(numberOfObjectTetrahedra, objectConnectTetrahedra, nodesOnTetrahedronEdges, objectConnectSegment);
//	m_nodeIdInTetrahedraIds.resize(numberOfObjectPoints);
    for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
    {
        ConnectTetrahedron_t currentConnectTetrahedron(objectConnectTetrahedra[objectTetrahedronId][0],
                objectConnectTetrahedra[objectTetrahedronId][1],
                objectConnectTetrahedra[objectTetrahedronId][2],
                objectConnectTetrahedra[objectTetrahedronId][3]);
        this->m_objectConnectTetrahedra.push_back(currentConnectTetrahedron);
//		for(int nodeId=0;nodeId<numberOfNodesTetrahedron;nodeId++)
//		{
//			m_nodeIdInTetrahedraIds[objectConnectTetrahedra[objectTetrahedronId][nodeId]].push_back(objectTetrahedronId);
//		}
    }

    // check which faces have neighbours
    std::vector<Ints_t> allSortedFaces;
    allSortedFaces.resize(numberOfObjectTetrahedra*4);
    m_globalFaceIdHasNeighbourWithglobalFaceId.resize(numberOfObjectTetrahedra*4);
    for(int globalFaceId=0;globalFaceId<numberOfObjectTetrahedra*4;globalFaceId++)
    {
        Ints_t currentFace;
        for(int nodeId=0;nodeId<numberOfNodesTriangle;nodeId++)
            currentFace.push_back(objectConnectTriangle[globalFaceId][nodeId]);
        std::sort(currentFace.begin(),currentFace.end());
        allSortedFaces[globalFaceId] = currentFace;
        m_globalFaceIdHasNeighbourWithglobalFaceId[globalFaceId] = -1;
    }
    for(int globalFaceId=0;globalFaceId<numberOfObjectTetrahedra*4;globalFaceId++)
    {
        bool neighbourNotAssigned = m_globalFaceIdHasNeighbourWithglobalFaceId[globalFaceId] == -1;
        int currentGlobalFaceId = globalFaceId+1;
        while(neighbourNotAssigned && currentGlobalFaceId<(numberOfObjectTetrahedra*4))
        {
            if(	(allSortedFaces[globalFaceId])[0]==(allSortedFaces[currentGlobalFaceId])[0]	&&
                (allSortedFaces[globalFaceId])[1]==(allSortedFaces[currentGlobalFaceId])[1]	&&
                (allSortedFaces[globalFaceId])[2]==(allSortedFaces[currentGlobalFaceId])[2])
            {
                m_globalFaceIdHasNeighbourWithglobalFaceId[globalFaceId] = currentGlobalFaceId;
                m_globalFaceIdHasNeighbourWithglobalFaceId[currentGlobalFaceId] = globalFaceId;
                neighbourNotAssigned = false;
            }
            currentGlobalFaceId++;
        }
    }


    // write the object as triangles and segments instead of tetrahedra
    std::vector<Triangle_t> objectTriangles;
    std::vector<Segment_t> objectSegments;
    getTrianglesFromConnectTriangleAndPointsDouble(numberOfObjectTetrahedra*4, objectConnectTriangle, objectPointsDouble, objectTriangles);
    getSegmentsFromConnectSegmentAndPointsDouble(numberOfObjectTetrahedra*6, objectConnectSegment, objectPointsDouble, objectSegments);
    this->m_objectTriangles = objectTriangles;

    // constructs AABB tree
    Tree_t objectAABBTree(objectTriangles.begin(),objectTriangles.end());





	/*GET CUT*/

	// define the reader, the grid, the points and the array, to receive information read from the vtk file
	reader = vtkUnstructuredGridReader::New();
	reader->SetFileName(CutFilename.c_str());
	reader->Update();
	myGrid = reader->GetOutput();
	myPoints = myGrid->GetPoints();
	vtkIdType numberOfCellNodes_vtkIdType = numberOfNodesTriangle;
	myArray = myGrid->GetCells();

	// get the connect of the cut in triangles
	numberOfCells = myGrid->GetNumberOfCells();
	int CutConnectTriangle[int(numberOfCells)][3];
	int	numberOfCutPoints;
	getTriangleConnectFromVTKCellArray(myArray, numberOfNodesTriangle, CutConnectTriangle, numberOfCells, numberOfCutPoints);
	int numberOfCutTriangles = int(numberOfCells);

	// get the points of the cut
	Points_t cutPoints;
	double cutGridPoints[numberOfCutPoints][3];
	getPointsFromVTKCellArray(numberOfCutPoints, myPoints, 1, cutGridPoints, cutPoints);

	// clean up
	reader->Delete();

	// write the triangle connect of the cut as segment connect
	int CutConnectSegment[numberOfCutTriangles*3][2];
	getConnectSegmentOfConnectTriangle(numberOfCutTriangles, CutConnectTriangle, CutConnectSegment);

	// write the cut as a vector of Triangle
	std::vector<Triangle_t> CutTriangles;
	std::vector<Segment_t> CutSegments;
	getTrianglesFromConnectTriangleAndPointsDouble(numberOfCutTriangles, CutConnectTriangle, cutGridPoints, CutTriangles);
	getSegmentsFromConnectSegmentAndPointsDouble(numberOfCutTriangles*3, CutConnectSegment, cutGridPoints, CutSegments);
	this->m_cutTriangles = CutTriangles;

	// get the neighbours of the cut cells
	int CutNeighbours[numberOfCutTriangles][3];
	for(Segments_t::iterator CutSegmentIterator1=CutSegments.begin();CutSegmentIterator1!=CutSegments.end();CutSegmentIterator1++)
	{
		Segment_t CutSegment1 = *CutSegmentIterator1;
		int CutSegment1Id = std::distance(CutSegments.begin(),CutSegmentIterator1);
		int CutTriangle1Id = CutSegment1Id/3;
		CutNeighbours[CutTriangle1Id][CutSegment1Id%3]=-1;
		for(Segments_t::iterator CutSegmentIterator2=CutSegments.begin();CutSegmentIterator2!=CutSegments.end();CutSegmentIterator2++)
		{
			Segment_t CutSegment2 = *CutSegmentIterator2;
			int CutSegment2Id = std::distance(CutSegments.begin(),CutSegmentIterator2);
			int CutTriangle2Id = CutSegment2Id/3;
			if((CutSegment1==CutSegment2 || CutSegment1.opposite()==CutSegment2) && CutSegment1Id!=CutSegment2Id)
			{
				CutNeighbours[CutTriangle1Id][CutSegment1Id%3] = CutTriangle2Id;
			}
		}
	}

	// construct AABB tree for the cut
	Tree_t CutAABBTree(CutTriangles.begin(),CutTriangles.end());





	/*CHECK WHICH NODES ARE ABOVE*/
	int pointsAbove[numberOfObjectPoints];
	isPointAbove(this->m_objectPoints, CutAABBTree, CutTriangles, pointsAbove, this->m_pointIdsAbove);
	int numberOfPointsAbove=std::distance(this->m_pointIdsAbove.begin(),this->m_pointIdsAbove.end());
	for(int nodeId=0;nodeId<numberOfObjectPoints;nodeId++)
		if(pointsAbove[nodeId]==-1)
			std::cout << "WATCH OUT the node " << nodeId << " lies on cut - code can't handle that!" << endl;



	/*CHECK WHICH CELLS ARE CUT AND HOW THEY ARE CUT*/
	int objectSegmentIsSplit[numberOfObjectTetrahedra][6];
	int objectTriangleIsSeparated[numberOfObjectTetrahedra][4];
	int objectTrianglePartiallySeparated[numberOfObjectTetrahedra][4];
	int objectTetrahedronIsSeparated[numberOfObjectTetrahedra];
	int objectTetrahedronPartiallySeparated[numberOfObjectTetrahedra];

	getInformationAboutPartialAndCompleteCut(objectSegments, CutAABBTree, edgesOnTetrahedronFaces,
			objectSegmentIsSplit, this->m_objectSplitSegmentIds,
			objectTriangleIsSeparated, objectTrianglePartiallySeparated,this->m_objectTriangleIdsThatIntersectWithCut,
			objectTetrahedronIsSeparated, objectTetrahedronPartiallySeparated, this->m_objectTetrahedraIdsThatIntersectWithCut);





	/*CHECK HOW THE OBJECT POINTS ARE ENRICHED*/
	bool objectPointIsSignEnriched[numberOfObjectPoints];
	bool objectPointIsBranchEnriched[numberOfObjectPoints];
	for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
	{
		for(int pointId=0;pointId<4;pointId++)
		{
			objectPointIsSignEnriched[objectConnectTetrahedra[objectTetrahedronId][pointId]] = false;
			objectPointIsBranchEnriched[objectConnectTetrahedra[objectTetrahedronId][pointId]] = false;
		}
	}
	for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
	{
		if(objectTetrahedronIsSeparated[objectTetrahedronId])
		{
			for(int pointId=0;pointId<4;pointId++)
			{
				objectPointIsSignEnriched[objectConnectTetrahedra[objectTetrahedronId][pointId]] = true;
			}
		}
		else if(objectTetrahedronPartiallySeparated[objectTetrahedronId])
		{
			for(int pointId=0;pointId<4;pointId++)
			{
				objectPointIsBranchEnriched[objectConnectTetrahedra[objectTetrahedronId][pointId]] = true;
			}
		}
	}
	for(int objectPointId=0;objectPointId<numberOfObjectPoints;objectPointId++)
	{
		if(objectPointIsBranchEnriched[objectPointId])
		{
			objectPointIsSignEnriched[objectPointId] = false;
			this->m_branchEnrichedObjectPoints.push_back(objectPointId);
		}
	}
	//std::cout << endl;
	for(int objectPointId=0;objectPointId<numberOfObjectPoints;objectPointId++)
	{
		if(objectPointIsSignEnriched[objectPointId])
		{
			this->m_signEnrichedObjectPoints.push_back(objectPointId);
		}
	}
	m_numberOfSignEnrichedNodes = std::distance(m_signEnrichedObjectPoints.begin(),m_signEnrichedObjectPoints.end());
	m_numberOfBranchEnrichedNodes = 4*std::distance(m_branchEnrichedObjectPoints.begin(),m_branchEnrichedObjectPoints.end());





	/* GET THE CONNECT WITH THE ENRICHED NODES*/
	for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
	{
		Ints_t currentConnect;
		for(int nodeId=0;nodeId<numberOfNodesTetrahedron;nodeId++)
		{
			currentConnect.push_back(objectConnectTetrahedra[objectTetrahedronId][nodeId]);
		}
		if(objectTetrahedronIsSeparated[objectTetrahedronId])
		{
			for(int nodeId=0;nodeId<numberOfNodesTetrahedron;nodeId++)
			{
				int objectPointId = objectConnectTetrahedra[objectTetrahedronId][nodeId];
				if(objectPointIsSignEnriched[objectPointId])
				{
					int Id;
					getVectorIdOfIntEntry(objectPointId,m_signEnrichedObjectPoints,Id);
					currentConnect.push_back(Id+numberOfObjectPoints);
				}
			}
		}
		for(int nodeId=0;nodeId<numberOfNodesTetrahedron;nodeId++)
		{
			int objectPointId = objectConnectTetrahedra[objectTetrahedronId][nodeId];
			if(objectPointIsBranchEnriched[objectPointId])
			{
				int Id;
				getVectorIdOfIntEntry(objectPointId,m_branchEnrichedObjectPoints,Id);
				for(int asymptoticCrackTipFunctionNumber=0;asymptoticCrackTipFunctionNumber<4;asymptoticCrackTipFunctionNumber++)
				{
					currentConnect.push_back(asymptoticCrackTipFunctionNumber+4*Id+numberOfObjectPoints+m_numberOfSignEnrichedNodes);
				}
			}
		}
		m_objectConnectTetrahedraWithEnrichedNodes.push_back(currentConnect);
	}





	/*GET THE INTERSECTIONS OF THE CUT SEGMENTS WITH THE OBJECT TRIANGLES*/

	Ints_t CutSegmentIdsThatIntersectWithObject;
	Ints_t objectTetrahedronIdOfCutSegmentIdsThatIntersectWithObject;

	// put the vertexes of the object tetrahedra into the intersection points, such that we can use them for the triangulation later
	for(int objectPointId=0;objectPointId<numberOfObjectPoints;objectPointId++)
	{
		this->m_intersectionPoints.push_back(this->m_objectPoints[objectPointId]);

		this->m_intersectionPointInObjectTetrahedronId.push_back(-1);
		this->m_intersectionPointInObjectTriangleId.push_back(-1);
		this->m_intersectionPointInObjectSegmentId.push_back(-1);
		this->m_intersectionPointInCutTriangleId.push_back(-1);
		this->m_intersectionPointInCutSegmentId.push_back(-1);
		this->m_intersectionPointAbove.push_back(pointsAbove[objectPointId]);
	}

	for(Segments_t::iterator CutSegmentIterator = CutSegments.begin();CutSegmentIterator != CutSegments.end();CutSegmentIterator++)
	{
		Segment_t CutCurrentSegment = *CutSegmentIterator;
		int CutSegmentId = std::distance(CutSegments.begin(),CutSegmentIterator);
		int CutTriangleId = CutSegmentId/3;

		if(objectAABBTree.do_intersect(CutCurrentSegment))
		{
			Intersections_t intersectionsOfCutCurrentSegmentWithObjectTriangles;
			objectAABBTree.all_intersections(CutCurrentSegment, std::back_inserter(intersectionsOfCutCurrentSegmentWithObjectTriangles));

			for(Intersections_t::iterator intersectionIterator = intersectionsOfCutCurrentSegmentWithObjectTriangles.begin();
					intersectionIterator != intersectionsOfCutCurrentSegmentWithObjectTriangles.end();intersectionIterator++)
			{
				Object_and_primitive_id_t currentIntersectionObjectAndPrimitiveId = *intersectionIterator;
				CGAL::Object currentIntersection = currentIntersectionObjectAndPrimitiveId.first;

				Point_t currentIntersectionPoint;
				if(CGAL::assign(currentIntersectionPoint,currentIntersection))
				{
					Primitive_id_t objectTrianglePrimitiveId = currentIntersectionObjectAndPrimitiveId.second;

					Triangles_t::iterator objectTriangleIterator = objectTrianglePrimitiveId;
					int objectTriangleId = std::distance(objectTriangles.begin(),objectTriangleIterator);

					if(!objectTriangleIsSeparated[objectTriangleId/4][objectTriangleId%4]&&!objectTrianglePartiallySeparated[objectTriangleId/4][objectTriangleId%4])
					{
						std::cout << endl << "Object Triangle " << objectTriangleId <<
						" intersects with cut but non of its edges are split - Code can't handle this, consider to work with smaller tetrahedra!" << endl;
					}

					int objectTetrahedronId = objectTriangleId/4;

					CutSegmentIdsThatIntersectWithObject.push_back(CutSegmentId);
					objectTetrahedronIdOfCutSegmentIdsThatIntersectWithObject.push_back(objectTetrahedronId);

					for(int above=0;above<2;above++)
					{
						this->m_intersectionPoints.push_back(currentIntersectionPoint);

						this->m_intersectionPointInObjectTetrahedronId.push_back(objectTetrahedronId);
						this->m_intersectionPointInObjectTriangleId.push_back(objectTriangleId);
						this->m_intersectionPointInObjectSegmentId.push_back(-1);
						this->m_intersectionPointInCutTriangleId.push_back(CutTriangleId);
						this->m_intersectionPointInCutSegmentId.push_back(CutSegmentId);
						this->m_intersectionPointAbove.push_back(above);
					}
				}
			}
		}
	}





	/*GET THE INTERSECTIONS OF THE OBJECT SEGMENTS WITH THE CUT TRIANGLES*/
	Ints_t CutTriangleSplitsObjectSegmentIds[numberOfCutTriangles];

	// the following loop improves the speed by only running over split segments of the object
	for(Ints_t::iterator objectSplitSegmentIterator=this->m_objectSplitSegmentIds.begin();objectSplitSegmentIterator!=this->m_objectSplitSegmentIds.end();objectSplitSegmentIterator++)
	{
		int objectSegmentId = *objectSplitSegmentIterator;
		Segment_t objectCurrentSegment = objectSegments[objectSegmentId];

		if(CutAABBTree.do_intersect(objectCurrentSegment))
		{
			Intersections_t intersectionsOfObjectCurrentSegmentWithCutTriangles;
			CutAABBTree.all_intersections(objectCurrentSegment, std::back_inserter(intersectionsOfObjectCurrentSegmentWithCutTriangles));

			for(Intersections_t::iterator intersectionIterator = intersectionsOfObjectCurrentSegmentWithCutTriangles.begin();
					intersectionIterator != intersectionsOfObjectCurrentSegmentWithCutTriangles.end();intersectionIterator++)
			{
				Object_and_primitive_id_t currentIntersectionObjectAndPrimitiveId = *intersectionIterator;
				CGAL::Object currentIntersection = currentIntersectionObjectAndPrimitiveId.first;
				Point_t currentIntersectionPoint;
				if(CGAL::assign(currentIntersectionPoint,currentIntersection))
				{
					Primitive_id_t CutTrianglePrimitiveId = currentIntersectionObjectAndPrimitiveId.second;

					Triangles_t::iterator CutTriangleIterator = CutTrianglePrimitiveId;
					int CutTriangleId = std::distance(CutTriangles.begin(),CutTriangleIterator);

					int objectTetrahedronId = objectSegmentId/6;

					CutTriangleSplitsObjectSegmentIds[CutTriangleId].push_back(objectSegmentId);
					for(int above=0;above<2;above++)
					{
						this->m_intersectionPoints.push_back(currentIntersectionPoint);

						this->m_intersectionPointInObjectTetrahedronId.push_back(objectTetrahedronId);
						this->m_intersectionPointInObjectTriangleId.push_back(-1);
						this->m_intersectionPointInObjectSegmentId.push_back(objectSegmentId);
						this->m_intersectionPointInCutTriangleId.push_back(CutTriangleId);
						this->m_intersectionPointInCutSegmentId.push_back(-1);
						this->m_intersectionPointAbove.push_back(above);
					}
				}
			}
		}
	}





	/* GET CUT POINTS IN THE OBJECT*/

	Points_t cutFront;

	for(Points_t::iterator cutPointsIterator=cutPoints.begin();cutPointsIterator!=cutPoints.end();cutPointsIterator++)
	{
		Point_t CutCurrentPoint = *cutPointsIterator;
		int currentPointId = std::distance(cutPoints.begin(),cutPointsIterator);

		int cutPointIsInObjectTetrahedronId;
		getTetrahedronThatSurroundsPoint(objectTriangles,CutCurrentPoint,cutPointIsInObjectTetrahedronId);

		// the following probably could be speeded up, by using another function in which we write the ids of the cuts that have a vertex with the current point id
		if(cutPointIsInObjectTetrahedronId!=-1)
		{
			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int pointId=0;pointId<3;pointId++)
				{
					if(currentPointId==CutConnectTriangle[CutTriangleId][pointId])
					{
						for(int above=0;above<2;above++)
						{
							this->m_intersectionPoints.push_back(CutCurrentPoint);

							this->m_intersectionPointInObjectTetrahedronId.push_back(cutPointIsInObjectTetrahedronId);
							this->m_intersectionPointInObjectTriangleId.push_back(-1);
							this->m_intersectionPointInObjectSegmentId.push_back(-1);
							this->m_intersectionPointInCutTriangleId.push_back(CutTriangleId);
							this->m_intersectionPointInCutSegmentId.push_back(-1);
							this->m_intersectionPointAbove.push_back(above);
						}
						for(int currentPointIdEdgeId=0;currentPointIdEdgeId<2;currentPointIdEdgeId++)
						{
							if(CutNeighbours[CutTriangleId][edgesOnTriangleNodes[pointId][currentPointIdEdgeId]]==-1)
							{
								cutFront.push_back(CutCurrentPoint);	// watch out - every node appears twice in the cut front
							}

							CutSegmentIdsThatIntersectWithObject
								.push_back(CutTriangleId*3 + edgesOnTriangleNodes[pointId][currentPointIdEdgeId]);
							objectTetrahedronIdOfCutSegmentIdsThatIntersectWithObject.push_back(cutPointIsInObjectTetrahedronId);
						}
					}
				}
			}
		}
	}





	/*GET THE PART OF THE CUT BOUNDARY WHICH IS IN THE CURRENT TETRAHEDRON*/
	for(Ints_t::iterator iterObjectTetrahedraIds=this->m_objectTetrahedraIdsThatIntersectWithCut.begin();
			iterObjectTetrahedraIds!=this->m_objectTetrahedraIdsThatIntersectWithCut.end();iterObjectTetrahedraIds++)
	{
		int objectTetrahedronId = *iterObjectTetrahedraIds;
		int internalId = std::distance(this->m_objectTetrahedraIdsThatIntersectWithCut.begin(),iterObjectTetrahedraIds);

		if(objectTetrahedronPartiallySeparated[objectTetrahedronId])
		{
			Ints_t Ids;
			getVectorIdsOfIntEntry(objectTetrahedronId,objectTetrahedronIdOfCutSegmentIdsThatIntersectWithObject,Ids);
			Ints_t currentObjectTetrahedronIntersectsWithCutSegmentIds;
			for(Ints_t::iterator iter=Ids.begin();iter!=Ids.end();iter++)
			{
				int currentId = *iter;
				currentObjectTetrahedronIntersectsWithCutSegmentIds.push_back(CutSegmentIdsThatIntersectWithObject[currentId]);
			}

			deleteMultipleInts(currentObjectTetrahedronIntersectsWithCutSegmentIds);
			for(Ints_t::iterator iterCutSegmentId=currentObjectTetrahedronIntersectsWithCutSegmentIds.begin();
					iterCutSegmentId!=currentObjectTetrahedronIntersectsWithCutSegmentIds.end();iterCutSegmentId++)
			{
				int globalCutSegmentId = *iterCutSegmentId;
				int cutTriangleId = globalCutSegmentId/3;
				int localSegmentId = globalCutSegmentId%3;
				if(CutNeighbours[cutTriangleId][localSegmentId]==-1 && globalCutSegmentId!=-1)
				{
                    this->m_cutBoundary.push_back(CutSegments[globalCutSegmentId]);
                    this->m_cutBoundaryTriangles.push_back(CutTriangles[cutTriangleId]);
                    this->m_cutBoundaryInObjectTetrahedronId.push_back(objectTetrahedronId);
					Matrix_3_3_t currentTransformationMatrix;
					int currentSegmentNodeIds[2] = {nodesOnTriangleEdges[localSegmentId][0],nodesOnTriangleEdges[localSegmentId][1]};
					getTransformationMatrixForOneElementOfCutboundary
						(currentSegmentNodeIds,CutTriangles[cutTriangleId],currentTransformationMatrix);

					this->m_cutBoundaryTransformationMatrices.push_back(currentTransformationMatrix);
				}
			}
		}
	}





	/* TRIANGULATION OF THE PARTS OF THE CUT THAT LIE IN EACH ELEMENT SEPARATELY
	 * Inner triangulation - the triangles that lie in the cut - NEARLY DONE, has to be added to the element above and below*/
	for(Ints_t::iterator objectTetrahedraIdsThatIntersectWithCutIterator = this->m_objectTetrahedraIdsThatIntersectWithCut.begin();
			objectTetrahedraIdsThatIntersectWithCutIterator!=this->m_objectTetrahedraIdsThatIntersectWithCut.end();
			objectTetrahedraIdsThatIntersectWithCutIterator++)
	{
		int objectTetrahedronId = *objectTetrahedraIdsThatIntersectWithCutIterator;

		for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
		{
			for(int above=0;above<2;above++)
			{
				Ints_t Ids;
				getVectorIdsOf6IntEntriesIn6Vectors(objectTetrahedronId,-2,-2,CutTriangleId,-2,above,
						this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
						this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
						this->m_intersectionPointAbove,Ids);
				Points_t currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId;
				for(Ints_t::iterator iter=Ids.begin();iter!=Ids.end();iter++)
				{
					int currentId = *iter;
					currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId.push_back(this->m_intersectionPoints[currentId]);
				}

				int numberOfCurrentPoints = std::distance(currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId.begin(),currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId.end());

				Ints_t IdsSortedPoints;
				// if we have more than three points then we need to order them (counter) clockwise,
				// s.t. we can build the triangulation and s.t. the triangle normals point in the same direction like the normal of the cut triangle
				if(numberOfCurrentPoints>3)
				{
					// push back the first point
					IdsSortedPoints.push_back(Ids[0]);

					// choose the second point
					int chosenPointId;
					Vector_t CutNormal;
					calculateNormalOfTriangle(CutTriangles[CutTriangleId],CutNormal);
					for(int pointId1=1;pointId1<numberOfCurrentPoints;pointId1++)
					{
						Vector_t v = currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[pointId1]-currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[0];
						Vector_t x = cross_product(v,CutNormal);
						int help=0;
						for(int pointId2=1;pointId2<numberOfCurrentPoints;pointId2++)
						{
							if(pointId1!=pointId2)
							{
								Vector_t w = currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[pointId2]-currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[0];
								help +=w*x>0;
							}
						}
						if(help==(numberOfCurrentPoints-2))
						{
							chosenPointId=pointId1;
						}
					}
					IdsSortedPoints.push_back(Ids[chosenPointId]);

					// use a set where the remaining points are saved
					Points_t helpPoints;
					Ints_t helpIds;
					for(int pointId=1;pointId<numberOfCurrentPoints;pointId++)
					{
						if(pointId!=chosenPointId)
						{
							helpPoints.push_back(currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[pointId]);
							helpIds.push_back(pointId);
						}
					}

					// attach the remaining points in clockwise order
					Vector_t v1 = currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[chosenPointId]-currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[0];
					double length_v1 = v1.squared_length();
					v1 = v1/sqrt(length_v1);

					std::vector<DoubleAndInt_t> angles;
					for(int pointId=0;pointId<numberOfCurrentPoints-2;pointId++)
					{
						DoubleAndInt_t currentAngleAndIndice;
						Vector_t v2 = helpPoints[pointId]-currentIntersectionPointsOfCutTriangleWithObjectTetrahedronId[0];
						double length_v2 = v2.squared_length();
						v2 = v2/sqrt(length_v2);
						currentAngleAndIndice.first = acos(v2*v1);
						currentAngleAndIndice.second = pointId;
						angles.push_back(currentAngleAndIndice);
					}
					std::sort(angles.begin(),angles.end(), this->compare_DoubleAndInt);

					for(std::vector<DoubleAndInt_t>::iterator it=angles.begin();it!=angles.end();it++)
					{
						DoubleAndInt_t currentAngleAndIndice = *it;
						int Id = currentAngleAndIndice.second;
						IdsSortedPoints.push_back(Ids[helpIds[Id]]);
					}
				}
				else
				{
					for(int i=0;i<std::distance(Ids.begin(),Ids.end());i++)
						IdsSortedPoints.push_back(Ids[i]);
				}

				// write the triangles in the structure
				for(int triangleId=0;triangleId<numberOfCurrentPoints-2;triangleId++)
				{
					ConnectTriangle_t currentConnectTriangleType(IdsSortedPoints[0],IdsSortedPoints[triangleId+1],IdsSortedPoints[triangleId+2]);

					// check whether the normal of the triangle points in the same direction - if not, turn it around
					secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
						(CutTriangles[CutTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

					this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
					this->m_TriangleOfTriangulationAbove.push_back(above);
					this->m_TriangleOfTriangulationInObjectTriangleId.push_back(-1);
					this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
				}
			}
		}
	}





	/* TRIANGULATION OF THE FACES OF THE OBJECT TETRAHEDRA */

	for(Ints_t::iterator objectTriangleIdsThatIntersectWithCutIterator=this->m_objectTriangleIdsThatIntersectWithCut.begin();
			objectTriangleIdsThatIntersectWithCutIterator!=this->m_objectTriangleIdsThatIntersectWithCut.end();
			objectTriangleIdsThatIntersectWithCutIterator++)
	{
		int objectTriangleId = *objectTriangleIdsThatIntersectWithCutIterator;
		int objectTetrahedronId = objectTriangleId/4;
		int objectLocalFaceId = objectTriangleId%4;

		// check which nodes lie above (which is the second entry) and below (first entry)
		Ints_t objectCurrentPointIdsBelowAndAbove[2];
		int numberOfNodesAbove=0;
		for(int pointId=0;pointId<3;pointId++)
		{
			int currentVertexId = objectConnectTriangle[objectTriangleId][pointId];
			objectCurrentPointIdsBelowAndAbove[pointsAbove[currentVertexId]].push_back(currentVertexId);
			numberOfNodesAbove += pointsAbove[currentVertexId];
		}

		Ints_t objectSegmentIds[numberOfCutTriangles];
		for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
		{
			// the following number are ALL split object segments, not only in the current object triangle
			int currentNumberOfSplitObjectSegments = std::distance(CutTriangleSplitsObjectSegmentIds[CutTriangleId].begin(),CutTriangleSplitsObjectSegmentIds[CutTriangleId].end());

			// if there are split object segments, get the number and id of split object segments of the current object triangle
			if(currentNumberOfSplitObjectSegments)
			{
				for(int edge=0;edge<3;edge++)
				{
					int currentObjectSegmentId = (objectTriangleId/4)*6+edgesOnTetrahedronFaces[objectTriangleId%4][edge];
					for(Ints_t::iterator iter=CutTriangleSplitsObjectSegmentIds[CutTriangleId].begin();iter!=CutTriangleSplitsObjectSegmentIds[CutTriangleId].end();iter++)
					{
						if(currentObjectSegmentId == *iter)
						{
							objectSegmentIds[CutTriangleId].push_back(currentObjectSegmentId);
						}
					}
				}
			}
		}

		// triangulize the completely cut triangles
		if(objectTriangleIsSeparated[objectTetrahedronId][objectLocalFaceId])
		{
			// this will be used if the a cut has two segments that intersect the object triangle
			Ints_t CurrentObjectSegmentIds;
			int useObjectSegment[2][2];
			for(int i=0;i<2;i++)
				for(int j=0;j<2;j++)
					useObjectSegment[i][j] = (i+j)%2;
			int useObjectSegmentId = 0;
			Ints_t CutTriangleIdThatCutsObjectSegment;

			// if the triangle is separated, then it has at least two split segments
			// for each of those intersections with a segment of the object we add a triangle
			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int above=0;above<2;above++)
				{
					int numberOfSplitEdgesByCut = std::distance(objectSegmentIds[CutTriangleId].begin(),objectSegmentIds[CutTriangleId].end());
					if(numberOfSplitEdgesByCut==1)
					{
						int currentObjectSegmentId = *objectSegmentIds[CutTriangleId].begin();
						if(CurrentObjectSegmentIds.empty() || CurrentObjectSegmentIds.back() != currentObjectSegmentId)
							CurrentObjectSegmentIds.push_back(currentObjectSegmentId);
						if(above==0)
							CutTriangleIdThatCutsObjectSegment.push_back(CutTriangleId);

						Ints_t Ids;
						getVectorIdsOf6IntEntriesIn6Vectors(objectTetrahedronId,-1,currentObjectSegmentId,CutTriangleId,-1,above,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids);

						Ints_t Ids2;
						getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleId,-2,above,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids2);

						// get the triangles
						for(int pointId=0;pointId<2;pointId++)
						{
							int currentVertexId = objectConnectSegment[currentObjectSegmentId][pointId];
							int currentVertexIsAbove = pointsAbove[currentVertexId];

							if(above==currentVertexIsAbove)
							{
								// new version
								ConnectTriangle_t currentConnectTriangleType(Ids2[0],currentVertexId,Ids[0]);

								// check whether the normal of the triangle points in the same direction - if not, turn it around
								secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
									(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

								this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
								this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
								this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
								this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
							}
						}
					}
					if(numberOfSplitEdgesByCut==2)
					{
						if(CurrentObjectSegmentIds.empty())
								CurrentObjectSegmentIds = objectSegmentIds[CutTriangleId];
						Ints_t Ids;
						for(int localSegmentId=0;localSegmentId<2;localSegmentId++)
						{
							getVectorIdsOf6IntEntriesIn6Vectors(-2,-1,CurrentObjectSegmentIds[localSegmentId],CutTriangleId,-1,above,
									this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
									this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
									this->m_intersectionPointAbove,Ids);
						}
						for(int pointId=0;pointId<2;pointId++)
						{
							int currentVertexId = objectConnectSegment[CurrentObjectSegmentIds[0]][pointId];	// "+" the id in CurrentObjectSegmentIds could also be 1, but then we would have to change the id in "* to 0
							int currentVertexIsAbove = pointsAbove[currentVertexId];

							if(above==currentVertexIsAbove)
							{
								// new, better, version
								ConnectTriangle_t currentConnectTriangleType(Ids[0],Ids[1],currentVertexId);

								// check whether the normal of the triangle points in the same direction - if not, turn it around
								secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
									(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

								this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
								this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
								this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
								this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
							}
						}

						if(above==(numberOfNodesAbove-1))
						{
							Ints_t Ids2;
							getVectorIdsOf6IntEntriesIn6Vectors(-2,-1,CurrentObjectSegmentIds[1],CutTriangleId,-1,above,// "*" the id in CurrentObjectSegmentIds could also be 0, but then we would have to change the id in "+" to 1
									this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
									this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
									this->m_intersectionPointAbove,Ids2);

							int objectPointId0 = *objectCurrentPointIdsBelowAndAbove[numberOfNodesAbove-1].begin();
							int objectPointId1 = *(objectCurrentPointIdsBelowAndAbove[numberOfNodesAbove-1].begin()+1);

							// new, better, version
							ConnectTriangle_t currentConnectTriangleType(Ids2[0],objectPointId0,objectPointId1);

							// check whether the normal of the triangle points in the same direction - if not, turn it around
							secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
								(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

							this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
							this->m_TriangleOfTriangulationAbove.push_back(numberOfNodesAbove-1);
							this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
							this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
						}
					}
				}
			}

			double area[2];
			area[0]=0;
			area[1]=0;
			for(int temporaryTriangulationId=0;temporaryTriangulationId<2;temporaryTriangulationId++)
			{
				for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
				{
					Ints_t Ids;
					getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleId,-2,numberOfNodesAbove-1,	// here we just need to make sure, that we don't count a node twice, that s why we set above = 0
							this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
							this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
							this->m_intersectionPointAbove,Ids);

					int numberOfCurrentIntersectionPoints = std::distance(Ids.begin(),Ids.end());
					if(numberOfCurrentIntersectionPoints==2)
					{
						int currentVertexId = -1;
						for(int pointId=0;pointId<2;pointId++)
						{
							int Id;
							getVectorIdOfIntEntry(objectConnectSegment[CurrentObjectSegmentIds[useObjectSegment[temporaryTriangulationId][0]]][pointId],
									objectCurrentPointIdsBelowAndAbove[numberOfNodesAbove-1],Id);
							if(Id!=-1)
							{
								currentVertexId = objectConnectSegment[CurrentObjectSegmentIds[useObjectSegment[temporaryTriangulationId][0]]][pointId];
							}
						}
						if(currentVertexId==-1)
							std::cout << "ERROR - wrong thoughts!" << endl;
						Triangle_t currentTriangle = Triangle_t(this->m_intersectionPoints[Ids[0]],
								this->m_intersectionPoints[Ids[1]],this->m_objectPoints[currentVertexId]);
						area[temporaryTriangulationId] += sqrt(currentTriangle.squared_area());
					}
					if(std::distance(CutTriangleIdThatCutsObjectSegment.begin(),CutTriangleIdThatCutsObjectSegment.end())==2)
					{
						Ints_t Ids2;
						getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleIdThatCutsObjectSegment[useObjectSegment[temporaryTriangulationId][1]],-2,numberOfNodesAbove-1,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids2);

						Ints_t currentPointIds = objectCurrentPointIdsBelowAndAbove[numberOfNodesAbove-1];
						Triangle_t currentTriangle = Triangle_t(this->m_intersectionPoints[currentPointIds[0]],
								this->m_intersectionPoints[currentPointIds[1]],this->m_intersectionPoints[Ids2[0]]);
						area[temporaryTriangulationId] += sqrt(currentTriangle.squared_area());
					}
				}
			}

			if(area[1]<area[0])
			{
				useObjectSegmentId = 1;
			}

			// get the triangulation which derives from cut triangles whos segments have two intersections with the object triangle
			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int above=0;above<2;above++)
				{
					Ints_t Ids;
					getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleId,-2,above,
							this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
							this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
							this->m_intersectionPointAbove,Ids);

					int numberOfCurrentIntersectionPoints = std::distance(Ids.begin(),Ids.end());
					if(numberOfCurrentIntersectionPoints==2)
					{
						for(int pointId=0;pointId<2;pointId++)
						{
							int currentVertexId = objectConnectSegment[CurrentObjectSegmentIds[useObjectSegment[useObjectSegmentId][0]]][pointId];
							int currentVertexIsAbove = pointsAbove[currentVertexId];

							if(above==currentVertexIsAbove)
							{
								// new, better, version
								ConnectTriangle_t currentConnectTriangleType(Ids[0],Ids[1],currentVertexId);

								// check whether the normal of the triangle points in the same direction - if not, turn it around
								secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
									(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

								this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
								this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
								this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
								this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
							}
						}
					}
				}
			}

			for(int above=0;above<2;above++)
			{
				if(above==numberOfNodesAbove-1)//if(above==m_intersectionPointAbove[currentPointIds[0]])//
				{
					if(std::distance(CutTriangleIdThatCutsObjectSegment.begin(),CutTriangleIdThatCutsObjectSegment.end())==2)
					{
						Ints_t Ids2;
						getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleIdThatCutsObjectSegment[useObjectSegment[useObjectSegmentId][1]],-2,above,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids2);

						Ints_t currentPointIds = objectCurrentPointIdsBelowAndAbove[numberOfNodesAbove-1];

						// new, better, version
						ConnectTriangle_t currentConnectTriangleType(currentPointIds[0],currentPointIds[1],Ids2[0]);

						// check whether the normal of the triangle points in the same direction - if not, turn it around
						secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
							(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

						this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
						this->m_TriangleOfTriangulationAbove.push_back(numberOfNodesAbove-1);
						this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
						this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
					}
				}
			}
		}


		// triangulize the partially cut triangles
		if(objectTrianglePartiallySeparated[objectTetrahedronId][objectLocalFaceId])
		{
			// this will be used if the a cut has two segments that intersect the object triangle
			Ints_t CurrentObjectSegmentIds;
			Ints_t CutTriangleIdThatCutsObjectSegment;

			// if the triangle is separated, then it has at least two split segments
			// for each of those intersections with a segment of the object we add a triangle
			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int above=0;above<2;above++)
				{
					int numberOfSplitEdgesByCut = std::distance(objectSegmentIds[CutTriangleId].begin(),objectSegmentIds[CutTriangleId].end());
					if(numberOfSplitEdgesByCut==1)
					{
						CurrentObjectSegmentIds = objectSegmentIds[CutTriangleId];
						int currentObjectSegmentId = CurrentObjectSegmentIds[0];
						CutTriangleIdThatCutsObjectSegment.push_back(CutTriangleId);

						Ints_t Ids;
						getVectorIdsOf6IntEntriesIn6Vectors(objectTetrahedronId,-1,currentObjectSegmentId,CutTriangleId,-1,above,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids);

						Ints_t Ids2;
						getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleId,-2,above,
								this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
								this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
								this->m_intersectionPointAbove,Ids2);

						// get the triangles
						for(int pointId=0;pointId<2;pointId++)
						{
							int currentVertexId = objectConnectSegment[currentObjectSegmentId][pointId];
							int currentVertexIsAbove = pointsAbove[currentVertexId];

							if(above==currentVertexIsAbove)
							{
								// new, better, version
								ConnectTriangle_t currentConnectTriangleType(Ids2[0],currentVertexId,Ids[0]);

								// check whether the normal of the triangle points in the same direction - if not, turn it around
								secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
									(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

								this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
								this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
								this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
								this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
							}
						}
					}
				}
			}

			// get the triangulation which derives from cut triangles whos segments have two intersections with the object triangle
			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int above=0;above<2;above++)
				{
					Ints_t Ids;
					getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,CutTriangleId,-2,above,
							this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
							this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
							this->m_intersectionPointAbove,Ids);

					int numberOfCurrentIntersectionPoints = std::distance(Ids.begin(),Ids.end());
					if(numberOfCurrentIntersectionPoints==2)
					{
						for(int pointId=0;pointId<2;pointId++)
						{
							int currentVertexId = objectConnectSegment[CurrentObjectSegmentIds[0]][pointId];
							int currentVertexIsAbove = pointsAbove[currentVertexId];

							if(above==currentVertexIsAbove)
							{
								// new, better, version
								ConnectTriangle_t currentConnectTriangleType(Ids[0],Ids[1],currentVertexId);

								// check whether the normal of the triangle points in the same direction - if not, turn it around
								secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
									(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

								this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
								this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
								this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
								this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
							}
						}
					}
				}
			}

			for(int CutTriangleId=0;CutTriangleId<numberOfCutTriangles;CutTriangleId++)
			{
				for(int above=0;above<2;above++)
				{
					int IdFirstPoint;
					bool theCurrentCutIsAtTheCutFront=false;
					for(int edgeId=0;edgeId<3;edgeId++)
					{
						if(CutNeighbours[CutTriangleId][edgeId]==-1)
						{
							Ints_t Ids;
							getVectorIdsOf6IntEntriesIn6Vectors(-2,objectTriangleId,-1,-2,CutTriangleId*3+edgeId,above,
									this->m_intersectionPointInObjectTetrahedronId,this->m_intersectionPointInObjectTriangleId,this->m_intersectionPointInObjectSegmentId,
									this->m_intersectionPointInCutTriangleId,this->m_intersectionPointInCutSegmentId,
									this->m_intersectionPointAbove,Ids);

							if(std::distance(Ids.begin(),Ids.end())>0)
							{
								IdFirstPoint = Ids[0];
								theCurrentCutIsAtTheCutFront = true;
								cutFront.push_back(this->m_intersectionPoints[Ids[0]]);	// watch out every node appears twice because of the above!
							}
						}
					}

					if(theCurrentCutIsAtTheCutFront)
					{
						int Id1;

						for(int edgeId=0;edgeId<3;edgeId++)
						{
							for(int iter=0;iter<std::distance(this->m_objectSplitSegmentIds.begin(),this->m_objectSplitSegmentIds.end());iter++)
							{
								if(((objectTriangleId/4)*6+edgesOnTetrahedronFaces[objectTriangleId%4][edgeId])==this->m_objectSplitSegmentIds[iter])
								{
									Id1=iter;
								}
							}
						}

						int globalSegmentNumber = this->m_objectSplitSegmentIds[Id1];
						int NodeId0 = objectConnectSegment[globalSegmentNumber][0];
						int NodeId1 = objectConnectSegment[globalSegmentNumber][1];

						int NodeIdSpecial;
						for(int nodeId=0;nodeId<3;nodeId++)
						{
							if(NodeId0 != objectConnectTriangle[objectTriangleId][nodeId] && NodeId1 != objectConnectTriangle[objectTriangleId][nodeId])
								NodeIdSpecial = objectConnectTriangle[objectTriangleId][nodeId];
						}

						int currentVertexIsAbove = pointsAbove[NodeId0];

						if(above==currentVertexIsAbove)
						{
							// new, better, version
							ConnectTriangle_t currentConnectTriangleType(IdFirstPoint,NodeId0,NodeIdSpecial);

							// check whether the normal of the triangle points in the same direction - if not, turn it around
							secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
								(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType);

							this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType);
							this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
							this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
							this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
						}

						currentVertexIsAbove = pointsAbove[NodeId1];

						if(above==currentVertexIsAbove)
						{
							// new, better, version
							ConnectTriangle_t currentConnectTriangleType2(IdFirstPoint,NodeId1,NodeIdSpecial);

							// check whether the normal of the triangle points in the same direction - if not, turn it around
							secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
								(objectTriangles[objectTriangleId],this->m_intersectionPoints,currentConnectTriangleType2);

							this->m_triangulationOfCutObjectTetrahedronFacesConnectTriangle.push_back(currentConnectTriangleType2);
							this->m_TriangleOfTriangulationAbove.push_back(currentVertexIsAbove);
							this->m_TriangleOfTriangulationInObjectTriangleId.push_back(objectTriangleId);
							this->m_TriangleOfTriangulationInObjectTetrahedronId.push_back(objectTetrahedronId);
						}
					}
				}
			}
		}
	}
}


void CutModelTopology_CGAL::getTetConnectFromVTKCellArray(vtkCellArray* myArray, int numberOfCellNodes, int ObjectConnectTetrahedra[][4], vtkIdType& numberOfCells, int &numberOfObjectNodes)
{
	vtkIdType* currentElement = new vtkIdType[numberOfCellNodes]; //noch zu löschen!!! da "new" verwendet wurde!
	vtkIdType numberOfCellNodes_vtkIdType = numberOfCellNodes;
	numberOfObjectNodes=0;
	for(vtkIdType elementId=0;elementId<numberOfCells*(numberOfCellNodes+1);elementId+=(numberOfCellNodes+1))
	{
		myArray->GetCell(elementId,numberOfCellNodes_vtkIdType,currentElement);
		for(int i=0;i<numberOfCellNodes;i++)
		{
			ObjectConnectTetrahedra[int(elementId/(numberOfCellNodes+1))][i] = currentElement[i];
			if(numberOfObjectNodes<currentElement[i]+1)
				numberOfObjectNodes = currentElement[i]+1;
		}
	}
}

// I would like to put the function above and below together, but I am having problems with the delivery of the output "connect"
void CutModelTopology_CGAL::getTriangleConnectFromVTKCellArray(vtkCellArray* myArray, int numberOfCellNodes, int ObjecttriangleConnect[][3], vtkIdType& numberOfCells, int &numberOfObjectNodes)
{
	vtkIdType* currentElement = new vtkIdType[numberOfCellNodes]; //noch zu löschen!!! da "new" verwendet wurde!
	vtkIdType numberOfCellNodes_vtkIdType = numberOfCellNodes;
	numberOfObjectNodes=0;
	for(vtkIdType elementId=0;elementId<numberOfCells*(numberOfCellNodes+1);elementId+=(numberOfCellNodes+1))
	{
		myArray->GetCell(elementId,numberOfCellNodes_vtkIdType,currentElement);
		for(int i=0;i<numberOfCellNodes;i++)
		{
			ObjecttriangleConnect[int(elementId/(numberOfCellNodes+1))][i] = currentElement[i];
			if(numberOfObjectNodes<currentElement[i]+1)
				numberOfObjectNodes = currentElement[i]+1;
		}
	}
}


void CutModelTopology_CGAL::getPointsFromVTKCellArray(int numberOfPoints, vtkPoints* myPoints, double factor, double pointsDouble[][3], Points_t& points)
{
	double currentPointDouble[3];
	for(unsigned int nodeId=0;nodeId<numberOfPoints;nodeId++)
	{
		myPoints->GetPoint(nodeId, currentPointDouble);
		for(int i=0;i<3;i++)
		{
			currentPointDouble[i]=currentPointDouble[i]*factor;	// this helps to visualize the points later as 3D glyphs in paraview
			pointsDouble[nodeId][i] = currentPointDouble[i];
		}
		points.push_back(Point_t(currentPointDouble[0],currentPointDouble[1],currentPointDouble[2]));
	}
}


void CutModelTopology_CGAL::getConnectTriangleOfConnectTet(int numberOfCells, int connectTetrahedra[][4], int nodesOnTetrahedronFaces[][3], int ConnectTriangle[][3])
{
	for(int elementId=0;elementId<numberOfCells;elementId++)
	{
		for(int faceId=0;faceId<4;faceId++)
		{
			for(int pointId=0;pointId<3;pointId++)
			{
				ConnectTriangle[elementId*4+faceId][pointId] = connectTetrahedra[elementId][nodesOnTetrahedronFaces[faceId][pointId]];
			}
		}
	}
}


void CutModelTopology_CGAL::getConnectSegmentOfConnectTet(int numberOfCells, int connectTetrahedra[][4], int nodesOnTetrahedronEdges[][2], int ConnectSegment[][2])
{
	int numberOfEdgeNodes = 2;
	int numberOfEdges = 6;
	for(int elementId=0;elementId<numberOfCells;elementId++)
	{
		for(int edgeId=0;edgeId<numberOfEdges;edgeId++)
		{
			for(int pointId=0;pointId<numberOfEdgeNodes;pointId++)
			{
				ConnectSegment[elementId*numberOfEdges+edgeId][pointId] =
						connectTetrahedra[elementId][nodesOnTetrahedronEdges[edgeId][pointId]];
			}
		}
	}
}


void CutModelTopology_CGAL::getConnectSegmentOfConnectTriangle(int numberOfCells, int ConnectTriangle[][3], int ConnectSegment[][2])
{
	int numberOfEdgeNodes = 2;
	int numberOfEdges = 3;
	int nodesOnTriangleEdges[3][2] = {{0,1},{1,2},{0,2}};
	for(int elementId=0;elementId<numberOfCells;elementId++)
	{
		for(int edgeId=0;edgeId<numberOfEdges;edgeId++)
		{
			for(int pointId=0;pointId<numberOfEdgeNodes;pointId++)
			{
				ConnectSegment[elementId*numberOfEdges+edgeId][pointId] =
						ConnectTriangle[elementId][nodesOnTriangleEdges[edgeId][pointId]];
			}
		}
	}
}


void CutModelTopology_CGAL::getTrianglesFromConnectTriangleAndPointsDouble(int numberOfCells, int ConnectTriangle[][3], double pointsDouble[][3], std::vector<Triangle_t> &triangles)
{
	double currentPointsDouble[3][3];
	for(int cellId=0;cellId<numberOfCells;cellId++)
	{
		for(int j=0;j<3;j++)
		{
			for(int i=0;i<3;i++)
			{
				currentPointsDouble[i][j] = pointsDouble[ConnectTriangle[cellId][j]][i];
			}
		}
		triangles.push_back(Triangle_t(	Point_t(currentPointsDouble[0][0],currentPointsDouble[1][0],currentPointsDouble[2][0]),
										Point_t(currentPointsDouble[0][1],currentPointsDouble[1][1],currentPointsDouble[2][1]),
										Point_t(currentPointsDouble[0][2],currentPointsDouble[1][2],currentPointsDouble[2][2])));
	}
}


void CutModelTopology_CGAL::getSegmentsFromConnectSegmentAndPointsDouble(int numberOfCells, int ConnectSegment[][2], double pointsDouble[][3], Segments_t &segments)
{
	double currentPointsDouble[3][2];
	for(int cellId=0;cellId<numberOfCells;cellId++)
	{
		for(int j=0;j<2;j++)
		{
			for(int i=0;i<3;i++)
			{
				currentPointsDouble[i][j] = pointsDouble[ConnectSegment[cellId][j]][i];
			}
		}
		segments.push_back(Segment_t(	Point_t(currentPointsDouble[0][0],currentPointsDouble[1][0],currentPointsDouble[2][0]),
										Point_t(currentPointsDouble[0][1],currentPointsDouble[1][1],currentPointsDouble[2][1])));
	}
}


void CutModelTopology_CGAL::getTetrahedronThatSurroundsPoint(Triangles_t &objectTriangles,Point_t point,int &pointIsInObjectTetrahedronId)
{
	Triangle_t objectCurrentTriangle;

	pointIsInObjectTetrahedronId = -1;

	int numberOfObjectTriangles = std::distance(objectTriangles.begin(),objectTriangles.end());
	int numberOfObjectTetrahedra = numberOfObjectTriangles/4;
	for(int objectTetrahedronId=0;objectTetrahedronId<numberOfObjectTetrahedra;objectTetrahedronId++)
	{
		int helper=0;
		for(int faceId=0;faceId<4;faceId++)
		{
			objectCurrentTriangle = objectTriangles[objectTetrahedronId*4+faceId];
			Vector_t normalOfCurrentFace;
			calculateNormalOfTriangle(objectCurrentTriangle,normalOfCurrentFace);

			Vector_t connectionVectorCutCurrentPointAndObjectCurrentTriangleVertex = point-objectTriangles[objectTetrahedronId*4+faceId].vertex(0);
			helper += (normalOfCurrentFace*connectionVectorCutCurrentPointAndObjectCurrentTriangleVertex)<=0;
		}
		if(helper==4)
		{
			pointIsInObjectTetrahedronId=objectTetrahedronId;
		}
	}
}


void CutModelTopology_CGAL::isPointAbove(Points_t &points, Tree_t &CutAABBTree, Triangles_t &CutTriangles, int* pointsAbove, Ints_t pointIdsAbove)
{
	for(Points_t::iterator pointsIterator=points.begin();pointsIterator!=points.end();++pointsIterator)
	{
		Point_t currentPoint = *pointsIterator;
		unsigned int ObjectPointId = std::distance(points.begin(),pointsIterator);

		Point_and_primitive_id_t closestPointOnCutANDidOfclosestTriangleOfCut = CutAABBTree.closest_point_and_primitive(currentPoint);
		Primitive_id_t PrimitiveIdOfclosestTriangleOfCut = closestPointOnCutANDidOfclosestTriangleOfCut.second; // closest primitive id
		Point_t closestPointOnCut = closestPointOnCutANDidOfclosestTriangleOfCut.first;
		std::vector<Triangle_t>::iterator IteratorOfclosestTriangleOfCut = PrimitiveIdOfclosestTriangleOfCut;
		unsigned int IdOfclosestTriangleOfCut = std::distance(CutTriangles.begin(),IteratorOfclosestTriangleOfCut);

		Vector_t currentNormal;
		Triangle_t currentTriangle = *PrimitiveIdOfclosestTriangleOfCut;
		calculateNormalOfTriangle(currentTriangle, currentNormal);

		Vector_t connectionVectorBetweenPointAndReferencePoint = currentPoint-closestPointOnCut;
		pointsAbove[ObjectPointId]  = ((currentNormal*connectionVectorBetweenPointAndReferencePoint)>0);
//			-((currentNormal*connectionVectorBetweenPointAndReferencePoint)==0);
		if(pointsAbove[ObjectPointId]==1)
			pointIdsAbove.push_back(ObjectPointId);
	}
}


int CutModelTopology_CGAL::pointAbove(Point_t &point)
{
	Tree_t CutAABBTree(m_cutTriangles.begin(),m_cutTriangles.end());

	Point_and_primitive_id_t closestPointOnCutANDidOfclosestTriangleOfCut = CutAABBTree.closest_point_and_primitive(point);
	Primitive_id_t PrimitiveIdOfclosestTriangleOfCut = closestPointOnCutANDidOfclosestTriangleOfCut.second; // closest primitive id
	Point_t closestPointOnCut = closestPointOnCutANDidOfclosestTriangleOfCut.first;
	std::vector<Triangle_t>::iterator IteratorOfclosestTriangleOfCut = PrimitiveIdOfclosestTriangleOfCut;
	unsigned int IdOfclosestTriangleOfCut = std::distance(m_cutTriangles.begin(),IteratorOfclosestTriangleOfCut);

	Vector_t currentNormal;
	Triangle_t currentTriangle = *PrimitiveIdOfclosestTriangleOfCut;
	calculateNormalOfTriangle(currentTriangle, currentNormal);

	Vector_t connectionVectorBetweenPointAndReferencePoint = point-closestPointOnCut;
	return ((currentNormal*connectionVectorBetweenPointAndReferencePoint)>0);
//	  -(currentNormal*connectionVectorBetweenPointAndReferencePoint)==0);
}


void CutModelTopology_CGAL::pointsAbove(Points_t &points,Ints_t &pointsAbove)
{
	Tree_t CutAABBTree(m_cutTriangles.begin(),m_cutTriangles.end());

	for(Points_t::iterator iter=points.begin();iter!=points.end();iter++)
	{
		Point_t currentPoint = *iter;
		int currentPointId = std::distance(points.begin(),iter);
		Point_and_primitive_id_t closestPointOnCutANDidOfclosestTriangleOfCut = CutAABBTree.closest_point_and_primitive(currentPoint);
		Primitive_id_t PrimitiveIdOfclosestTriangleOfCut = closestPointOnCutANDidOfclosestTriangleOfCut.second; // closest primitive id
		Point_t closestPointOnCut = closestPointOnCutANDidOfclosestTriangleOfCut.first;
		std::vector<Triangle_t>::iterator IteratorOfclosestTriangleOfCut = PrimitiveIdOfclosestTriangleOfCut;
		unsigned int IdOfclosestTriangleOfCut = std::distance(m_cutTriangles.begin(),IteratorOfclosestTriangleOfCut);

		Vector_t currentNormal;
		Triangle_t currentTriangle = *PrimitiveIdOfclosestTriangleOfCut;
		calculateNormalOfTriangle(currentTriangle, currentNormal);
		double squared_length_currentNormal = currentNormal.squared_length();
		Vector_t currentNormalWithLengthOne = currentNormal/sqrt(squared_length_currentNormal);

		Vector_t connectionVectorBetweenPointAndReferencePoint = currentPoint-closestPointOnCut;
		pointsAbove.push_back(((currentNormal*connectionVectorBetweenPointAndReferencePoint)>0));
//		  -((currentNormalWithLengthOne*connectionVectorBetweenPointAndReferencePoint<2e-008)*
//				  (currentNormalWithLengthOne*connectionVectorBetweenPointAndReferencePoint>-2e-008)*5));
	}
}

// THE FOLLOWING FUNCTION CAN BE DELETED
void CutModelTopology_CGAL::pointsInObjectTetrahedronNumber(std::string RefGridFilename, Ints_t &pointsReferenceSolutionInObjectTetrahedronId, Ints_t &pointsReferenceSolutionAbove, Doubles_t &pointsReferenceSolution_r, Doubles_t &pointsReferenceSolution_phi)
{
	// define the reader, the grid, the points and the array, to receive information read from the vtk file
	vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
	reader->SetFileName(RefGridFilename.c_str());
	reader->Update();
	vtkUnstructuredGrid* myGrid = reader->GetOutput();
	vtkPoints* myPoints = myGrid->GetPoints();
	vtkCellArray* myArray = myGrid->GetCells();

	// get the connect of the tetrahedra grid
	vtkIdType numberOfCells = myGrid->GetNumberOfCells();
	vtkIdType numberOfPoints = myPoints->GetNumberOfPoints();
	std::vector<Ints_t> tetrahedraConnectedToPointId;
	tetrahedraConnectedToPointId.resize((int)numberOfPoints);
	int objectConnectTetrahedra[int(numberOfCells)][4];
	int numberOfObjectPoints;
	vtkIdType* currentElement = new vtkIdType[4];
	vtkIdType numberOfCellNodes = 4;
	for(vtkIdType elementId=0;elementId<numberOfCells*(numberOfCellNodes+1);elementId+=(numberOfCellNodes+1))
	{
		myArray->GetCell(elementId,numberOfCellNodes,currentElement);
		for(int i=0;i<numberOfCellNodes;i++)
		{
			int currentNodeId = (int)currentElement[i];
			objectConnectTetrahedra[int(elementId/(numberOfCellNodes+1))][i] = currentNodeId;
			tetrahedraConnectedToPointId[currentNodeId].push_back((int)elementId/(numberOfCellNodes+1));
		}
	}
	
	Points_t points;
	points.resize(numberOfPoints);
	pointsReferenceSolution_r.resize(numberOfPoints);
	pointsReferenceSolution_phi.resize(numberOfPoints);
	
	double currentPointDouble[3];
	for(unsigned int nodeId=0;nodeId<numberOfPoints;nodeId++)
	{
		myPoints->GetPoint(nodeId, currentPointDouble);
		Point_t currentPoint(currentPointDouble[0],currentPointDouble[1],currentPointDouble[2]);
		points[nodeId] = currentPoint;
		getRAndPhiOfPointDependentOnCutFront(currentPoint,pointsReferenceSolution_r[nodeId],pointsReferenceSolution_phi[nodeId]);
	}
	
	
	pointsAbove(points,pointsReferenceSolutionAbove);
	
	for(int pointId=0;pointId<numberOfPoints;pointId++)
	{
		if(pointsReferenceSolutionAbove[pointId]<-1 && (points[pointId])[0]<=0.015)
		{
		  for(int i=0;i<4;i++)
		  {
		    int currentElementId = (tetrahedraConnectedToPointId[pointId])[0];
		    int currentNodeId = objectConnectTetrahedra[currentElementId][i];
		    int currentPointOfTetIsAbove = pointsReferenceSolutionAbove[currentNodeId];
		    if(currentPointOfTetIsAbove>-1)
		      pointsReferenceSolutionAbove[pointId] = currentPointOfTetIsAbove;
		  }
		  pointsReferenceSolution_phi[pointId] = atan(1.0)*4*(pointsReferenceSolutionAbove[pointId]*2-1);
		}
		else if(pointsReferenceSolutionAbove[pointId]<-1)
		{
			pointsReferenceSolutionAbove[pointId] = pointsReferenceSolution_phi[pointId]>0;
		}
	}

    pointsReferenceSolutionInObjectTetrahedronId.resize(numberOfPoints);
	for(unsigned int currentPointId=0;currentPointId<numberOfPoints;currentPointId++)
	{
		Point_t currentPoint = points[currentPointId];

		int currentPointIsInObjectTetrahedronId;
		getTetrahedronThatSurroundsPoint(this->m_objectTriangles,currentPoint,currentPointIsInObjectTetrahedronId);

		if(currentPointIsInObjectTetrahedronId<0)
			std::cout << "ERROR in CutModelTopology, the node id " << currentPointId << " can't be related to a tetrahedron in the current object!" << endl;
		if(pointsReferenceSolutionAbove[currentPointId]<-1)
			std::cout << "ERROR in CutModelTopology, the node " << currentPointId << " lies ON cut!" << endl;
		pointsReferenceSolutionInObjectTetrahedronId[currentPointId] = currentPointIsInObjectTetrahedronId;
	}
}


double CutModelTopology_CGAL::getROfPoint(Point_t point)
{
	// we could also choose to check the distance to the tree that derives from the cut direction - maybe we even have
	// to build a triangle tree of all the triangles involved in the element
	TreeSegment_t CutBoundaryAABBTree(m_cutBoundary.begin(),m_cutBoundary.end());
	TreeSegment_Point_and_primitive_id_t closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary = CutBoundaryAABBTree.closest_point_and_primitive(point);
	Point_t closestPoint = closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary.first;

	Vector_t connectionVectorBetweenPointAndReferencePoint = point-closestPoint;
	double squaredLengthOfVector = connectionVectorBetweenPointAndReferencePoint.squared_length();
	double r = sqrt(squaredLengthOfVector);
	return r;
}


void CutModelTopology_CGAL::getTransformationMatrixForOneElementOfCutboundary
		(int* cutBoundarySegmentNodeIdsOfCutBoundaryTriangle,Triangle_t &cutBoundaryTriangle,Matrix_3_3_t &transformationMatrix)
{
	// compare to the matlab code of Christoph Paulus, in which the coordinate system is mostly called a,b,c
	Vector_t b;
	calculateNormalOfTriangle(cutBoundaryTriangle, b);
	double squaredLengthOfVector = b.squared_length();
	b = b/sqrt(squaredLengthOfVector);

	Vector_t c = cutBoundaryTriangle.vertex(cutBoundarySegmentNodeIdsOfCutBoundaryTriangle[0])-cutBoundaryTriangle.vertex(cutBoundarySegmentNodeIdsOfCutBoundaryTriangle[1]);
	squaredLengthOfVector = c.squared_length();
	c = c/sqrt(squaredLengthOfVector);
	Vector_t a = cross_product(b,c);
	int nodeIdThatIsNotInCutBoundary = 3-(cutBoundarySegmentNodeIdsOfCutBoundaryTriangle[0]+cutBoundarySegmentNodeIdsOfCutBoundaryTriangle[1]);
	Vector_t vectorFromPointInTriangleToPointOfSegment = cutBoundaryTriangle.vertex(cutBoundarySegmentNodeIdsOfCutBoundaryTriangle[0])-cutBoundaryTriangle.vertex(nodeIdThatIsNotInCutBoundary);
	int aPointsInRightDirection = 2*((vectorFromPointInTriangleToPointOfSegment*a)>0)-1;
	if(aPointsInRightDirection == -1)
		std::cout << "WATCH OUT: one of the cut triangles is not constructed in the way it is supposed to be!\n"
				"For further information please contact Christoph Paulus\n";

	for(int i=0;i<3;i++)
	{
		// transformationMatrix[i][0] = aPointsInRightDirection*a[i];
		// transformationMatrix[i][1] = b[i];
		// transformationMatrix[i][2] = aPointsInRightDirection*c[i];
		transformationMatrix(i,0) = aPointsInRightDirection*a[i];
		transformationMatrix(i,1) = b[i];
		transformationMatrix(i,2) = aPointsInRightDirection*c[i];
	}
}


void CutModelTopology_CGAL::getR_Phi_XAndYOfPointsDependentOnCutFront(Points_t &points,
		Doubles_t &r, Doubles_t &phi, Doubles_t &x, Doubles_t &y, Ints_t &idClosestCutBoundarySegment)
{
	for(Points_t::iterator iter=points.begin();iter!=points.end();iter++)
	{
		int currentPointId = std::distance(points.begin(),iter);
		Point_t point = *iter;
			
		// we could also choose to check the distance to the tree that derives from the cut direction - maybe we even have
		// to build a triangle tree of all the triangles involved in the element
		TreeSegment_t CutBoundaryAABBTree(m_cutBoundary.begin(),m_cutBoundary.end());
		TreeSegment_Point_and_primitive_id_t closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary = CutBoundaryAABBTree.closest_point_and_primitive(point);
        Point_t closestPoint = closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary.first;

        Vector_t connectionVectorBetweenPointAndReferencePoint = point-closestPoint;
		double squaredLengthOfVector = connectionVectorBetweenPointAndReferencePoint.squared_length();
		r[currentPointId] = sqrt(squaredLengthOfVector);

		TreeSegment_Primitive_id_t PrimitiveIdOfClosestSegmentOfCutBoundary = closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary.second;
		Segments_t::iterator IteratorOfClosestSegmentOfCutBoundary = PrimitiveIdOfClosestSegmentOfCutBoundary;
		unsigned int IdOfClosestSegmentOfCutBoundary = std::distance(m_cutBoundary.begin(),IteratorOfClosestSegmentOfCutBoundary);

        idClosestCutBoundarySegment[currentPointId] = IdOfClosestSegmentOfCutBoundary;

		double currentX = 0;
		for(int i=0;i<3;i++)
			currentX += connectionVectorBetweenPointAndReferencePoint[i]
				 // *(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])[i][0];
					*(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])(i,0);
		x[currentPointId] = currentX;

		double currentY = 0;
		for(int i=0;i<3;i++)
			currentY += connectionVectorBetweenPointAndReferencePoint[i]
				 // *(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])[i][1];
					*(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])(i,1);
		y[currentPointId] = currentY;

        if(currentX != 0)
            phi[currentPointId] = (double) (currentX <  0 ) * (atan(currentY/currentX)+(-(currentY<0)+(currentY>0)+(currentY==0)*(currentX<0))*atan(1.0)*4)	//	atan(1.0)*4 = pi = 3.1415...
                              + (currentX >  0 ) * atan(currentY/currentX);
        else
        {
            if(currentY != 0)
                phi[currentPointId] = (double) (currentX == 0 ) * ( ( ( currentY > 0 ) - ( currentY < 0 ) ) *atan(1.0)*2 );
            else
                phi[currentPointId] = 0 ;
                //std::cerr << "There exist integration points that lie ON the cut front" << std::endl;
        }
	}
}


void CutModelTopology_CGAL::getRAndPhiOfPointDependentOnCutFront(Point_t &point,double &r, double &phi)
{
	// we could also choose to check the distance to the tree that derives from the cut direction - maybe we even have
	// to build a triangle tree of all the triangles involved in the element
	TreeSegment_t CutBoundaryAABBTree(m_cutBoundary.begin(),m_cutBoundary.end());
	TreeSegment_Point_and_primitive_id_t closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary = CutBoundaryAABBTree.closest_point_and_primitive(point);
	Point_t closestPoint = closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary.first;

	Vector_t connectionVectorBetweenPointAndReferencePoint = point-closestPoint;
	double squaredLengthOfVector = connectionVectorBetweenPointAndReferencePoint.squared_length();
	r = sqrt(squaredLengthOfVector);

	TreeSegment_Primitive_id_t PrimitiveIdOfClosestSegmentOfCutBoundary = closestPointToCutBoundaryANDidOfClosestSegmentOfCutBoundary.second;
	Segments_t::iterator IteratorOfClosestSegmentOfCutBoundary = PrimitiveIdOfClosestSegmentOfCutBoundary;
	unsigned int IdOfClosestSegmentOfCutBoundary = std::distance(m_cutBoundary.begin(),IteratorOfClosestSegmentOfCutBoundary);

	Matrix_3_3_t currentTransformationMatrix = m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary];
	double x = 0;
	for(int i=0;i<3;i++)
		// x += connectionVectorBetweenPointAndReferencePoint[i]
		// 	 *(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])[i][0];
        x += connectionVectorBetweenPointAndReferencePoint[i]*currentTransformationMatrix(i,0);
		
	double y = 0;
	for(int i=0;i<3;i++)
		// y += connectionVectorBetweenPointAndReferencePoint[i]
		// 	 *(m_cutBoundaryTransformationMatrices[IdOfClosestSegmentOfCutBoundary])[i][1];
        y += connectionVectorBetweenPointAndReferencePoint[i]*currentTransformationMatrix(i,1);

    if(x != 0)
        phi = (double) (x <  0 ) * (atan(y/x)+(-(y<0)+(y>0)+(y==0)*(x<0))*atan(1.0)*4)	//	atan(1.0)*4 = pi = 3.1415...
                          + (x >  0 ) * atan(y/x);
    else
    {
        if(y != 0)
            phi = (double) (x == 0 ) * ( ( ( y > 0 ) - ( y < 0 ) ) *atan(1.0)*2 );
        else
            phi = 0 ;
            //std::cerr << "There exist integration points that lie ON the cut front" << std::endl;
    }
}


void CutModelTopology_CGAL::calculateNormalOfTriangle(Triangle_t triangle, Vector_t &normalOfTriangle)
{
	Vector_t vectorPoint0ToPoint1(triangle.vertex(0),triangle.vertex(1));
	Vector_t vectorPoint1ToPoint2(triangle.vertex(1),triangle.vertex(2));

	normalOfTriangle = cross_product(vectorPoint0ToPoint1,vectorPoint1ToPoint2);
}


void CutModelTopology_CGAL::secureThatTriangleNormalsInSameDirection(Triangle_t triangle1, Triangle_t &triangle2) // the second triangle will be changed if the normal doesn't point in the right direction
{
	Vector_t normalOfTriangle1;
	calculateNormalOfTriangle(triangle1,normalOfTriangle1);
	Vector_t normalOfTriangle2;
	calculateNormalOfTriangle(triangle2,normalOfTriangle2);

	if(normalOfTriangle1*normalOfTriangle2<0)
		triangle2 = Triangle_t(triangle2.vertex(0),triangle2.vertex(2),triangle2.vertex(1));
}


void CutModelTopology_CGAL::secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
(Triangle_t triangle1,Points_t &points, ConnectTriangle_t currentConnectTriangle)
{
	Vector_t normalOfTriangle1;
	calculateNormalOfTriangle(triangle1,normalOfTriangle1);

	Triangle_t triangle2(	points[currentConnectTriangle.get<0>()],
						points[currentConnectTriangle.get<1>()],
						points[currentConnectTriangle.get<2>()]);
	Vector_t normalOfTriangle2;
	calculateNormalOfTriangle(triangle2,normalOfTriangle2);

	if(normalOfTriangle1*normalOfTriangle2<0)
	{
		currentConnectTriangle =
			ConnectTriangle_t(currentConnectTriangle.get<0>(),currentConnectTriangle.get<2>(),currentConnectTriangle.get<1>());
	}
}


void CutModelTopology_CGAL::getInformationAboutPartialAndCompleteCut(Segments_t &segments, Tree_t &CutAABBTree, int edgesOnTetrahedronFaces[][3],
			int ObjectSegmentIsSplit[][6], Ints_t &ObjectSplitSegmentIds, int ObjectTriangleIsSeparated[][4], int ObjectTrianglePartiallySeparated[][4],
			Ints_t &ObjectTriangleIdsThatIntersectWithCut,
			int* ObjectTetrahedronIsSeparated, int* ObjectTetrahedronPartiallySeparated,
			Ints_t &ObjectTetrahedraIdsThatIntersectWithCut)
{
	for(Segments_t::iterator segmentIterator=segments.begin(); segmentIterator!=segments.end();++segmentIterator)
	{
		Segment_t currentSegment = *segmentIterator;
		int globalEdgeId = std::distance(segments.begin(),segmentIterator);
		int elementId = globalEdgeId/6;
		ObjectSegmentIsSplit[elementId][globalEdgeId-elementId*6] = CutAABBTree.number_of_intersected_primitives(currentSegment)>0;
		if(CutAABBTree.number_of_intersected_primitives(currentSegment)>1)
        {
            std::cout << "ObjectSegmentId " << globalEdgeId << " is intersected " << CutAABBTree.number_of_intersected_primitives(currentSegment) << " times, the programm can't handle that!!"<< endl;
            SHOW(currentSegment.point(0));
            SHOW(currentSegment.point(1));
        }
		if(ObjectSegmentIsSplit[elementId][globalEdgeId-elementId*6])
			ObjectSplitSegmentIds.push_back(globalEdgeId);
		if(!((globalEdgeId-5)%6)) // every time when we checked six segments of an element, we check whether the object element and its faces are intersecting with the cut
		{
			int helpEdge=0;
			for(int edgeId=0;edgeId<6;edgeId++)
			{
				helpEdge += ObjectSegmentIsSplit[elementId][edgeId];
			}
			if(helpEdge)
			{
				int helpElement=0;
				for(int faceId=0;faceId<4;faceId++)
				{
					int helpFace=0;
					for(int edgeIdFace=0;edgeIdFace<3;edgeIdFace++)
					{
						helpFace += ObjectSegmentIsSplit[elementId][edgesOnTetrahedronFaces[faceId][edgeIdFace]];
					}
					ObjectTriangleIsSeparated[elementId][faceId] = helpFace>1;
					if(helpFace>0)
					{
						ObjectTriangleIdsThatIntersectWithCut.push_back(elementId*4+faceId);
					}
					ObjectTrianglePartiallySeparated[elementId][faceId] = helpFace==1;
					helpElement += ObjectTriangleIsSeparated[elementId][faceId];
				}
				ObjectTetrahedronIsSeparated[elementId] = helpElement>2;
				ObjectTetrahedronPartiallySeparated[elementId] = helpElement<3;
				ObjectTetrahedraIdsThatIntersectWithCut.push_back(elementId);
			}
			else
			{
				ObjectTetrahedronIsSeparated[elementId] = 0;
				ObjectTetrahedronPartiallySeparated[elementId] = 0;
				for(int faceId=0;faceId<4;faceId++)
				{
					ObjectTriangleIsSeparated[elementId][faceId] = 0;
					ObjectTrianglePartiallySeparated[elementId][faceId] = 0;
				}
			}
		}
	}
}


void CutModelTopology_CGAL::getVectorIdOfIntEntry(int Entry, Ints_t Vector, int &Id)
{
	int maxId = std::distance(Vector.begin(),Vector.end());
	Id=-1; // this is for the case, that the entry is not in the Vector -> when we plug that Id into the vector, we receive an error and can see what the reason is!
	for(int currentId=0;currentId<maxId;currentId++)
	{
		if(Entry==Vector[currentId])
		{
			Id=currentId;
		}
	}
}


void CutModelTopology_CGAL::getVectorIdsOfIntEntry(int Entry, Ints_t Vector, Ints_t &Ids)
{
	int maxId=std::distance(Vector.begin(),Vector.end());
	for(int currentId=0;currentId<maxId;currentId++)
	{
		if(Entry==Vector[currentId])
		{
			Ids.push_back(currentId);
		}
	}
}


void CutModelTopology_CGAL::getVectorIdsOf3IntEntriesIn3Vectors(int Entry0, int Entry1, int Entry2,
		Ints_t Vector0, Ints_t Vector1, Ints_t Vector2, Ints_t &Ids)
{
	int maxId=std::distance(Vector0.begin(),Vector0.end());
	for(int currentId=0;currentId<maxId;currentId++)
	{
		if(	(Entry0==Vector0[currentId]	|| Entry0==-2)&&
			(Entry1==Vector1[currentId]	|| Entry1==-2)&&
			(Entry2==Vector2[currentId]	|| Entry2==-2))
		{
			Ids.push_back(currentId);
		}
	}
}


void CutModelTopology_CGAL::getVectorIdsOf6IntEntriesIn6Vectors(int Entry0, int Entry1, int Entry2, int Entry3, int Entry4, int Entry5,
		Ints_t &Vector0, Ints_t &Vector1, Ints_t &Vector2, Ints_t &Vector3, Ints_t &Vector4, Ints_t &Vector5, Ints_t &Ids)
{
	int maxId=std::distance(Vector0.begin(),Vector0.end());
	for(int currentId=0;currentId<maxId;currentId++)
	{
		if(	(Entry0==Vector0[currentId]	|| Entry0==-2)&&
			(Entry1==Vector1[currentId]	|| Entry1==-2)&&
			(Entry2==Vector2[currentId]	|| Entry2==-2)&&
			(Entry3==Vector3[currentId]	|| Entry3==-2)&&
			(Entry4==Vector4[currentId]	|| Entry4==-2)&&
			(Entry5==Vector5[currentId]	|| Entry5==-2))
		{
			Ids.push_back(currentId);
		}
	}
}


void CutModelTopology_CGAL::deleteMultiplePointsAndAdaptConnectTriangle(int numberOfTriangles, Points_t &points, int connect[][3])
{
	int numberOfPoints = std::distance(points.begin(),points.end());
	if(numberOfPoints)
	{
		int indices[numberOfPoints];
		deleteMultiplePoints(points,indices);
		adaptConnectTriangle(numberOfTriangles,indices,connect);
	}
}


void CutModelTopology_CGAL::adaptConnectTriangle(int numberOfTriangles, int* indices, int connect[][3])
{
	for(int triangleId=0;triangleId<numberOfTriangles;triangleId++)
	{
		for(int pointId=0;pointId<3;pointId++)
		{
			connect[triangleId][pointId] = indices[connect[triangleId][pointId]];
		}
	}
}


void CutModelTopology_CGAL::deleteMultiplePoints(Points_t &points, int *indices)
{
	int numberOfPoints = std::distance(points.begin(),points.end());
	int indicesAfterSort[numberOfPoints];
	sortPoints(points,indicesAfterSort);

	Points_t temporaryPoints;
	int indicesAfterDelete[numberOfPoints];
	int newPointId=0;
	temporaryPoints.push_back(points[0]);
	indicesAfterDelete[0]=0;
	for(int pointId=1;pointId<numberOfPoints;pointId++)
	{
		if(points[pointId-1]!=points[pointId])
		{
			temporaryPoints.push_back(points[pointId]);
			newPointId++;
		}
		indicesAfterDelete[pointId] = newPointId;
	}

	points = temporaryPoints;
	for(int pointId=0;pointId<numberOfPoints;pointId++)
		indices[pointId] = indicesAfterDelete[indicesAfterSort[pointId]];
}


void CutModelTopology_CGAL::deleteMultipleInts(Ints_t &Ids)
{
	std::sort(Ids.begin(),Ids.end());

	Ints_t temporaryIds;
	int lastId = *Ids.begin();
	temporaryIds.push_back(lastId);
	for(Ints_t::iterator iter=Ids.begin();iter!=Ids.end();iter++)
	{
		int currentId = *iter;
		if(currentId!=lastId)
		{
			temporaryIds.push_back(currentId);
		}
		lastId = *iter;
	}

	Ids = temporaryIds;
}


void CutModelTopology_CGAL::sortPoints(Points_t &points, int* indices)	// both inputs are outputs as well and will be overwritten
{
	std::vector<PointAndInt_t> PointsAndIndices;
	int numberOfPoints = std::distance(points.begin(),points.end());
	for(int pointId=0;pointId<numberOfPoints;pointId++)
	{
		PointAndInt_t currentPointAndIndice;
		currentPointAndIndice.first = points[pointId];
		currentPointAndIndice.second = pointId;

		PointsAndIndices.push_back(currentPointAndIndice);
	}

	std::sort(PointsAndIndices.begin(),PointsAndIndices.end(), this->compare_PointAndInt);

	Points_t temporaryPoints;
	for(int pointId=0;pointId<numberOfPoints;pointId++)
	{
		temporaryPoints.push_back(PointsAndIndices[pointId].first);
		indices[PointsAndIndices[pointId].second]=pointId;
	}

	points = temporaryPoints;
}


bool CutModelTopology_CGAL::compare_PointAndInt(PointAndInt_t A, PointAndInt_t B)
{
	return compare_Point(A.first,B.first);
}


bool CutModelTopology_CGAL::compare_Point(Point_t P,Point_t Q)
{
	return P[0]<Q[0] + (P[0]==Q[0])*(P[1]<Q[1]) + (P[0]==Q[0])*(P[1]==Q[1])*(P[2]<Q[2]);
}


bool CutModelTopology_CGAL::compare_DoubleAndInt( const DoubleAndInt_t& l, const DoubleAndInt_t& r)
{
	return l.first < r.first;
}


void CutModelTopology_CGAL::outputFromVTKContainersToVTKFile(vtkPoints* pointContainerTopoDebug, vtkCellArray* cellContainerTopoDebug, int numberOfCellNodes, const char* filenameTopoDebug)
{
	vtkUnstructuredGrid* myGridTopoDebug = vtkUnstructuredGrid::New();
	myGridTopoDebug->SetPoints(pointContainerTopoDebug);
	if(numberOfCellNodes==3)
	{
		myGridTopoDebug->SetCells(VTK_TRIANGLE, cellContainerTopoDebug);
	}
	else if (numberOfCellNodes==2)
	{
		myGridTopoDebug->SetCells(VTK_LINE, cellContainerTopoDebug);
	}

	vtkUnstructuredGridWriter* writerTopoDebug = vtkUnstructuredGridWriter::New();

	writerTopoDebug->SetInputData(myGridTopoDebug);
	writerTopoDebug->SetFileName(filenameTopoDebug);
	writerTopoDebug->Write();

	pointContainerTopoDebug->Delete();
	cellContainerTopoDebug->Delete();
	myGridTopoDebug->Delete();
	writerTopoDebug->Delete();
}


void CutModelTopology_CGAL::outputFromTrianglesToVTKFile(Triangles_t &triangles, const char* filename)
{
	int numberOfCellNodes = 3;
	int numberOfCells = std::distance(triangles.begin(),triangles.end());
	Points_t points;
	TrianglesToPoints(triangles, points);
	vtkPoints* pointContainerVTK = vtkPoints::New();
	PointsToVTKPointContainer(points, pointContainerVTK);
	vtkCellArray* triangleContainerVTK = vtkCellArray::New();
	CellsToVTKCellContainer(numberOfCells, numberOfCellNodes, triangleContainerVTK);

	outputFromVTKContainersToVTKFile(pointContainerVTK, triangleContainerVTK, numberOfCellNodes, filename);
}


void CutModelTopology_CGAL::outputFromSegmentsToVTKFile(Segments_t &segments, const char* filename)
{
	int numberOfCellNodes = 2;
	int numberOfCells = std::distance(segments.begin(),segments.end());
	Points_t points;
	SegmentsToPoints(segments, points);
	vtkPoints* pointContainerVTK = vtkPoints::New();
	PointsToVTKPointContainer(points, pointContainerVTK);
	vtkCellArray* lineContainerVTK = vtkCellArray::New();
	CellsToVTKCellContainer(numberOfCells, numberOfCellNodes, lineContainerVTK);

	outputFromVTKContainersToVTKFile(pointContainerVTK, lineContainerVTK, numberOfCellNodes, filename);
}


void CutModelTopology_CGAL::outputFromPointsToVTKFile(Points_t &points, const char* filenameTopoDebug)
{
	vtkPoints* pointContainerVTK = vtkPoints::New();
	PointsToVTKPointContainer(points, pointContainerVTK);

	vtkPolyData* myGridTopoDebug = vtkPolyData::New();
	myGridTopoDebug->SetPoints(pointContainerVTK);

	vtkPolyDataWriter* writerTopoDebug = vtkPolyDataWriter::New();

	writerTopoDebug->SetInputData(myGridTopoDebug);
	writerTopoDebug->SetFileName(filenameTopoDebug);
	writerTopoDebug->Write();

	pointContainerVTK->Delete();
	myGridTopoDebug->Delete();
	writerTopoDebug->Delete();
}


void CutModelTopology_CGAL::TrianglesToPoints(Triangles_t &triangles, Points_t &points)
{
	for(Triangles_t::iterator triangleIterator=triangles.begin(); triangleIterator!=triangles.end();++triangleIterator)
	{
		Triangle_t currentTriangle = *triangleIterator;
		for(int pointId=0; pointId<3;pointId++)
		{
			Point_t trianglePoint = currentTriangle.vertex(pointId);
			points.push_back(trianglePoint);
		}
	}
}


void CutModelTopology_CGAL::SegmentsToPoints(Segments_t &segments, Points_t &points)
{
	for(Segments_t::iterator segmentIterator=segments.begin(); segmentIterator!=segments.end();++segmentIterator)
	{
		Segment_t currentSegment = *segmentIterator;
		int segmentId = std::distance(segments.begin(),segmentIterator);
		for(int pointId=0; pointId<2;pointId++)
		{
			Point_t segmentPoint = currentSegment.vertex(pointId);
			points.push_back(segmentPoint);
		}
	}
}


void CutModelTopology_CGAL::PointsToVTKPointContainer(Points_t &points, vtkPoints* pointContainerVTK)
{
	int numberOfPoints = std::distance(points.begin(),points.end());
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


void CutModelTopology_CGAL::CellsToVTKCellContainer(int numberOfCells, int numberOfCellNodes, vtkCellArray* cellContainerVTK)
{
	vtkIdType* currentCell;
	currentCell = new vtkIdType[numberOfCellNodes];
	for(int cellId=0; cellId<numberOfCells; cellId++)
	{
		for(int nodeId=0; nodeId<numberOfCellNodes; nodeId++)
		{
			currentCell[nodeId] = cellId*numberOfCellNodes+nodeId;
		}

		cellContainerVTK->InsertNextCell(numberOfCellNodes, currentCell);
	}
}
