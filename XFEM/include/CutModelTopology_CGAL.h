/* CutModelTopology_CGAL.h
 *
 * 		Created On: 14.07.2013
 * 		Author: 	Christoph Joachim Paulus
 *
 * */



//defines
#ifndef CutModelTopology_CGAL_H_
#define CutModelTopology_CGAL_H_

// includes

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/squared_distance_3.h>

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkType.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"

#include <math.h>
#include <time.h>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <CGAL/Simple_cartesian.h>

#include <Eigen/Core>

#include "utils.h"

/*
Class CutModelTopology_CGAL:
- loading of the meshs of the object (tetrahedra) and the cut (triangles)
- representation of those meshes in another form: segments and triangles
- check whether the nodes of the object mesh are above or below the triangulated cut
- check how the elements are cut (not, partially or completely cut)
- receive the kind of the enrichment of the nodes (not, sign or branch enriched)
- adjust the connects of the elements with enriched nodes
- calculate the interesections of the object tetrahedral mesh with the cut triangle mesh:
	- receive the intersections of cut segments with object triangles
	- receive the intersections of object segments with cut triangles
	- identify the points of the mesh of the cut that are in the object
- obtain the part of the cut boundary that are in the cut elements (these are necessary for the calculation of the asymptotic crack tip functions)
- triangulation of the new surfaces:
	- the triangles of the cut that lie completely inside an element
	- the faces of the old tetrahedra

- this class also contains a method that calculates the positions of a mesh with another resolution in our current object mesh (one of the public functions)
*/

class CutModelTopology_CGAL
{
public:
	// typedefs
	// CGAL typedefs
	typedef CGAL::Simple_cartesian<double> K_t;
	typedef K_t::Point_3 Point_t;
	typedef std::vector<Point_t> Points_t;
	typedef K_t::Vector_3 Vector_t;
	typedef K_t::Segment_3 Segment_t;
	typedef std::vector<Segment_t> Segments_t;
	typedef K_t::Triangle_3 Triangle_t;
	typedef std::vector<Triangle_t> Triangles_t;
	typedef std::vector<Triangle_t>::iterator Iterator_t;
	typedef CGAL::AABB_triangle_primitive<K_t,Iterator_t> Primitive_t;
	typedef CGAL::AABB_traits<K_t, Primitive_t> Traits_t;
	typedef CGAL::AABB_tree<Traits_t> Tree_t;
	typedef Tree_t::Object_and_primitive_id Object_and_primitive_id_t;
	typedef std::vector<Object_and_primitive_id_t> Intersections_t;
	typedef Tree_t::Point_and_primitive_id Point_and_primitive_id_t;
	typedef Tree_t::Primitive_id Primitive_id_t;

	typedef std::vector<Segment_t>::iterator IteratorSegment_t;
	typedef CGAL::AABB_segment_primitive<K_t,IteratorSegment_t> PrimitiveSegment_t;
	typedef CGAL::AABB_traits<K_t, PrimitiveSegment_t> TraitsSegment_t;
	typedef CGAL::AABB_tree<TraitsSegment_t> TreeSegment_t;
	typedef TreeSegment_t::Point_and_primitive_id TreeSegment_Point_and_primitive_id_t;
	typedef TreeSegment_t::Primitive_id TreeSegment_Primitive_id_t;

	// general typedefs
	typedef std::vector<int> Ints_t;
	typedef std::vector<double> Doubles_t;
	typedef std::pair<double,int> DoubleAndInt_t;
	typedef std::pair<Point_t,int> PointAndInt_t;

	typedef boost::tuple<int,int,int> ConnectTriangle_t;
	typedef std::vector<ConnectTriangle_t> ConnectTriangles_t;
	typedef boost::tuple<int,int,int,int> ConnectTetrahedron_t;
	typedef std::vector<ConnectTetrahedron_t> ConnectTetrahedra_t;

	typedef Eigen::Matrix<double,3,3> Matrix_3_3_t;


	// functions
	CutModelTopology_CGAL(void);
	virtual ~CutModelTopology_CGAL();

	// initialization of the cut model, this function will be executed at the beginning of the simulation to load the meshes of the object and the cut in order to obtain the information for the enriched elements
	void init(std::string ObjectFilename, std::string CutFilename);
    void loadObject();
	void getVectorIdOfIntEntry(int Entry, Ints_t Vector, int &Id);
	void getVectorIdsOfIntEntry(int Entry, Ints_t Vector, Ints_t &Ids);
	void pointsAbove(Points_t &points,Ints_t &pointsAbove);	// which of the points are above/below the cut?
    void pointsInObjectTetrahedronNumber(std::string RefGridFilename, Ints_t &pointsReferenceSolutionInObjectTetrahedronId, Ints_t &pointsReferenceSolutionAbove, Doubles_t &pointsReferenceSolution_r, Doubles_t &pointsReferenceSolution_phi);	// in order to check the convergence, this function need to be called
	void getR_Phi_XAndYOfPointsDependentOnCutFront(Points_t &points,
			Doubles_t &r, Doubles_t &phi, Doubles_t &x, Doubles_t &y,Ints_t &idClosestCutBoundarySegment);
	void getRAndPhiOfPointDependentOnCutFront(Point_t &point,double &r, double &phi);
	double getROfPoint(Point_t point);
	int pointAbove(Point_t &point);

private:
	void getTetConnectFromVTKCellArray(vtkCellArray* myArray, int numberOfCellNodes, int objectConnectTetrahedra[][4], vtkIdType& numberOfCells, int &numberOfObjectNodes);
	void getTriangleConnectFromVTKCellArray(vtkCellArray* myArray, int numberOfCellNodes, int objecttriangleConnect[][3], vtkIdType& numberOfCells, int &numberOfObjectNodes);
	void getPointsFromVTKCellArray(int numberOfPoints, vtkPoints* myPoints, double factor, double pointsDouble[][3], Points_t& points);

	void getConnectTriangleOfConnectTet(int numberOfCells, int connectTetrahedra[][4], int nodesOnTetrahedronFaces[][3], int ConnectTriangle[][3]);
	void getConnectSegmentOfConnectTet(int numberOfCells, int connectTetrahedra[][4], int nodesOnTetrahedronEdges[][2], int ConnectSegment[][2]);
	void getConnectSegmentOfConnectTriangle(int numberOfCutTriangles, int ConnectTriangle[][3], int ConnectSegment[][2]);

	void getTrianglesFromConnectTriangleAndPointsDouble(int numberOfCells, int ConnectTriangle[][3], double pointsDouble[][3], std::vector<Triangle_t> &triangles);
	void getSegmentsFromConnectSegmentAndPointsDouble(int numberOfCells, int ConnectSegment[][2], double pointsDouble[][3], Segments_t &segments);

	void getTetrahedronThatSurroundsPoint(Triangles_t &objectTriangles,Point_t point,int &pointIsInObjectTetrahedronId);
	void isPointAbove(Points_t &points, Tree_t &CutAABBTree, Triangles_t &CutTriangles, int* pointsAbove, Ints_t pointIdsAbove);
	void getTransformationMatrixForOneElementOfCutboundary
		(int* cutBoundarySegmentNodeIdsOfCutBoundaryTriangle,Triangle_t &cutBoundaryTriangle,Matrix_3_3_t &transformationMatrix);

	void calculateNormalOfTriangle(Triangle_t triangle, Vector_t &normalOfTriangle);
	void secureThatTriangleNormalsInSameDirection(Triangle_t triangle1, Triangle_t &triangle2);
	void secureThatTriangleNormalOfConnectTriangleTypeGoesInSameDirectionLikeTheNormalOfTriangle
	(Triangle_t triangle1,Points_t &points, ConnectTriangle_t currentConnectTriangle);


	void getInformationAboutPartialAndCompleteCut(Segments_t &segments, Tree_t &CutAABBTree, int edgesOnTetrahedronFaces[][3],
				int objectSegmentIsSplit[][6], Ints_t &objectSplitSegmentIds, int objectTriangleIsSeparated[][4], int objectTrianglePartiallySeparated[][4],
				Ints_t &objectTriangleIdsThatIntersectWithCut,
				int* objectTetrahedronIsSeparated, int* objectTetrahedronPartiallySeparated,
				Ints_t &objectTetrahedraIdsThatIntersectWithCut);

	void getVectorIdsOf3IntEntriesIn3Vectors(int Entry0, int Entry1, int Entry2,
			Ints_t Vector0, Ints_t Vector1, Ints_t Vector2, Ints_t &Ids);
	void getVectorIdsOf6IntEntriesIn6Vectors(int Entry0, int Entry1, int Entry2, int Entry3, int Entry4, int Entry5,
			Ints_t &Vector0, Ints_t &Vector1, Ints_t &Vector2, Ints_t &Vector3, Ints_t &Vector4, Ints_t &Vector5, Ints_t &Ids);
//	void getVectorIdsOf5IntEntriesIn5Vectors(int Entry0, int Entry1, int Entry2, int Entry3, int Entry4,
//			Ints_t Vector0, Ints_t Vector1, Ints_t Vector2, Ints_t Vector3, Ints_t Vector4, Ints_t &Ids);

	void deleteMultiplePointsAndAdaptConnectTriangle(int numberOfTriangles, Points_t &points, int connect[][3]);
	void adaptConnectTriangle(int numberOfTriangles, int* indices, int connect[][3]);
	void deleteMultiplePoints(Points_t &points, int *indices);
	void deleteMultipleInts(Ints_t &Ids);
	void sortPoints(Points_t &points, int* indices);

	static bool compare_PointAndInt(PointAndInt_t A, PointAndInt_t B);
	static bool compare_Point(Point_t P,Point_t Q);
	static bool compare_DoubleAndInt( const DoubleAndInt_t& l, const DoubleAndInt_t& r);

	void outputFromVTKContainersToVTKFile(vtkPoints* pointContainerTopoDebug, vtkCellArray* cellContainerTopoDebug, int numberOfCellNodes, const char* filenameTopoDebug);
	void outputFromTrianglesToVTKFile(Triangles_t &triangles, const char* filename);
	void outputFromSegmentsToVTKFile(Segments_t &segments, const char* filename);
	void outputFromPointsToVTKFile(Points_t &points, const char* filenameTopoDebug);

	void TrianglesToPoints(Triangles_t &triangles, Points_t &points);
	void SegmentsToPoints(Segments_t &segments, Points_t &points);
	void PointsToVTKPointContainer(Points_t &points, vtkPoints* pointContainerVTK);
	void CellsToVTKCellContainer(int numberOfCells, int numberOfCellNodes, vtkCellArray* cellContainerVTK);

public:
	// declaration of class variables
	Points_t m_objectPoints;
	ConnectTetrahedra_t m_objectConnectTetrahedra;
	std::vector<Ints_t>	m_objectConnectTetrahedraWithEnrichedNodes;
//	std::vector<Ints_t> m_nodeIdInTetrahedraIds;
	Ints_t m_globalFaceIdHasNeighbourWithglobalFaceId;
	Ints_t m_pointIdsAbove;
	Triangles_t m_objectTriangles;

	Triangles_t m_cutTriangles;
	Segments_t m_cutBoundary;
	Triangles_t m_cutBoundaryTriangles;
	Ints_t m_cutBoundaryInObjectTetrahedronId;
	std::vector<Matrix_3_3_t> m_cutBoundaryTransformationMatrices;

	std::vector<int> m_objectTetrahedraIdsThatIntersectWithCut;
	std::vector<int> m_objectTriangleIdsThatIntersectWithCut;
	std::vector<int> m_objectSplitSegmentIds;

	Ints_t m_signEnrichedObjectPoints;
	int m_numberOfSignEnrichedNodes;
	Ints_t m_branchEnrichedObjectPoints;
	int m_numberOfBranchEnrichedNodes;

	Points_t m_intersectionPoints;
	Ints_t m_intersectionPointInObjectTetrahedronId;
	Ints_t m_intersectionPointInObjectTriangleId;
	Ints_t m_intersectionPointInObjectSegmentId;
	Ints_t m_intersectionPointInCutTriangleId;
	Ints_t m_intersectionPointInCutSegmentId;
	Ints_t m_intersectionPointAbove;

	ConnectTriangles_t m_triangulationOfCutObjectTetrahedronFacesConnectTriangle;
	Ints_t m_TriangleOfTriangulationAbove;
	Ints_t m_TriangleOfTriangulationInObjectTriangleId;
    Ints_t m_TriangleOfTriangulationInObjectTetrahedronId;
};

#endif /* CutModelTopology_CGAL_H_ */
