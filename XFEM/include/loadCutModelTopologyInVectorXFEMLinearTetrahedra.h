/* loadCutModelTopologyInVectorXFEMLinearTetrahedra.h
 *
 * 		Created On: 19.07.2013
 * 		Author: 	Christoph Joachim Paulus
 *
 * */



//defines
#ifndef loadCutModelTopologyInVectorXFEMLinearTetrahedra_H_
#define loadCutModelTopologyInVectorXFEMLinearTetrahedra_H_

// includes
#include "CutModelTopology_CGAL.h"
#include "XFemLinearTetrahedron.h"

/*
class loadCutModelTopologyInVectorXFEMLinearTetrahedra
- initialization of the cut model topology
- transfer of the obtained information to a vector of xfem linear tetrahedra
*/

class loadCutModelTopologyInVectorXFEMLinearTetrahedra
{
public:
	//typedefs
	typedef std::vector<int> Ints_t;
	typedef std::vector<double> Doubles_t;
	// functions
    loadCutModelTopologyInVectorXFEMLinearTetrahedra();
	virtual ~loadCutModelTopologyInVectorXFEMLinearTetrahedra();

	void init(std::string ObjectFilename, std::string CutFilename,	std::string RefGridFilename, 
			double lambda,double mu,int numberOfIntegrationPoints,
			Ints_t &pointsReferenceSolutionInObjectTetrahedronId, Ints_t &pointsReferenceSolutionAbove,Doubles_t &pointsReferenceSolution_r, Doubles_t &pointsReferenceSolution_phi, 
	   std::vector<XFemLinearTetrahedron> &VectorXFemLinearTetrahedra);	// output

	void informationAboutPointsOnTetrahedron(int tetrahedronId,
			Ints_t &pointIsAbove, Ints_t &signEnrichedVertexIds, Ints_t &branchEnrichedVertexIds, Ints_t &pointIds, XFemLinearTetrahedron::Points_t &points);
	void getConnectTriangle(int tetrahedronId,Ints_t &pointIds,Ints_t &pointIsAbove,XFemLinearTetrahedron::ConnectTriangles_t &connectTriangle, Ints_t &triangleInFaceId, Ints_t &triangleAbove);
	bool faceIntersectsWithCut(int globalFaceId);
	int getNumberOfSignEnrichedPoints();
	int getNumberOfBranchEnrichedPoints();
	int getNumberOfPoints();
//	void outputFromVTKContainersToVTKFile(vtkPoints* pointContainer, vtkCellArray* cellContainer, int numberOfCellNodes, const char* filename);

//private:
	CutModelTopology_CGAL m_cutModel;
};

#endif /* loadCutModelTopologyInVectorXFEMLinearTetrahedra_H_ */
