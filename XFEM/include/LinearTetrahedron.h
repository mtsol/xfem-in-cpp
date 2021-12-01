/*
 * linearTetrahedron.h
 *
 *  Created on: March 28, 2013
 *      Author: Christoph Paulus
 */

//defines
#ifndef LINEARTETRAHEDRON_H_
#define LINEARTETRAHEDRON_H_

#define NUMBER_OF_ELEMENT_NODES 4
#define NUMBER_OF_INTEGRATION_POINTS_TET 1

//#ifndef NUMBER_OF_ELEMENT_NODES	// FIXME at the start of Tetrahedron.h!
//#define NUMBER_OF_ELEMENT_NODES 4

// includes
#include "FEM_Elementtype.h"
#include <math.h>	// for floor function

#include <Eigen/Core>
#include <Eigen/LU>

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkType.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"

/*
class LinearTetrahedron : provides all functions and variables that are necessary to handle a linear tetrahedron
- shape functions
- derivatives of shape functions
- the transformation between the tetrahedron and the reference tetrahedron
- the derivative of the transformation mentioned above
- integration points (equally distributed)
- calculation of the stiffness matrix
- a public function that pushes the displaced DOFs into the structure of LinearTetrahedron
- the displaced points, based on the displacements of the degrees of freedom and the shape functions
- methods that can output the meshes to a vtk file
*/

class LinearTetrahedron : public FEM_Elementtype
{
public:
	//typedefs
	typedef std::vector<int> Ints_t;
	typedef std::vector<double> Doubles_t;

	typedef Eigen::Matrix<double,DIM,1> Point_t;
	typedef std::vector<Point_t> Points_t;
	typedef Eigen::Matrix<int,3,1> ConnectTriangle_t;
	typedef std::vector<ConnectTriangle_t> ConnectTriangles_t;

	typedef Eigen::Matrix<double,3,3> Matrix_3_3_t;
	typedef Doubles_t ShapeFct_t;
	typedef std::vector<ShapeFct_t> ShapeFcts_t;
	typedef Points_t ShapeFctDeriv_t;	// could also be defined as std::vector<ShapeFct_t>
	typedef Eigen::MatrixXd StiffnessMatrix_t;
		// if there is a dynamic matrix, then we could avoid having so big matrices for uncut elements



	// functions
	LinearTetrahedron(void);
	virtual ~LinearTetrahedron();
	virtual void init(int tetrahedronId,Ints_t connect,int neighboursOfFaces[][2],Points_t &points,double lambda,double mu,int numberOfIntegrationPoints);
	virtual void getGlobalToLocal();
	virtual void getShapeFunctionsOf_m_points();

	virtual void getShapeFunctionDependentOnLocalCoords(Point_t xi, ShapeFct_t &ShapeFunction);
	virtual void getShapeFunction(Point_t x, ShapeFct_t &ShapeFunction);
	virtual void transformLocalToGlobalCoordinatesByUsingLinearShapeFunctions(Points_t &TetCoords, Point_t xi, Point_t &x);
	virtual void transformGlobalToLocalCoordinatesByUsingLinearShapeFunctions(Point_t x, Point_t &xi);
	virtual void getShapeFunctionLocalDerivatives(Point_t xi,
			ShapeFctDeriv_t &ShapeFunctionLocalDerivatives);

	virtual void getJacobian(ShapeFctDeriv_t &ShapeFunctionLocalDerivatives, Matrix_3_3_t &jacobian);
    virtual double getJacobianDeterminant(Point_t xi = Point_t());
	virtual void getShapeFunctionDerivatives(Point_t xi, ShapeFctDeriv_t &ShapeFunctionDerivative);
	virtual void getShapeFunctionDerivativesDependentOnIntegrationPointNumber(int integrationPointNumber, ShapeFctDeriv_t &ShapeFunctionDerivative);
	virtual void getNumberOfEquallyDistributedIntegrationPoints(int numberOfIntegrationPointsWanted);
	virtual void getGaussIntegrationPoints(Points_t &integrationPoints, Doubles_t &integrationWeights);
	virtual void getGlobalIntegrationPointsIn_Cube_CubeWithoutTet_Tet(
			int walkX,int walkY,int walkZ,
			Eigen::Matrix<double,8,3> &CoordsReferenceCube, int tetrahedronConnectOfCube8[][4],
			int numberOfTetrahedra);
	virtual void getEquallyDistributedIntegrationPoints();
	//virtual void addTetrahedronVertexesToIntegrationPoints();
	virtual void stiffnessBlock();
	virtual void stiffnessBlockDependentOnShapeFunctionDerivatives(ShapeFctDeriv_t &SFDs,StiffnessMatrix_t &stiffnessBlock);
	virtual void stiffnessBlock2();
	virtual void stiffnessBlock3();
    virtual void writeDisplacementInTetrahedronData(Points_t &displacement);
//    virtual void displacePoint(ShapeFct_t& SF, Point_t &x_displaced);
    virtual Point_t displacePoint(ShapeFct_t& SF, Point_t& initialPosition);
    virtual void displaceTetrahedronPoints();
    virtual void displaceReferencePoints();
	virtual void PointsToVTKPointContainer(Points_t &points, vtkPoints* pointContainerVTK, int& numberOfPoints);
	virtual void cellsToVTKCellContainer(Ints_t &cellIds,int numberOfAlreadyExistingPoints,vtkCellArray* cellContainerVTK);
	virtual void getCellIdsForOutput(Ints_t &cellIds);

	virtual void tetrahedronToVTKContainer(int plotDisplaced, Ints_t &cellIds, int& numberOfAlreadyExistingPoints, // this is in- and output
			vtkPoints* pointContainerVTK,vtkCellArray* cellContainerVTK);
	virtual void outputFromVTKContainersToVTKFile(int plotDisplaced);

    template<class Element>
    void writeVectorInRow(std::vector<Element> & vector, std::ofstream & myfile);
    template<class Element>
    void writeMatrixInFile(std::vector<  std::vector<Element>   > & matrix, std::string filename);
    void writeShapeFunctionDerivativeInFile(ShapeFctDeriv_t & dSF, std::string filename);
    void writeMatrixInFile(StiffnessMatrix_t & matrix, int sizeOfMatrix, std::string filename);
    template<class Element>
    void writeVectorInRowOfFile(std::vector<Element> & vector, std::string filename);

	// variables
	StiffnessMatrix_t m_stiffnessMatrix;
	ConnectTriangles_t m_connectTriangle;	// for the visualization
	int m_numberOfElementNodes;
	Ints_t m_connect;

	int m_tetrahedronIsSeparated;

    Points_t m_referencePointsInTet;
    Ints_t m_referencePointIdsInTet;
    ShapeFcts_t m_referencePointShapeFunctions;
    Points_t m_displacedReferencePointsInTet;

protected:
	int m_tetrahedronId;
	int m_neighboursOfFaces[4][2];
	Ints_t m_triangleInFaceId;
	Points_t m_points;
	ShapeFcts_t m_shapeFunctions; // here we calculate the shape function of all points
	Matrix_3_3_t m_globalToLocal;

	Points_t m_GaussIntPoints;
	Doubles_t m_GaussIntWeights;
	int m_numberOfIntegrationPoints;
	int m_numberOfSegments;
	Points_t m_localIntegrationPoints;
	Points_t m_globalIntegrationPoints;
	Doubles_t m_integrationWeights;

	double m_lambda;
    double m_mu;

	Points_t m_displacement;
	Points_t m_pointsDisplaced;
};

#endif /* linearTetrahedron_H_ */
