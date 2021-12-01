/*
 * XFemLinearTetrahedron.h
 *
 *  Created on: July 10, 2013
 *      Author: Christoph Paulus
 */

//defines
#ifndef XFemLinearTetrahedron_H_
#define XFemLinearTetrahedron_H_

// includes

#include "CutModelTopology_CGAL.h"
#include "LinearTetrahedron.h"

#include <Eigen/Core>
#include <Eigen/LU>

#include <math.h>	// for floor function

/*
class XFemLinearTetrahedron : provides all functions and variables that are necessary, on top of the functions in LinearTetrahedron, to handle enriched elements 
- shape functions of the enriched degrees of freedom
- derivatives of the shape functions of the enriched degrees of freedom
*/

class XFemLinearTetrahedron : public LinearTetrahedron
{
public:
	//typedefs - see LinearTetrahedron.h

	// functions
	XFemLinearTetrahedron(void);
	virtual ~XFemLinearTetrahedron();
	virtual void init(int tetrahedronId,Ints_t connect,int neighboursOfFaces[][2],Points_t &points, Ints_t pointIsAbove, Ints_t &signEnrichedVertexIds, Ints_t &branchEnrichedVertexIds,
			CutModelTopology_CGAL *cutModel,ConnectTriangles_t &connectTriangle, Ints_t &triangleInFaceId, Ints_t &triangleAbove,
			double lambda,double mu,int numberOfIntegrationPoints);
	virtual void getInformationAboutPartialAndCompleteCut();
	virtual void getShapeFunctionsOf_m_points();
	virtual bool isEnriched();
	virtual void getShapeFunctionSignEnrichedNode(Point_t x, ShapeFct_t &ShapeFunction);
	virtual void getShapeFunctionSignEnrichedNode2(int pointIsAbove, ShapeFct_t ShapeFunctionIn, ShapeFct_t &ShapeFunctionOut);
	virtual void getShapeFunctionBranchEnrichedNode(double r, double phi, ShapeFct_t ShapeFunctionIn, ShapeFct_t &ShapeFunctionOut);
	virtual void getAsymptoticCrackTipFunctions(double r, double phi, Doubles_t &ACTFs);
	virtual void getShapeFunction(Point_t x, int pointIsAbove, double r, double phi, ShapeFct_t &ShapeFunction);
	virtual void derivativeACTFsWithRespectToRAndPhi(double r,double phi,std::vector<Doubles_t> &dACTFs_drdphi);
	virtual void getDerivativeACTFs(int integrationPointNumber, ShapeFctDeriv_t &dACTFs);
	virtual void getShapeFunctionDerivativesDependentOnIntegrationPointNumber(int integrationPointNumber, ShapeFctDeriv_t &ShapeFunctionDerivative);
	virtual void getCellIdsForOutput(int above, Ints_t &cellIds);

	int m_tetrahedronIsSeparated;
	int m_tetrahedronIsPartiallySeparated;

//protected:
	Ints_t m_pointIsAbove;
	Doubles_t m_points_r;
	Doubles_t m_points_phi;
	Ints_t m_signEnrichedVertexIds;
	Ints_t m_branchEnrichedVertexIds;

	int m_edgeIsSplit[6];
	int m_faceIsSeparated[4];
	int m_faceIsPartiallySeparated[4];

	Ints_t m_triangleAbove;

    CutModelTopology_CGAL *m_cutModel;  // GET RID OF THIS!!!! DESIGN A STRUCT WITH A VECTOR OF XFEMLinear... and the cut model!

	Ints_t m_integrationPointIsAbove;
	Doubles_t m_integrationPoint_r;
	Doubles_t m_integrationPoint_phi;
	Doubles_t m_integrationPoint_x_dependentOnCutFront;
	Doubles_t m_integrationPoint_y_dependentOnCutFront;
	Ints_t m_integrationPointClosestToCutBoundaryElementId;
};

#endif /* XFemLinearTetrahedron_H_ */
