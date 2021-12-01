// XFEM

#include <iostream>
#include <cstdio>  // for the removement of a file
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
#include <ctime>

// for tiny xml
#include <tinyxml.h>
#include <vector>
#include <ostream>

// includes Christoph PAULUS
#include "loadCutModelTopologyInVectorXFEMLinearTetrahedra.h"
#include "FEM_Elementtype.h"
#include "utils.h"

// # if (EIGEN_WORLD_VERSION < 3) || (EIGEN_MAJOR_VERSION < 2) || (EIGEN_MINOR_VERSION < 0)
// # ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
// # pragma message ("If the version of eigen is smaller 3.2.0, we define the makro EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET in order to allow for a compilation for example with the eigen version 3.0.5")
// # define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
// # endif
// # endif

// # pragma message (EIGEN_WORLD_VERSION)
// # pragma message (EIGEN_MAJOR_VERSION)
// # pragma message (EIGEN_MINOR_VERSION)

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
// #include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>



typedef Eigen::Triplet<double> T;
typedef LinearTetrahedron::Point_t Point_t; // REPLACE THE REPRESENTATION OF A POINT AS IT HAS BEEN DONE IN Cutmodel....
typedef LinearTetrahedron::Points_t Points_t;
typedef LinearTetrahedron::ShapeFct_t ShapeFct_t;

struct Plane
{
    CutModelTopology_CGAL::Point_t normal;
    CutModelTopology_CGAL::Point_t pointOfPlane;
};

std::ostream& operator<< ( std::ostream& os, const Plane & plane )
{
    os << "normal " << plane.normal << ", point " << plane.pointOfPlane;
    return os;
}

struct ROI
{
    std::vector<Plane> restrictingPlanes;
    std::vector<int> affectedPoints;
};

//std::ostream& operator<< ( std::ostream& os, const ROI & roi )
//{
//    for
//    os << "n " << plane.normal << ", p " << plane.pointOfPlane;
//    return os;
//}

struct BC
{
    std::vector<double> value;
    std::vector<Plane> restrictingPlanes;
    std::vector<int> affectedPoints;
};

std::ostream& operator<< ( std::ostream& os, const BC & bc ) 
{
    os << "BC value ";
    for (unsigned i = 0 ; i < bc.value.size() ; ++i )
        os << bc.value[i] << " ";
    os << "\nRestrictingPlanes ";
    std::vector<Plane> planes = bc.restrictingPlanes;
    //std::vector< std::vector< std::vector< double > > > planes = bc.restrictingPlanes;
    for (unsigned i = 0 ; i < planes.size() ; ++i )
    {
        /*
        for (unsigned j = 0 ; j < planes[i].size() ; ++j )
        {
            for (unsigned k = 0 ; k < (planes[i])[j].size() ; ++k )
                os << ((planes[i])[j])[k] << " ";
            if (j< (planes[i].size()-1) )
                os << ", ";
        }
        os << "\n";
        */
        os << "\n" << planes[i];
    }
    if (bc.affectedPoints.size())
        os << "AffectedPoints ";
    for (unsigned i = 0 ; i < bc.affectedPoints.size() ; ++i )
    {
        os << bc.affectedPoints[i] << " ";
    }
    return os;
}

struct Object
{
    std::string filename;
    double youngsModulus,poissonsRatio/*,mass_density*/;
    int nIntPts;
    //std::vector<double> gravity;
};

std::ostream& operator<< ( std::ostream& os, const Object & object ) 
{
    os << "Object: \nfilename " << object.filename << "\nyoungsModulus " << object.youngsModulus << ", poissonsRatio " << object.poissonsRatio /*<< ", mass_density " << object.mass_density*/;
    /*os << ", gravity : ";
    for (unsigned i = 0 ; i < object.gravity.size() ; ++i )
        os << object.gravity[i] << " ";*/
    os << ", nIntegrationPoints " << object.nIntPts;
    return os;
}

struct OutputOptions
{
    std::string outputFolder;
    int outputConsole;      // is int since it can not be read otherwise! -> search for QueryIntValue to understand that
    int outputVariablesInFile;  // is int since it can not be read otherwise! -> search for QueryIntValue to understand that
};

std::ostream& operator<< ( std::ostream& os, const OutputOptions & output ) 
{
    os << "Output: \noutputFolder " << output.outputFolder << "\noutputConsole " << output.outputConsole << ", outputVariablesInFile " << output.outputVariablesInFile ; //  << #object
}

struct CutProblem
{
    Object object;
    std::string filenameReferenceObject;
    std::string filenameCut;
    std::vector<BC> displacementBCs;
    std::vector<BC> forceBCs;
    OutputOptions outputOptions;

};

std::ostream& operator<< ( std::ostream& os, const CutProblem & cp ) 
{
    os << cp.object << "\nCut:\n " << cp.filenameCut << "\nReferenceObject:\n " << cp.filenameReferenceObject << cp.outputOptions << "\nDisplacementBCs: " << cp.displacementBCs.size() << "\n" ;
    for ( unsigned i = 0 ; i < cp.displacementBCs.size() ; ++i )
        os << cp.displacementBCs[i] << "\n";
    os << "ForceBCs: " << cp.forceBCs.size() << "\n";
    for ( unsigned j = 0 ; j < cp.forceBCs.size() ; ++j )
        os << cp.forceBCs[j]  << "\n";
    return os;
}

std::vector<double> stringToVecDouble(std::string & input)
{
    std::stringstream oss(input);
    std::vector<double> result;
    double temp = 0.0;
    while (!oss.eof())
    {
        oss >> temp;
        result.push_back(temp);
    }
    return result;
}

CutModelTopology_CGAL::Point_t stringToPoint(std::string & input)
{
    std::vector<double> const & res_vec = stringToVecDouble(input);
    return CutModelTopology_CGAL::Point_t(res_vec[0],res_vec[1],res_vec[2]);
}

Plane getPlane(TiXmlAttribute* pAttrib)
{
    Plane result;
//    DEBUG;
    while ( pAttrib )
    {
        std::string name = pAttrib->Name();
        std::string value_str = pAttrib->Value();
        //std::cout << name << " " << value_str << std::endl;
        if ( name == "normal")
            result.normal = stringToPoint(value_str);
        else if ( name == "point" )
            result.pointOfPlane = stringToPoint(value_str);
        else
            std::cerr << "Unrecognized attribute in plane: " << name << "\n";

        pAttrib = pAttrib->Next();
    }
//    DEBUG;

    if( result.normal == CutModelTopology_CGAL::Point_t(0.0,0.0,0.0) )
        std::cerr << "A plane has a normal == 0!\n";

    return result;
}

std::vector< Plane > getRestrictingPlanes(TiXmlNode* pChild)
{
    std::vector < Plane > result;
    while ( pChild )
    {
        TiXmlElement* curElement = pChild->ToElement();
        TiXmlAttribute* pAttrib = curElement->FirstAttribute();
        //std::vector< std::vector<double> > curPlane;
        Plane curPlane = getPlane(pAttrib);
        /*
        while ( pAttrib )
        {
            std::string value_str = pAttrib->Value();
            std::vector<double> const & vec = stringToVecDouble(value_str);
            curPlane.push_back(vec);
            pAttrib = pAttrib->Next();
        }
        */
//        DEBUG;
        //std::cout << curPlane << std::endl;
        result.push_back(curPlane);
        pChild = pChild->NextSibling();
    }
    return result;
}

Object getObject(TiXmlElement* pElement) 
{
    Object result;

    TiXmlAttribute* pAttrib=pElement->FirstAttribute();
    while ( pAttrib )
    {
        std::string name = pAttrib->Name();

        if ( name == "filename" )
            result.filename = pAttrib->Value(); 
        if ( name == "poissonsRatio" )
            pAttrib->QueryDoubleValue(&result.poissonsRatio);
        if ( name == "youngsModulus" )
            pAttrib->QueryDoubleValue(&result.youngsModulus);
        if ( name == "nIntPts" )
            pAttrib->QueryIntValue(&result.nIntPts);
        /*if ( name == "mass_density" )
            pAttrib->QueryDoubleValue(&result.mass_density);
        if ( name == "gravity" )
        {
            std::string gravity = pAttrib->Value();
            result.gravity = stringToVecDouble(gravity);
        }*/
        pAttrib = pAttrib->Next();
    }

    return result;
}

OutputOptions getOutputOptions(TiXmlElement* pElement) 
{
    OutputOptions result;

    TiXmlAttribute* pAttrib=pElement->FirstAttribute();
    while ( pAttrib )
    {
        std::string name = pAttrib->Name();

        if ( name == "outputFolder" )
            result.outputFolder = pAttrib->Value(); 
        if ( name == "outputConsole" )
            pAttrib->QueryIntValue(&result.outputConsole); 
        if ( name == "outputVariablesInFile" )
            pAttrib->QueryIntValue(&result.outputVariablesInFile);

        pAttrib = pAttrib->Next();
    }

    return result;
}

std::string getFirstAttribute(TiXmlElement* pElement)
{
    std::string result;
    TiXmlAttribute* pAttrib=pElement->FirstAttribute();
    result = pAttrib->Value();
    return result;
}

BC getBC(TiXmlNode* pChild)
{
    BC result;
    //std::string elementName = pChild->Value();

    TiXmlElement* curElement = pChild->ToElement();
    TiXmlAttribute* pAttrib=curElement->FirstAttribute();
    std::string a = pAttrib->Value();
    result.value = stringToVecDouble(a);
    result.restrictingPlanes = getRestrictingPlanes(pChild->FirstChild());
    return result;
}

CutProblem loadXML(const char* pFilename)
{
    CutProblem result;
    TiXmlDocument doc(pFilename);
    bool loadOkay = doc.LoadFile();
    if (loadOkay)
    {
//        printf("\n%s:\n", pFilename);
        TiXmlNode* pChild;
        TiXmlDocument* doc_ptr = &doc;
        for ( pChild = doc_ptr->FirstChild(); pChild != 0 ; pChild = pChild->NextSibling())
        {
            std::string pChildName = pChild->Value();
            TiXmlElement* curElement = pChild->ToElement();
//            SHOW(pChildName);
            //TiXmlAttribute* pAttrib=curElement->FirstAttribute();
            if ( pChildName == "Object" )
                result.object = getObject(curElement);
            else if ( pChildName == "Cut" )
                result.filenameCut = getFirstAttribute(curElement);
            else if ( pChildName == "ReferenceObject" )
                result.filenameReferenceObject = getFirstAttribute(curElement);
                //result.filenameReferenceObject = pAttrib->Value();
            else if ( pChildName == "Output" )
                result.outputOptions = getOutputOptions(curElement);                
            else if ( pChildName == "DisplacementBC" )
                result.displacementBCs.push_back( getBC ( pChild ) ); 
            else if ( pChildName == "ForceBC" )
                result.forceBCs.push_back( getBC ( pChild ) );
        }
    }
    else
    {
        printf("Failed to load file \"%s\"\n", pFilename);
    }

    return result;
}

void fillGlobalSparseStiffnessMatrix(std::vector<int> & curConnect, LinearTetrahedron::StiffnessMatrix_t & k, Eigen::SparseMatrix<double> & K)   {
    int mn, ij, MN, IJ;
    int numberOfNodes = curConnect.size();

    for(int m=0;m<numberOfNodes;m++)
        for(int n=0;n<3;n++)
            for(int i=0;i<numberOfNodes;i++)
                for(int j=0;j<3;j++)
                {
                    mn = 3*m+n;
                    ij = 3*i+j;

                    MN = 3*curConnect[m]+n;
                    IJ = 3*curConnect[i]+j;

                    K.coeffRef(IJ,MN) += k(ij,mn);
                }
}

template<class Element>
void writeVectorInRow(std::vector<Element> & vector, std::ofstream & myfile) {
    for(unsigned i=0;i<vector.size();++i)
        myfile << vector[i] << " ";
    myfile << "\n";
}

void writeSparseMatrixInFile(Eigen::SparseMatrix<double> & matrix, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::trunc);

    for(unsigned i = 0 ; i < matrix.cols(); ++i)
    {
        for(unsigned j = 0 ; j < matrix.rows(); ++j)
            myfile << setprecision(20) << matrix.coeffRef(i,j) << " ";
        myfile << "\n";
    }

    myfile.close();
}

void writeSofaMatrixInFile(LinearTetrahedron::StiffnessMatrix_t & matrix, int sizeOfMatrix, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::trunc);
	for(unsigned i=0;i<sizeOfMatrix;++i)
    {
	    for(unsigned j=0;j<sizeOfMatrix;++j)
	        myfile << setprecision(20) << matrix(i,j) << " ";
	    myfile << "\n";
	}
    myfile.close();
}

template<class Element>
void writeVectorInRowOfFile(std::vector<Element> & vector, std::string filename) {
    std::ofstream myfile;
    myfile.open (filename.c_str(),std::ios::app);
    writeVectorInRow(vector, myfile);
    myfile.close();
}

// add the following function in the three functions defined after it!!!!!!!
//bool isPointAbovePlane(Point_t & pt, Plane & plane)
//{
//    Point_t planeToPt = pt-plane.pointOfPlane;
//    double scalarProduct = planeToPt*plane.normal;
//    return scalarProduct>0;
//}

void getPointsAbovePlaneInit(CutModelTopology_CGAL::Points_t & pts, std::vector<int> & ptIds, CutModelTopology_CGAL::Point_t normalOfPlane, CutModelTopology_CGAL::Point_t ptOfPlane){
    std::vector <int> ptIdsOut;
    for(unsigned int curPtId = 0 ; curPtId < pts.size() ; ++curPtId)    // change for loop in order to obtain the possibility to
    {
        double scalProd = 0.0;
        for(unsigned j = 0 ; j < 3 ; ++j)
            scalProd += ((pts[curPtId])[j] - ptOfPlane[j]) * normalOfPlane[j];

        if(scalProd > 0.0)
            ptIdsOut.push_back(curPtId);
    }
    ptIds = ptIdsOut;
}

void getPointsAbovePlane(CutModelTopology_CGAL::Points_t & pts, std::vector<int> & ptIds, CutModelTopology_CGAL::Point_t normalOfPlane, CutModelTopology_CGAL::Point_t ptOfPlane){
    std::vector <int> ptIdsOut;
    for(unsigned int i = 0 ; i < ptIds.size() ; ++i)    // change for loop in order to obtain the possibility to
    {
        int curPtId = ptIds[i];
        double scalProd = 0.0;
        for(unsigned j = 0 ; j < 3 ; ++j)
            scalProd += ((pts[curPtId])[j] - ptOfPlane[j]) * normalOfPlane[j];

        if(scalProd > 0.0)
            ptIdsOut.push_back(curPtId);
    }
    ptIds = ptIdsOut;
}

std::vector<int> getPtIdsInRegionRestrictedByPlanes(CutModelTopology_CGAL::Points_t & pts, std::vector< Plane > & restrictingPlanes)
{
    std::vector<int> ptIdsInRegion;
    for (unsigned planeId = 0 ; planeId < restrictingPlanes.size() ; ++planeId)
    {
        Plane curPlane = restrictingPlanes[planeId];
        if ( planeId == 0 )
            getPointsAbovePlaneInit(pts,ptIdsInRegion,curPlane.normal,curPlane.pointOfPlane);
        else
            getPointsAbovePlane(pts,ptIdsInRegion,curPlane.normal,curPlane.pointOfPlane);
    }
    return ptIdsInRegion;
}

void changeVectorEntries(std::vector<int> & ptIds, std::vector<double> & newEntries, bool add, Eigen::VectorXd & f) { 
	// MAYBE REPLACE newEntries by a Eigen::VectorXd as well
	for(unsigned it = 0 ; it < ptIds.size() ; ++it)
	{
		int curPtId = ptIds[it];

        for(unsigned k = 0 ; k < newEntries.size() ; ++k)
			f(curPtId*newEntries.size() + k) = f(curPtId*newEntries.size() + k)*add + newEntries[k];
	}
}

void changeSparseMatrixEntries(std::vector<int> & ptIds, Eigen::SparseMatrix<double> & A)	{

    for(unsigned it = 0 ; it < ptIds.size() ; ++it)
    {
        unsigned curPtId = ptIds[it];
        for(unsigned i = 0 ; i < A.rows() ; ++i)
            for(unsigned j = 0 ; j < 3 ; ++j)
            {
                A.coeffRef(curPtId*3+j,i) = 0.0;
                A.coeffRef(i,curPtId*3+j) = 0.0;
            }
    }
    for(unsigned it = 0 ; it < ptIds.size() ; ++it)
    {
        unsigned curPtId = ptIds[it];
        for(unsigned j = 0 ; j < 3 ; ++j)
            A.coeffRef(curPtId*3+j,curPtId*3+j) = 1.0;
    }
}

void assignDisplBC(BC & displBC,std::vector<int> & enforcedBC)
{

}

vtkUnstructuredGrid* getUnstructuredGrid(std::string & filename)
{
    // define the reader, the grid, the points and the array, to receive information read from the vtk file
    vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
    reader->SetFileName(filename.c_str());
    reader->Update();
    vtkUnstructuredGrid* myGrid = reader->GetOutput();

    return myGrid;
}

std::vector<unsigned> calculateShapeFunctionsOfReferenceMesh(vtkUnstructuredGrid* referenceMesh, std::vector<XFemLinearTetrahedron> & simulatedObject)
{
    vtkPoints* referencePoints = referenceMesh->GetPoints();
    int numberOfPoints = int(referencePoints->GetNumberOfPoints());
    int numberOfCells = int(referenceMesh->GetNumberOfCells());

    std::vector<unsigned> result;
    result.resize(numberOfPoints);
    std::set<unsigned> assignedPtIds;
    double ptGlobalDouble[3];
    Point_t ptGlobal;
    Point_t ptLocal;
    ROI smallRegionAroundCut;
    double epsilon = 0.0000001d;
    Plane belowCut;
    belowCut.normal = CutModelTopology_CGAL::Point_t(0,0,-1);
    belowCut.pointOfPlane = CutModelTopology_CGAL::Point_t(0,0,0.1+epsilon);
    smallRegionAroundCut.restrictingPlanes.push_back(belowCut);
    Plane aboveCut;
    aboveCut.normal = CutModelTopology_CGAL::Point_t(0,0,1);
    aboveCut.pointOfPlane = CutModelTopology_CGAL::Point_t(0,0,0.1-epsilon);
    smallRegionAroundCut.restrictingPlanes.push_back(aboveCut);
    Plane cutFront;
    cutFront.normal = CutModelTopology_CGAL::Point_t(-1,0,0);
    cutFront.pointOfPlane = CutModelTopology_CGAL::Point_t(0.015-epsilon,0,0.1);
    smallRegionAroundCut.restrictingPlanes.push_back(cutFront);
    CutModelTopology_CGAL::Points_t pts;
    pts.resize(numberOfPoints);
    for (unsigned ptId = 0 ; ptId < numberOfPoints ; ++ptId )
    {
        referencePoints->GetPoint(ptId, ptGlobalDouble);
        pts[ptId] = CutModelTopology_CGAL::Point_t(ptGlobalDouble[0],ptGlobalDouble[1],ptGlobalDouble[2]);
    }
    smallRegionAroundCut.affectedPoints = getPtIdsInRegionRestrictedByPlanes(pts,smallRegionAroundCut.restrictingPlanes);
    vtkCellArray* referenceCellArray = referenceMesh->GetCells();

    vtkIdType numberOfCellNodes = 4;
    vtkIdType* currentElement = new vtkIdType[4];
    for (unsigned i = 0 ; i < smallRegionAroundCut.affectedPoints.size() ; ++i)
    {
        int curPtId = smallRegionAroundCut.affectedPoints[i];
        CutModelTopology_CGAL::Point_t curPt = pts[curPtId];
        bool foundElementThatContainsPtId = false;
        for(    vtkIdType elementId=0;
                elementId<numberOfCells*(numberOfCellNodes+1) && !foundElementThatContainsPtId;
                elementId+=(numberOfCellNodes+1))
        {
            referenceCellArray->GetCell(elementId,numberOfCellNodes,currentElement);
            for(int i=0;i<numberOfCellNodes;i++)
            {
                int currentNodeId = (int)currentElement[i];
                if (currentNodeId == curPtId)
                {
                    foundElementThatContainsPtId = true;
                    int j = (i+1)%numberOfCellNodes;
                    double zEntries = 0.0;
                    while ( j!= i )
                    {
                        int otherNodeIdInCell = (int)currentElement[j];
                        zEntries += (pts[otherNodeIdInCell])[2];
                        j = (j+1)%numberOfCellNodes;
                    }
                    zEntries /= 3.0;
                    bool isCellAboveCut = zEntries > curPt[2] ;
                    CutModelTopology_CGAL::Point_t displaceCurPt(0,0,2*(isCellAboveCut-0.5)*epsilon);
                    CutModelTopology_CGAL::Point_t newPosOfCurPt(curPt[0]+displaceCurPt[0],curPt[1]+displaceCurPt[1],curPt[2]+displaceCurPt[2]);
                    referencePoints->SetPoint(curPtId,newPosOfCurPt[0],newPosOfCurPt[1],newPosOfCurPt[2]);
                }
            }
        }
    }
    for (unsigned elId = 0 ; elId<simulatedObject.size() ; ++elId )
    {
        XFemLinearTetrahedron & currentElement = simulatedObject[elId];
        for (unsigned ptId = 0 ; ptId < numberOfPoints ; ++ptId )
        {
            if(     (assignedPtIds.size() == 0)
                ||  (assignedPtIds.find(ptId) == assignedPtIds.end())) // check whether we already found the element for this point id
            {
                referencePoints->GetPoint(ptId, ptGlobalDouble);
                ptGlobal = Point_t(ptGlobalDouble[0],ptGlobalDouble[1],ptGlobalDouble[2]);
                currentElement.transformGlobalToLocalCoordinatesByUsingLinearShapeFunctions(ptGlobal,ptLocal);
                bool ptIsInTetrahedron = (std::abs(ptLocal[0])+std::abs(ptLocal[1])+std::abs(ptLocal[2])) <= 1+epsilon;
                for (unsigned i = 0 ; i < DIM && ptIsInTetrahedron; ++i)
                    if ( ptLocal[i] < -epsilon)
                        ptIsInTetrahedron = false;
                if (ptIsInTetrahedron)
                {
                    assignedPtIds.insert(ptId);
                    currentElement.m_referencePointIdsInTet.push_back(ptId);
                    CutModelTopology_CGAL::Point_t currentPt(ptGlobalDouble[0],ptGlobalDouble[1],ptGlobalDouble[2]);
                    int pointIsAbove = currentElement.m_cutModel->pointAbove(currentPt);
                    double r,phi;
                    if(currentElement.m_branchEnrichedVertexIds.size())
                        currentElement.m_cutModel->getRAndPhiOfPointDependentOnCutFront(currentPt,r,phi);
                    ShapeFct_t currentShapeFunction;
                    currentElement.getShapeFunction(ptGlobal,pointIsAbove,r,phi,currentShapeFunction);
                    currentElement.m_referencePointsInTet.push_back(ptGlobal);
                    currentElement.m_referencePointShapeFunctions.push_back(currentShapeFunction);
                    result[ptId] = elId;
                }
            }
        }
    }

    return result;
}

void deformReferenceMeshAndWriteItToFile(vtkUnstructuredGrid* referenceMesh, std::vector<unsigned> & refPtInObjElId, std::vector<XFemLinearTetrahedron> & simulatedObject, std::string filename)
{
    vtkPoints* refPoints = referenceMesh->GetPoints();
    int numberOfPoints = int(refPoints->GetNumberOfPoints());
    double currentPointDouble[3];
    for (unsigned ptId = 0 ; ptId < numberOfPoints ; ++ptId )
    {
        int elId = refPtInObjElId[ptId];
        XFemLinearTetrahedron curEl = simulatedObject[ elId ];
        unsigned id = 0;
        for ( ; id < curEl.m_referencePointIdsInTet.size() && curEl.m_referencePointIdsInTet[id] != ptId ; ++id) {}
        Point_t displacedPoint = curEl.m_displacedReferencePointsInTet[id];
        for(int i=0;i<3;i++)
            currentPointDouble[i] = displacedPoint[i];
        refPoints->SetPoint(ptId,displacedPoint[0],displacedPoint[1],displacedPoint[2]);
    }

    vtkUnstructuredGridWriter* writerTopoDebug = vtkUnstructuredGridWriter::New();
    writerTopoDebug->SetInputData(referenceMesh);
    writerTopoDebug->SetFileName(filename.c_str());
    writerTopoDebug->Write();
}

void lanceScene(CutProblem& curCutProblem) 
{
    std::vector<double> times;
    std::vector<unsigned> timesBetweenLines;
    boost::timer timer;
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    std::string m_ModelFilename = curCutProblem.object.filename;
    std::string m_CutFilename = curCutProblem.filenameCut;
    std::string referenceSolutionFilename = curCutProblem.filenameReferenceObject;

    vtkUnstructuredGrid* referenceMesh;
    if (referenceSolutionFilename != std::string() )
        referenceMesh = getUnstructuredGrid(referenceSolutionFilename);

    int numberOfIntegrationPointsXFEMTetrahedron = curCutProblem.object.nIntPts;
    double E = curCutProblem.object.youngsModulus;
    double nu = curCutProblem.object.poissonsRatio;
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2*(1+nu));

    loadCutModelTopologyInVectorXFEMLinearTetrahedra m_loadCutModelInVectorXFEMLinearTetrahedra;
    std::vector<XFemLinearTetrahedron> m_VectorXFemLinearTetrahedra;
    std::vector<int> m_pointsReferenceSolutionInObjectTetrahedronId;
    std::vector<int> m_pointsReferenceSolutionAbove;
    std::vector<double> m_pointsReferenceSolution_r;
    std::vector<double> m_pointsReferenceSolution_phi;
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    m_loadCutModelInVectorXFEMLinearTetrahedra.init(m_ModelFilename,m_CutFilename,referenceSolutionFilename,
                lambda,mu,numberOfIntegrationPointsXFEMTetrahedron,
                m_pointsReferenceSolutionInObjectTetrahedronId,m_pointsReferenceSolutionAbove,
                m_pointsReferenceSolution_r, m_pointsReferenceSolution_phi,
                m_VectorXFemLinearTetrahedra);
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);


    CutModelTopology_CGAL::Points_t m_objectPoints = m_loadCutModelInVectorXFEMLinearTetrahedra.m_cutModel.m_objectPoints;
    std::vector<unsigned> m_refPtInObjElId;
    if (referenceSolutionFilename != std::string() )
        m_refPtInObjElId = calculateShapeFunctionsOfReferenceMesh(referenceMesh,m_VectorXFemLinearTetrahedra);

    // std::cout << m_pointsReferenceSolutionInObjectTetrahedronId.size() << std::endl;
    // std::cout << m_VectorXFemLinearTetrahedra.size() << std::endl;
    int numberOfOriginalDOFs = 3*m_loadCutModelInVectorXFEMLinearTetrahedra.getNumberOfPoints();
    int numberOfSignDOFs = 3*m_loadCutModelInVectorXFEMLinearTetrahedra.getNumberOfSignEnrichedPoints();
    int numberOfBranchDOFs = 3*m_loadCutModelInVectorXFEMLinearTetrahedra.getNumberOfBranchEnrichedPoints();
    int numberOfDOFs = numberOfOriginalDOFs + numberOfSignDOFs + numberOfBranchDOFs ;

    Eigen::SparseMatrix<double> K(numberOfOriginalDOFs,numberOfOriginalDOFs);

    if (curCutProblem.outputOptions.outputConsole) // HAS TO STAY AS OUTPUT!
    {
        SHOW(numberOfOriginalDOFs);
        SHOW(numberOfSignDOFs);
        SHOW(numberOfBranchDOFs);
        SHOW(numberOfDOFs);
    }

    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

//    if (curCutProblem.outputOptions.outputConsole)
//        std::cout << "Compute stiffness matrix " << std::endl;
    std::string filenameConnect = curCutProblem.outputOptions.outputFolder + "connect";
    if (curCutProblem.outputOptions.outputVariablesInFile)
        SHOW( std::remove (filenameConnect.c_str() ) == 0 );
    for(unsigned int i = 0 ; i < m_VectorXFemLinearTetrahedra.size() ; ++i)
    {
        XFemLinearTetrahedron & curTet = m_VectorXFemLinearTetrahedra[i];
        int sizeStiffnessmatrix = curTet.m_connect.size()*3;
        
        if (curCutProblem.outputOptions.outputVariablesInFile)
            writeVectorInRowOfFile(curTet.m_connect,filenameConnect);

        curTet.stiffnessBlock2();
//        if (curCutProblem.outputOptions.outputConsole)
//        {
//            std::cout << i << " : " << curTet.m_signEnrichedVertexIds.size() << " " << 4*curTet.m_branchEnrichedVertexIds.size() << " c "; // comment protected in XFEMLinearT...
//            for(unsigned j=0; j<curTet.m_connect.size() ; ++j)
//                std::cout << " " << curTet.m_connect[j];
//            std::cout << std::endl;
//        }

        std::string filenameLmnStiffness = curCutProblem.outputOptions.outputFolder + "stiffness/" + boost::lexical_cast<std::string>(i);
        if (curCutProblem.outputOptions.outputVariablesInFile)
            writeSofaMatrixInFile(curTet.m_stiffnessMatrix,sizeStiffnessmatrix,filenameLmnStiffness);

        fillGlobalSparseStiffnessMatrix(curTet.m_connect, curTet.m_stiffnessMatrix, K);
    }
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

    std::string fNGM2 = curCutProblem.outputOptions.outputFolder + "stiffness/global";
    if (curCutProblem.outputOptions.outputVariablesInFile)
        writeSparseMatrixInFile(K,fNGM2);

    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

    Eigen::VectorXd u;
    Eigen::VectorXd u_f;
    Eigen::VectorXd f;
    u.resize(numberOfOriginalDOFs);
    u_f.resize(numberOfOriginalDOFs);
    f.resize(numberOfOriginalDOFs);
    for(unsigned i = 0 ; i < f.size() ; ++i)
    {
        f(i) = 0.0;
        u(i) = 0.0;
        u_f(i) = 0.0;
    }

    // calculate the points being affected by the boundary conditions
    //std::set<int> ptsWithDisplBCs;
    std::vector<int> ptsWithBCs;

    Eigen::VectorXd u_f_new(numberOfOriginalDOFs);
    Eigen::VectorXd f_new(numberOfOriginalDOFs);
    for(unsigned i = 0 ; i < f.size() ; ++i)
    {
        f_new(i) = 0.0;
        u_f_new(i) = 0.0;
    }
    // the following two for loops can be put into the same function!!!!! just with another input: u_f_new
    for ( unsigned displBCId = 0 ; displBCId<curCutProblem.displacementBCs.size() ; ++displBCId)
    {
        BC & curDisplacementBC = curCutProblem.displacementBCs[displBCId];
        std::vector<int> curAffectedPoints = getPtIdsInRegionRestrictedByPlanes(m_objectPoints,curDisplacementBC.restrictingPlanes);
        curDisplacementBC.affectedPoints = curAffectedPoints;
        changeVectorEntries(curAffectedPoints, curDisplacementBC.value, true, u_f_new); // maybe put false!
        ptsWithBCs.insert(ptsWithBCs.end(),curAffectedPoints.begin(),curAffectedPoints.end());
    }
    f_new = -K*u_f_new;
    std::vector < double > newEntries;
    newEntries.push_back(0.0);
    newEntries.push_back(0.0);
    newEntries.push_back(0.0);
    changeVectorEntries(ptsWithBCs, newEntries, false, f_new);
    for ( unsigned forceBCId = 0 ; forceBCId<curCutProblem.forceBCs.size() ; ++forceBCId )
    {
        BC & curForceBC = curCutProblem.forceBCs[forceBCId];
        std::vector<int> curAffectedPoints = getPtIdsInRegionRestrictedByPlanes(m_objectPoints,curForceBC.restrictingPlanes);
        curForceBC.affectedPoints = curAffectedPoints;

        changeVectorEntries(curAffectedPoints, curForceBC.value, true, f_new);
        ptsWithBCs.insert(ptsWithBCs.end(),curAffectedPoints.begin(),curAffectedPoints.end());
    }

    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    if (curCutProblem.outputOptions.outputConsole)
    {
        std::cout << "Point ids with boundary conditions: ";
        for (unsigned i = 0 ; i < ptsWithBCs.size() ; ++i)
            std::cout << ptsWithBCs[i] << " ";
        std::cout << std::endl;
    }
    // TD: cause error if there are two identical entries in a vector
    changeSparseMatrixEntries(ptsWithBCs,K);
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    if (curCutProblem.outputOptions.outputVariablesInFile)
        writeSparseMatrixInFile(K,curCutProblem.outputOptions.outputFolder + "stiffness/global_BC");

    // solve the linear system using the conjugate gradient
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    Eigen::MatrixXd K2 = K;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(K2);
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    u = solver.solve(f_new+u_f_new);

    std::vector<double> f_vec;
    std::vector<double> BC_vec;
    for(unsigned i = 0 ; i < f_new.size() ; ++i)
    {
        f_vec.push_back(f_new[i]);
        BC_vec.push_back(f_new[i]+u_f_new[i]);
    }
    if (curCutProblem.outputOptions.outputVariablesInFile)
    {
        std::string filenameRHS = curCutProblem.outputOptions.outputFolder + "RHS";
        SHOW( std::remove (filenameRHS.c_str() ) == 0 );
        writeVectorInRowOfFile(f_vec, filenameRHS);
        std::string filenameBC = curCutProblem.outputOptions.outputFolder + "BC";
        SHOW( std::remove (filenameBC.c_str() ) == 0 );
        writeVectorInRowOfFile(BC_vec, filenameBC);
    }

//    for (unsigned i = 0 ; i < u.size() ; ++i )
//        std::cout << u[i] << " " << u_f_new[i] << " " << f_new[i] << std::endl;
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);
    for ( unsigned i = 0 ; i < u_f.size() && false; ++i )
        if (u_f[i])
            u[i] = u_f[i];
    std::vector<double> u_vec;
    for(unsigned i = 0 ; i < u.size() ; ++i)
        u_vec.push_back(u[i]);

    if (curCutProblem.outputOptions.outputVariablesInFile)
    {
        std::string filenameSolution = curCutProblem.outputOptions.outputFolder + "solution";
        SHOW( std::remove (filenameSolution.c_str() ) == 0 );
        writeVectorInRowOfFile(u_vec, filenameSolution  );
    }

    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

    int  nElements = m_VectorXFemLinearTetrahedra.size();
    //DISPLACE THE ELEMENTS
    for(unsigned int iterElement = 0 ; iterElement < m_VectorXFemLinearTetrahedra.size() ; ++iterElement)
    {
        XFemLinearTetrahedron & curTet = m_VectorXFemLinearTetrahedra[iterElement];
        int numberOfPoints = curTet.m_numberOfElementNodes;
        LinearTetrahedron::Points_t currentDisplacements;

        LinearTetrahedron::Ints_t curConnect = curTet.m_connect;

        LinearTetrahedron::Point_t currentDisp;
        for(int i=0; i<numberOfPoints; i++)
        {
            currentDisp[0] = u[3*curConnect[i] + 0];
            currentDisp[1] = u[3*curConnect[i] + 1];
            currentDisp[2] = u[3*curConnect[i] + 2];
            currentDisplacements.push_back(currentDisp);
        }

        curTet.writeDisplacementInTetrahedronData(currentDisplacements);
        curTet.displaceTetrahedronPoints();
        curTet.displaceReferencePoints();
    }

    std::string referenceMeshFilename = curCutProblem.outputOptions.outputFolder + "displacedObject/HighRes.vtk";

    if (referenceSolutionFilename != std::string() )
        deformReferenceMeshAndWriteItToFile(referenceMesh, m_refPtInObjElId, m_VectorXFemLinearTetrahedra, referenceMeshFilename);
    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

    // WRITE THE ELEMENTS INTO FILES
    //write output for displaced object
    vtkPoints* pointContainerVTK = vtkPoints::New();
    vtkCellArray* triangleContainerVTK = vtkCellArray::New();
    vtkPoints* pointContainerVTK_uncut = vtkPoints::New();
    vtkCellArray* triangleContainerVTK_uncut = vtkCellArray::New();
    vtkPoints* pointContainerVTK_completeCutBelow = vtkPoints::New();
    vtkCellArray* triangleContainerVTK_completeCutBelow = vtkCellArray::New();
    vtkPoints* pointContainerVTK_completeCutAbove = vtkPoints::New();
    vtkCellArray* triangleContainerVTK_completeCutAbove = vtkCellArray::New();
    vtkPoints* pointContainerVTK_partialCutBelow = vtkPoints::New();
    vtkCellArray* triangleContainerVTK_partialCutBelow = vtkCellArray::New();
    vtkPoints* pointContainerVTK_partialCutAbove = vtkPoints::New();
    vtkCellArray* triangleContainerVTK_partialCutAbove = vtkCellArray::New();


    int numberOfAlreadyExistingPoints = 0;
    int numberOfAlreadyExistingPoints_uncut = 0;
    int numberOfAlreadyExistingPoints_cc_b = 0;
    int numberOfAlreadyExistingPoints_cc_a = 0;
    int numberOfAlreadyExistingPoints_pc_b = 0;
    int numberOfAlreadyExistingPoints_pc_a = 0;
    int plotdisplaced = 1;
    for(unsigned int iterElement=0;iterElement<nElements;iterElement++)
    {
        if( m_VectorXFemLinearTetrahedra[iterElement].m_tetrahedronIsSeparated==0 &&
            m_VectorXFemLinearTetrahedra[iterElement].m_tetrahedronIsPartiallySeparated==0)
        {
            LinearTetrahedron::Ints_t cellIds;
            m_VectorXFemLinearTetrahedra[iterElement].getCellIdsForOutput(-1,cellIds);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds,numberOfAlreadyExistingPoints_uncut,
                                pointContainerVTK_uncut,triangleContainerVTK_uncut);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds,numberOfAlreadyExistingPoints,pointContainerVTK,triangleContainerVTK);
        }
        else if (m_VectorXFemLinearTetrahedra[iterElement].m_tetrahedronIsSeparated==1)
        {
            LinearTetrahedron::Ints_t cellIds_b;
            m_VectorXFemLinearTetrahedra[iterElement].getCellIdsForOutput(0,cellIds_b);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_b,numberOfAlreadyExistingPoints_cc_b,
                    pointContainerVTK_completeCutBelow,triangleContainerVTK_completeCutBelow);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_b,numberOfAlreadyExistingPoints,pointContainerVTK,triangleContainerVTK);
            LinearTetrahedron::Ints_t cellIds_a;
            m_VectorXFemLinearTetrahedra[iterElement].getCellIdsForOutput(1,cellIds_a);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_a,numberOfAlreadyExistingPoints_cc_a,
                    pointContainerVTK_completeCutAbove,triangleContainerVTK_completeCutAbove);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_a,numberOfAlreadyExistingPoints,pointContainerVTK,triangleContainerVTK);
        }
        else if (m_VectorXFemLinearTetrahedra[iterElement].m_tetrahedronIsPartiallySeparated==1)
        {
            LinearTetrahedron::Ints_t cellIds_b;
            m_VectorXFemLinearTetrahedra[iterElement].getCellIdsForOutput(0,cellIds_b);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_b,numberOfAlreadyExistingPoints_pc_b,
                    pointContainerVTK_partialCutBelow,triangleContainerVTK_partialCutBelow);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_b,numberOfAlreadyExistingPoints,pointContainerVTK,triangleContainerVTK);
            LinearTetrahedron::Ints_t cellIds_a;
            m_VectorXFemLinearTetrahedra[iterElement].getCellIdsForOutput(1,cellIds_a);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_a,numberOfAlreadyExistingPoints_pc_a,
                    pointContainerVTK_partialCutAbove,triangleContainerVTK_partialCutAbove);
            m_VectorXFemLinearTetrahedra[iterElement].tetrahedronToVTKContainer(plotdisplaced /*plotdisplaced */,cellIds_a,numberOfAlreadyExistingPoints,pointContainerVTK,triangleContainerVTK);
        }
    }

    int counter = 0;
    // define output filename
    std::string baseFilename = curCutProblem.outputOptions.outputFolder + "displacedObject/";
    vtkUnstructuredGrid* myGridTopoDebug = vtkUnstructuredGrid::New();
    vtkUnstructuredGridWriter* writerTopoDebug = vtkUnstructuredGridWriter::New();

    // all elements
    std::string completeFilename(baseFilename);
    completeFilename += "XFEM_Visualization";
//    completeFilename += boost::lexical_cast<std::string>(counter);
    completeFilename += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename.c_str());
    writerTopoDebug->Write();
    // uncut elements
    std::string completeFilename_uncut(baseFilename);
    completeFilename_uncut += "XFEM_Visualization_uncutElements";//_";
//    completeFilename += boost::lexical_cast<std::string>(counter);
    completeFilename_uncut += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK_uncut);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK_uncut);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename_uncut.c_str());
    writerTopoDebug->Write();
    // completely cut elements below the cut
    std::string completeFilename_cc_b(baseFilename);
    completeFilename_cc_b+= "XFEM_Visualization_completelyCut_below";//_";
//    completeFilename_cc_b += boost::lexical_cast<std::string>(counter);
    completeFilename_cc_b += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK_completeCutBelow);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK_completeCutBelow);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename_cc_b.c_str());
    writerTopoDebug->Write();
    // completely cut elements above the cut
    std::string completeFilename_cc_a(baseFilename);
    completeFilename_cc_a += "XFEM_Visualization_completelyCut_above";//_";
//    completeFilename_cc_a += boost::lexical_cast<std::string>(counter);
    completeFilename_cc_a += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK_completeCutAbove);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK_completeCutAbove);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename_cc_a.c_str());
    writerTopoDebug->Write();
    // partially cut elements below the cut
    std::string completeFilename_pc_b(baseFilename);
    completeFilename_pc_b+= "XFEM_Visualization_partiallyCut_below";//_";
//    completeFilename_pc_b += boost::lexical_cast<std::string>(counter);
    completeFilename_pc_b += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK_partialCutBelow);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK_partialCutBelow);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename_pc_b.c_str());
    writerTopoDebug->Write();
    // partially cut elements above the cut
    std::string completeFilename_pc_a(baseFilename);
    completeFilename_pc_a+= "XFEM_Visualization_partiallyCut_above";//_";
//    completeFilename_pc_a += boost::lexical_cast<std::string>(counter);
    completeFilename_pc_a += ".vtk";
    myGridTopoDebug->SetPoints(pointContainerVTK_partialCutAbove);
    myGridTopoDebug->SetCells(VTK_TRIANGLE, triangleContainerVTK_partialCutAbove);
    writerTopoDebug->SetInputData(myGridTopoDebug);
    writerTopoDebug->SetFileName(completeFilename_pc_a.c_str());
    writerTopoDebug->Write();

    pointContainerVTK_uncut->Delete();
    triangleContainerVTK_uncut->Delete();
    myGridTopoDebug->Delete();
    writerTopoDebug->Delete();

    times.push_back(timer.elapsed());
    timesBetweenLines.push_back(__LINE__);

    if (curCutProblem.outputOptions.outputConsole)
    {
        std::cout << "computation times between lines :\n";
        for(unsigned i = 1 ; i < times.size() ; ++i)
            std::cout << timesBetweenLines[i-1] << " - " << timesBetweenLines[i] << " : " << times[i] - times[i-1] << std::endl;
        std::cout << "overall computation time : " << times[times.size()-1] - times[0] << std::endl;
    }

    //std::cout << std::endl << "Salut XFEM!" << std::endl;
}

int main (int argc, char *argv[])
{
    for (int i=1; i<argc; i++)
    {
        std::cout << argv[i] << std::endl;
        CutProblem curCutProblem = loadXML(argv[i]);
        if (curCutProblem.outputOptions.outputConsole)
        {
            std::cout << "_______________________________________________________________" << std::endl;
            std::cout << curCutProblem << std::endl;
        }
        lanceScene(curCutProblem);
        if (curCutProblem.outputOptions.outputConsole)
            std::cout << "_______________________________________________________________" << std::endl;
    }
	return 0;
}
