#ifndef MESH_HPP
#define MESH_HPP

#include <cstdio>
#include <vector>
#include <string>

// The element class contains all information relative to an element type in the mesh. Its constructor allows to
// create everything required for the DG solver to work. Namely, it builds the mass and stiffness matrices,
// while also computing the neighbours of the elements that belong to a specific element type. Two parameters
// only are required for the initialization to occur: the gauss type of the integration performed on the element
// , and the element type.

class Frontier;

class Element{

    protected:

    std::string m_name; // Name of the element
    int m_type; // Type of the element
    int m_dim; // dimension of the element
    std::vector<std::size_t> m_tag; // Tag of the elements of the given type.
    std::vector<std::size_t> m_nodeTags; // Node tags of the element types.
    int m_numNodes; // Number of nodes on an element of the given type.
    int m_order; // element order
    std::vector<double> m_massMatrixInverse; // Mass matrix of the element
    std::vector<double> m_stiffnessMatrixX; // stiffness matrix of the element along X.
    std::vector<double> m_stiffnessMatrixY; // stiffness matrix of the element along Y.
    std::vector<double> m_stiffnessMatrixZ; // stiffness matrix of the element along Z.
    std::vector<double> m_gaussPoints; // Gauss points of the element.
    std::vector<double> m_gaussWeight; // Gauss points of the element.
    int m_gaussPointsNumber; // NUmber of Gauss points on the reference element.
    int m_totalGaussPointsNumber; // Total number of gauss points on the mesh.
    std::vector<double> m_jacobian; // Jacobian of the element at each gauss points.
    std::vector<double> m_jacobianInverse; // Jacobian inverse of the element at each gauss points.
    std::vector<std::pair<int,int>> m_neighbours; // Dim and tag of the element neighbour.
    std::vector<std::size_t> m_frontierNodes; // Nodes at the frontier of the element.
    int m_numFrontier; // number of frontiers on the element.
    int m_numNodesPerFrontier; // Number of nodes per frontier of an element.
    int m_numFrontierNodes; // Total number of nodes on all frontier of an element.
    std::vector<double> m_basisFunctions; // Basis functions of the element.
    int m_basisFunctionsNumber; // Number of basis functions.
    int m_basisFunctionsCompo; // Number of components of basis functions.
    std::vector<double> m_gradBasisFunctions; // Gradient of basis functions of the element.
    std::vector<double> m_gradRealBasisFunctionsX; // Gradient of basis functions of the element in the real coordinate in X.
    std::vector<double> m_gradRealBasisFunctionsY; // Gradient of basis functions of the element in the real coordinate in Y.
    std::vector<double> m_gradRealBasisFunctionsZ; // Gradient of basis functions of the element in the real coordinate in Z.
    int m_gradBasisFunctionsCompo; // Number of components of gradient of basis functions.
    std::vector<double> m_jacobianMatrix; // Jacobian matrix of an element.
    Frontier * m_frontierElements = NULL; // Contains the frontier elements created.
    std::vector<std::pair<int,int>> m_nodeCorrespondance; // Retains the indices of the nodes (not the tags) of the main elements that correspond to specific nodes on the frontier.
    
    public :

    Element(const std::string gaussType, int type, int entityTag);
    ~Element();
    

    std::vector<double> getBasisFunctions()
    {
        return m_basisFunctions;
    }

    int getType()
    {
        return m_type;
    }

    std::vector<double> getGaussWeight()
    {
        return m_gaussWeight;
    }

    std::vector<std::pair<int,int>> getCorrespondance()
    {
        return m_nodeCorrespondance;
    }

    std::vector<double> getMassMatrixInverse()
    {
        return m_massMatrixInverse;
    }

    double getStiffnessMatrixX(int i)
    {
        return m_stiffnessMatrixX[i];
    }

    double getStiffnessMatrixY(int i)
    {
        return m_stiffnessMatrixY[i];
    }

    double getStiffnessMatrixZ(int i)
    {
        return m_stiffnessMatrixZ[i];
    }

    std::size_t getGaussPointNumber()
    {
        return m_gaussPointsNumber;
    }

    std::size_t getNumNodes()
    {
        return m_numNodes;
    }

    std::size_t getTotalGaussPointNumber()
    {
        return m_totalGaussPointsNumber;
    }

    std::size_t getTotalNumNodes()
    {
        return m_nodeTags.size();
    }

    Frontier* getFrontierElements()
    {
        return m_frontierElements;
    }

    std::vector<std::pair<int,int>> getNeighbours()
    {
        return m_neighbours;
    }

    std::vector<std::size_t> getNodeTags()
    {
        return m_nodeTags;
    }

    std::vector<double> getJacobians()
    {
        return m_jacobian;
    }

    void frontierAndNeighbouring(std::string gaussType);

};

class Frontier : public Element{

    protected:

    std::vector<double> m_normals;

    public:

    Frontier(const std::string gaussType, int type, int entityTag);
    
    void getNormals(std::vector<double> mainJacobian,
                    int mainNumGp, 
                    std::vector<std::pair<int,int>> neighbours);

    std::vector<double> normals()
    {
        return m_normals;
    }

};

#endif