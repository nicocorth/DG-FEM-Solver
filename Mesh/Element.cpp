#include "mesh.hpp"
#include <Eigen/Eigen>
#include <gmsh.h>
#include <iostream>
#include <string.h>

// This is the constructor of an element. It gathers all required data for the dgfem resolution.
Element::Element(const std::string gaussType, int type, int entityTag)
{

    std::vector<double> bin;
    std::size_t i, j, k, l;

    m_type = type;

    std::vector<std::size_t> tmp_sortedNodes(m_frontierNodes.size());// temporary storage of the sorted nodes.
    std::vector<std::pair<std::size_t, std::size_t>> tmp_neighbours(m_tag.size() * m_numFrontierNodes); // Vector storing the neighbours temporarily.
    std::vector<std::size_t> sortedNodes; // Final vector of sorted nodes.

    // Get the global properties of the given element type.
    gmsh::logger::write("Element properties...");
    gmsh::model::mesh::getElementProperties(m_type, 
                                            m_name,
                                            m_dim,
                                            m_order,
                                            m_numNodes,
                                            bin);

    gmsh::logger::write("Element by Type...");
    // Gets the tags of the elements and the tags of their nodes
    gmsh::model::mesh::getElementsByType(m_type,
                                         m_tag,
                                         m_nodeTags,
                                         entityTag);

    gmsh::logger::write("Integration points...");
    // Getting the integration points of the element.
    gmsh::model::mesh::getIntegrationPoints(m_type,
                                            gaussType,
                                            m_gaussPoints,
                                            m_gaussWeight);

    // Gauss points number of the element.
    m_gaussPointsNumber = m_gaussPoints.size()/3;
    m_totalGaussPointsNumber = m_gaussPointsNumber * m_tag.size();
    
    gmsh::logger::write("Jacobians...");
    // Getting the jacobian of all elements at the Gauss points.
    gmsh::model::mesh::getJacobians(m_type,
                                    m_gaussPoints,
                                    m_jacobianMatrix,
                                    m_jacobian, 
                                    m_gaussCoordinates,
                                    entityTag);
 

    // Inverting the jacobian matrices
    gmsh::logger::write("Jacobians Inversion...");
    std::vector<Eigen::MatrixXd> jacobianMatrix(m_jacobianMatrix.size()/9);

    m_jacobianInverse.resize(m_jacobianMatrix.size());

    for(i = 0; i < jacobianMatrix.size(); ++i)
    {
        jacobianMatrix[i].setZero(3,3);
    
        for(j = 0; j < 3; ++j)
        {
            for(k = 0; k < 3; ++k)
            {
                jacobianMatrix[i](j, k) = m_jacobianMatrix[i * 9 + j * 3 + k];
            }
        }
    }

    for(i = 0; i < jacobianMatrix.size(); ++i)
    {
        Eigen::MatrixXd inverse = jacobianMatrix[i].inverse();

        for(j = 0; j < 3; ++j)
        {
            for(k = 0; k < 3; ++k)
            {
                m_jacobianInverse[i * 9 + j * 3 + k] = inverse(j, k);
            }
        }
    }

    gmsh::logger::write("Basis functions...");
    // Gets the basis functions of the element at the Gauss points.
    gmsh::model::mesh::getBasisFunctions(m_type,
                                         m_gaussPoints,
                                         "Lagrange",
                                         m_basisFunctionsCompo,
                                         m_basisFunctions);

    // Number of basis functions in the element.
    m_basisFunctionsNumber = m_basisFunctions.size() / (m_basisFunctionsCompo * m_gaussPointsNumber);


    gmsh::logger::write("Gradient Basis Functions...");
    // Gets the gradient of the basis functions of the element at the Gauss points.
    gmsh::model::mesh::getBasisFunctions(m_type,
                                         m_gaussPoints,
                                         "GradLagrange",
                                         m_gradBasisFunctionsCompo,
                                         m_gradBasisFunctions); 
 
    // Computation of the real components of the gradient of the basis functions.

    m_gradRealBasisFunctionsX.resize(m_gradBasisFunctions.size() * m_tag.size()/3, 0);
    m_gradRealBasisFunctionsY.resize(m_gradBasisFunctions.size() * m_tag.size()/3, 0);
    m_gradRealBasisFunctionsZ.resize(m_gradBasisFunctions.size() * m_tag.size()/3, 0);

    gmsh::logger::write("Real gradient of Basis Functions...");
    for(i = 0; i < m_tag.size(); ++i) // loop over the element tags.
        for(j = 0; j < m_gaussPointsNumber; ++j) // loop over the gauss points.
            for(k = 0; k < m_numNodes; ++k) // loop over the nodes.
            {
                int realIndex = i * m_gaussPointsNumber * m_numNodes + j * m_numNodes + k;

                for(l = 0; l < 3; ++l) // loop over the lines of the jacobian.
                {   
                    int jacobIndexX = i * m_gaussPointsNumber * 9 + j * 9 + l;
                    int jacobIndexY = i * m_gaussPointsNumber * 9 + j * 9 + 3 + l;
                    int jacobIndexZ = i * m_gaussPointsNumber * 9 + j * 9 + 6 + l;
                    int paramIndex = j * m_numNodes * 3 + k * 3 + l;

                    // Multiplies the corresponding row of the jacobian with the isoparam gradient to obtain the real gradient.
                    m_gradRealBasisFunctionsX[realIndex] += m_jacobianInverse[jacobIndexX]
                                                         * m_gradBasisFunctions[paramIndex];
                    m_gradRealBasisFunctionsY[realIndex] += m_jacobianInverse[jacobIndexY]
                                                         * m_gradBasisFunctions[paramIndex];
                    m_gradRealBasisFunctionsZ[realIndex] += m_jacobianInverse[jacobIndexZ]
                                                         * m_gradBasisFunctions[paramIndex];
                }
            }

    gmsh::logger::write("Element frontier nodes...");
    // Get the node at the edge or at the surface of an element.
    if(m_dim < 3)
    {
        // Get all the nodes at the edge.
        gmsh::model::mesh::getElementEdgeNodes(m_type,
                                               m_frontierNodes,
                                               entityTag);
    }

    else
    {
        int faceType;
        if(m_name.find("Tetrahedron") != std::string::npos)
        {
            faceType = 3;
        }
        else if (m_name.find("Hexahedron") != std::string::npos)
        {
            faceType = 4;
        }

        gmsh::model::mesh::getElementFaceNodes(m_type, 
                                               faceType,
                                               m_frontierNodes,
                                               entityTag);

    }

    m_numFrontierNodes = m_frontierNodes.size()/m_tag.size();

    // Number of frontier per element type.
    if(m_name.find("Point") != std::string::npos)
    {
        m_numFrontier = 1;
    } 
    else if(m_name.find("Line") != std::string::npos)
    {
        m_numFrontier = 2;
    }
    else if(m_name.find("Triangle") != std::string::npos)
    {
        m_numFrontier = 3;
    }
    else if(m_name.find("Quadrilateral") != std::string::npos 
            || m_name.find("Tetrahedron") != std::string::npos)
    {
        m_numFrontier = 4;
    }

    else if(m_name.find("Hexahedron") != std::string::npos)
    {
        m_numFrontier = 6;
    }
    
    else
    {
        gmsh::logger::write("The element is unsupported by gmsh or by the solver.","error");
    }

    // Number of nodes per frontier
    m_numNodesPerFrontier = m_frontierNodes.size()/(m_tag.size() * m_numFrontier);

    gmsh::logger::write("Matrices computation...");
    // Mass and stiffness matrices computation.
    std::vector<Eigen::MatrixXd> massMatrix(m_tag.size());

    m_massMatrixInverse.resize(m_basisFunctionsNumber * m_basisFunctionsNumber * m_tag.size(), 0);
    m_stiffnessMatrixX.resize(m_basisFunctionsNumber * m_basisFunctionsNumber * m_tag.size(), 0);
    m_stiffnessMatrixY.resize(m_basisFunctionsNumber * m_basisFunctionsNumber * m_tag.size(), 0);
    m_stiffnessMatrixZ.resize(m_basisFunctionsNumber * m_basisFunctionsNumber * m_tag.size(), 0);

    for(i = 0; i < massMatrix.size(); ++i)
    {
        massMatrix[i].setZero(m_basisFunctionsNumber, 
                              m_basisFunctionsNumber);
    }

    for(i = 0; i < m_tag.size(); ++i)
    {

        for(j = 0; j < m_basisFunctionsNumber; ++j)
        {

            for(k = 0; k < m_basisFunctionsNumber; ++k)
            {
                
                for(l = 0; l < m_gaussPointsNumber; ++l)
                {

                    int matrixIndex = i * m_basisFunctionsNumber * m_basisFunctionsNumber
                                    + j * m_basisFunctionsNumber
                                    + k;

                    int realGradIndex = i * m_gaussPointsNumber * m_basisFunctionsNumber 
                                      + l * m_basisFunctionsNumber 
                                      + j;

                    massMatrix[i](j,k) += m_gaussWeight[l] 
                                        * m_basisFunctions[l * m_basisFunctionsNumber + j] 
                                        * m_basisFunctions[l * m_basisFunctionsNumber + k] 
                                        * m_jacobian[i * m_gaussPointsNumber + l];
                
                    m_stiffnessMatrixX[matrixIndex] += m_gaussWeight[l] 
                                                     * m_basisFunctions[l * m_basisFunctionsNumber + k] 
                                                     * m_gradRealBasisFunctionsX[realGradIndex] 
                                                     * m_jacobian[i * m_gaussPointsNumber + l];

                    m_stiffnessMatrixY[matrixIndex] += m_gaussWeight[l] 
                                                     * m_basisFunctions[l * m_basisFunctionsNumber + k] 
                                                     * m_gradRealBasisFunctionsY[realGradIndex] 
                                                     * m_jacobian[i * m_gaussPointsNumber + l];

                    m_stiffnessMatrixZ[matrixIndex] += m_gaussWeight[l] 
                                                     * m_basisFunctions[l * m_basisFunctionsNumber + k] 
                                                     * m_gradRealBasisFunctionsZ[realGradIndex]
                                                     * m_jacobian[i * m_gaussPointsNumber + l];

                }
                
            }

        }

        massMatrix[i] = massMatrix[i].inverse();
        
        // Copy the eigen matrix into a c++ vector.
        for(j = 0; j < m_basisFunctionsNumber; ++j)
        {

            for(k = 0; k < m_basisFunctionsNumber; ++k)
            {
                int matrixIndex = i * m_basisFunctionsNumber * m_basisFunctionsNumber
                                  + j * m_basisFunctionsNumber
                                  + k;

                m_massMatrixInverse[matrixIndex] = massMatrix[i](j,k);  

            }

        }

    }

}

// This fuctions elmiminates the 
void Element::frontierAndNeighbouring(std::string gaussType)
{
    std::size_t i, j, k, l;

    std::vector<std::size_t> tmp_sortedNodes(m_frontierNodes.size());// temporary storage of the sorted nodes.
    std::vector<std::pair<std::size_t, std::size_t>> tmp_neighbours(m_tag.size()* m_numFrontier); // Vector storing the neighbours temporarily.
    std::vector<std::size_t> sortedNodes; // Final vector of sorted nodes.

    std::size_t sortedSize = 0; // Integer that always keeps the position of the bloc in tmp_sortedNodes.
    std::size_t neighboursSize = 0;
    int nodeCount; // Counter of matching nodes bewtenn two blocks

    std::vector<int> frontierType(1);
    std::size_t frontierTag;
    std::string frontierName;

    // Initialisation of tmp_sortedNodes. It first contains the frontier of the first element.
    for(i = 0; i < m_numFrontierNodes; ++i)
    {
        tmp_sortedNodes[sortedSize + i] = m_frontierNodes[i];
    }

    // Initialization of the first entries of neighbours vector.
    for(i = 0; i < m_numFrontier; ++i)
    {
        tmp_neighbours[i].first = 0;
    }

    // Initialize everything to 0 for the second neighbour.
    for(i = 0; i < tmp_neighbours.size(); ++i)
    {
        tmp_neighbours[i].second = -1;
    }

    sortedSize += m_numFrontierNodes;
    neighboursSize += m_numFrontier;

    // Construction of the sorted nodes and the neighbours.
    for(i = m_numFrontierNodes; i < m_frontierNodes.size(); i += m_numNodesPerFrontier)
    {
        nodeCount = 0;

        for(j = 0; j < sortedSize && nodeCount != m_numNodesPerFrontier; j += m_numNodesPerFrontier)
        {
            nodeCount = 0;

            for(k = 0; k < m_numNodesPerFrontier; ++k)
            {
                for(l = 0; l < m_numNodesPerFrontier; ++l)
                {
                    if(m_frontierNodes[i + k] == tmp_sortedNodes[j + l])
                    {
                        ++nodeCount;
                    }
                }
            }
        }

        if(nodeCount != m_numNodesPerFrontier)
        {
            for(k = 0; k < m_numNodesPerFrontier; ++k)
            {
                tmp_sortedNodes[sortedSize + k] = m_frontierNodes[i + k];
            }

            tmp_neighbours[neighboursSize].first = i/m_numFrontierNodes;
            tmp_neighbours[neighboursSize].second = -1;

            sortedSize += m_numNodesPerFrontier;
            ++neighboursSize;

        }

        else if(m_tag[i/m_numFrontierNodes] != tmp_neighbours[j/m_numNodesPerFrontier].first)
        {
            tmp_neighbours[j/m_numNodesPerFrontier - 1].second = i/m_numFrontierNodes;
        }
        
    }

    m_neighbours.resize(neighboursSize);
    sortedNodes.resize(sortedSize);

    for(i = 0; i < neighboursSize; ++i)
    {
       m_neighbours[i] = tmp_neighbours[i];
    }

    for(i = 0; i < sortedSize; ++i)
    {
        sortedNodes[i] = tmp_sortedNodes[i];
    }

    m_nodeCorrespondance.resize(sortedNodes.size());

    for(i = 0; i < sortedNodes.size(); ++i)
    {
        for(j = 0; j < m_numNodes; ++j)
        {
            
            if(m_nodeTags[m_neighbours[i/m_numNodesPerFrontier].first * m_numNodes + j] == sortedNodes[i])
            {
                m_nodeCorrespondance[i].first = m_neighbours[i/m_numNodesPerFrontier].first * m_numNodes + j;
            }

            if(m_neighbours[i/m_numNodesPerFrontier].second >= 0)
            {
                if(m_nodeTags[m_neighbours[i/m_numNodesPerFrontier].second * m_numNodes + j] == sortedNodes[i])
                {
                    m_nodeCorrespondance[i].second = m_neighbours[i/m_numNodesPerFrontier].second * m_numNodes + j;
                }
            }

            else
            {
                m_nodeCorrespondance[i].second = -1;
            }
            
        }
    }
    
    if(gmsh::model::getDimension() == 3)
    {
        if(m_name.find("Tetrahedron") != std::string::npos)
        {

            frontierName = "triangle";

        }

        else if(m_name.find("Hexahedron") != std::string::npos)
        {

            frontierName = "quadrangle";

        }
    }

    else if(gmsh::model::getDimension() == 2)
    {

        frontierName = "line";

    }

    else
    {

        frontierName = "point";

    }

    frontierTag = gmsh::model::addDiscreteEntity(gmsh::model::getDimension() - 1);

    frontierType[0] = gmsh::model::mesh::getElementType(frontierName, m_order);
    
    gmsh::model::mesh::addElementsByType(frontierTag, 
                                         frontierType[0],
                                         {},
                                         sortedNodes);


    m_frontierElements = new Frontier(gaussType, frontierType[0], frontierTag);

    m_frontierElements->getNormals(m_jacobianInverse, 
                                   m_gaussPointsNumber,
                                   m_neighbours);

}

Frontier::Frontier(std::string gaussType, int type, int entityTag) 
    : Element(gaussType, type, entityTag)
{

}

void Frontier::getNormals(const std::vector<double> & mainJacobianInverse,
                          int mainNumGp, 
                          const std::vector<std::pair<int,int>> & neighbours)
{
    std::size_t i, j, k;

    m_normals.resize(m_tag.size() * m_gaussPointsNumber * 3);

    for(i = 0; i < m_tag.size(); ++i) // Run through the elements
    {
        for(j = 0; j < m_gaussPointsNumber; ++j) // run through the gauss points of a given element.
        {
            
            int jacobianIndex = i *  m_gaussPointsNumber * 9 + j * 9;

            int frontierIndex = i * m_gaussPointsNumber * 3 + j * 3;

            int mainJacobIndex = neighbours[i].first * mainNumGp * 9; 

            double compoX = 0, compoY = 0, compoZ = 0;

            double sign = 1;

            double norm;

            if(m_dim == 0)
            {
                compoX = 1;

                if(!i && !j)
                {
                    sign = -1;
                }
                     
            }

            if(m_dim == 1)
            {
                compoX = m_jacobianInverse[jacobianIndex + 1];
                compoY = m_jacobianInverse[jacobianIndex + 4];

                sign = mainJacobianInverse[mainJacobIndex + 8];
            }

            else if(m_dim == 2)
            {
                compoX = m_jacobianInverse[jacobianIndex + 2];
                compoY = m_jacobianInverse[jacobianIndex + 5];
                compoZ = m_jacobianInverse[jacobianIndex + 8];
            }

            norm = sqrt(compoX * compoX + compoY * compoY + compoZ * compoZ);
            
            m_normals[frontierIndex] = sign * compoX/norm;
            m_normals[frontierIndex + 1] = sign * compoY/norm;
            m_normals[frontierIndex + 2] = sign * compoZ/norm;

        }

    }

}

Element::~Element()
{
    if(m_frontierElements)
    {

        delete m_frontierElements;

    }
}