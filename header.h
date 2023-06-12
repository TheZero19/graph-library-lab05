#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/degree_centrality.hpp>
#include <boost/graph/clustering_coefficient.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;

//function prototype:
int calculateDiameter(Graph graph);
double calculateAverageClusteringCoefficientIdontknowWhyitFails(Graph graph); //I tried to implement this on my own but it failed miserably working for only few graphs
double calculateAverageClusteringCoefficient(const Graph& graph); //I used the library's implementation instead for clustering coefficient
void calculateDegreeDistribution(const Graph& graph, std::unordered_map<int, int>& degreeDistribution);
void printDegreeDistribution(const std::unordered_map<int, int>& degreeDistribution, int numNodes);

void calculations(Graph graph){
    // Calculate the number of nodes and edges:
    int numberOfNodes = num_vertices(graph);
    std::cout<<"Number of Nodes: " << numberOfNodes << std::endl;
    int numberOfEdges = num_edges(graph);
    std::cout<<"Number of edges: " << numberOfEdges << std::endl;

    // Calculate the average degree
    double totalDegree = 0;
    auto vertexRange = vertices(graph);
    for (auto vertexIt = vertexRange.first; vertexIt != vertexRange.second; ++vertexIt) {
        totalDegree += degree(*vertexIt, graph);
    }
    double averageDegree = totalDegree / num_vertices(graph);
    std::cout << "Average Degree: " << averageDegree <<std::endl;

    //Calculate density
    double density;
    if(numberOfNodes == 0){
        density = 0;
    } 
    else{
        density = (2*(double)numberOfEdges)/(numberOfNodes*(numberOfNodes - 1));
    } 
    std::cout << "Density: " << density << std::endl;
    
    std::cout << "Diameter: " << calculateDiameter(graph) << std::endl;

    //clustering coefficient
    double averageClusteringCoefficient = calculateAverageClusteringCoefficient(graph);
    cout << "Clustering Coefficient: " << averageClusteringCoefficient << endl;

    // Calculate the degree distribution
    std::unordered_map<int, int> degreeDistribution;
    calculateDegreeDistribution(graph, degreeDistribution);
    printDegreeDistribution(degreeDistribution, numberOfNodes);

    
}


int calculateDiameter(Graph graph){
    int diameter = 0;

    auto vertexRange = vertices(graph);
    for (auto vertexIt = vertexRange.first; vertexIt != vertexRange.second; ++vertexIt) {
        int maxDistance = 0;

        // Perform a BFS starting from the current vertex
        vector<bool> visited(num_vertices(graph), false);
        vector<int> distance(num_vertices(graph), 0);
        queue<int> q;

        visited[*vertexIt] = true;
        q.push(*vertexIt);

        while (!q.empty()) {
            int currentVertex = q.front();
            q.pop();

            // Update the maximum distance reached
            maxDistance = distance[currentVertex];

            // Visit all neighboring vertices
            auto neighborRange = adjacent_vertices(currentVertex, graph);
            for (auto neighborIt = neighborRange.first; neighborIt != neighborRange.second; ++neighborIt) {
                if (!visited[*neighborIt]) {
                    visited[*neighborIt] = true;
                    q.push(*neighborIt);
                    distance[*neighborIt] = distance[currentVertex] + 1;
                }
            }
        }

        // Update the diameter of the graph
        diameter = max(diameter, maxDistance);

    }
    return diameter;
}

double calculateAverageClusteringCoefficientIdontknowWhyitFails(Graph graph) {
    double sumClusteringCoefficient = 0.0;

    auto vertexRange = vertices(graph);
    for (auto vertexIt = vertexRange.first; vertexIt != vertexRange.second; ++vertexIt) {
        int ki = degree(*vertexIt, graph);
        if (ki <= 1) {
            continue;
        }

        int ei = 0;
        auto neighborRange = adjacent_vertices(*vertexIt, graph);
        for (auto neighborIt = neighborRange.first; neighborIt != neighborRange.second; ++neighborIt) {
            auto neighbor = *neighborIt;
            auto neighborNeighborRange = adjacent_vertices(neighbor, graph);
            for (auto neighborNeighborIt = neighborNeighborRange.first; neighborNeighborIt != neighborNeighborRange.second; ++neighborNeighborIt) {
                auto neighborNeighbor = *neighborNeighborIt;
                if (neighborNeighbor != *vertexIt && edge(neighborNeighbor, neighbor, graph).second) {
                    ++ei;
                }
            }
        }

        double clusteringCoefficient = 2.0 * ei / (ki * (ki - 1));
        sumClusteringCoefficient += clusteringCoefficient;
    }

    double averageClusteringCoefficient = sumClusteringCoefficient / num_vertices(graph);
    return averageClusteringCoefficient;
}


void calculateDegreeDistribution(const Graph& graph, std::unordered_map<int, int>& degreeDistribution) {
    auto vertexRange = vertices(graph);
    for (auto vertexIt = vertexRange.first; vertexIt != vertexRange.second; ++vertexIt) {
        int degree = boost::out_degree(*vertexIt, graph);
        ++degreeDistribution[degree];
    }
}

void printDegreeDistribution(const std::unordered_map<int, int>& degreeDistribution, int numNodes) {
    for (const auto& entry : degreeDistribution) {
        int degree = entry.first;
        int count = entry.second;
        double fraction = static_cast<double>(count) / numNodes;
        cout << "Degree " << degree << ": " << fraction << endl;
    }
}


double calculateAverageClusteringCoefficient(const Graph& graph) {
    typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
    std::vector<double> clusteringCoefficients(num_vertices(graph));

    VertexIterator vertexIt, vertexEnd;
    boost::tie(vertexIt, vertexEnd) = vertices(graph);

    for (; vertexIt != vertexEnd; ++vertexIt) {
        clusteringCoefficients[*vertexIt] = boost::clustering_coefficient(graph, *vertexIt);
    }

    double sumClusteringCoefficient = 0.0;
    for (const auto& cc : clusteringCoefficients) {
        sumClusteringCoefficient += cc;
    }

    double averageClusteringCoefficient = sumClusteringCoefficient / num_vertices(graph);
    return averageClusteringCoefficient;
}
