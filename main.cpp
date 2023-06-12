#include "header.h"

int main() {
    // Load the edge CSV file
    ifstream edgeFile;
    edgeFile.open("light.edges");

    // Create a Boost Graph using adjacency_list
    Graph graph;

    string line;
    
    // Parse the edge CSV file and add edges to the graph
    while (getline(edgeFile, line)) {
        istringstream iss(line);
        int source, target;
        // Assuming each line contains a pair of source and target node IDs
        if (iss >> source >> target) {
            add_edge(source, target, graph);
        }
    }

    //Calculate all the network properties of small graph and print it
    calculations(graph);

    return 0;
}
