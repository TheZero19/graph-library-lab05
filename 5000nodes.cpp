#include "header.h"
#include <string>


void doCalculations(std::string fileName){

    ifstream file;
    file.open(fileName);

    Graph graph;
    
    string line;

    // Parse the edge CSV file and add edges to the graph
    while (getline(file, line)) {
        istringstream iss(line);
        int source, target;
        // Assuming each line contains a pair of source and target node IDs
        if (iss >> source >> target) {
            add_edge(source, target, graph);
        }
    }

    calculations(graph);
    std::cout << std::endl;

}


int main(){

    // Load the edge CSV file
    
    doCalculations("1.mtx");
    doCalculations("2.mtx");
    doCalculations("3.mtx");
    doCalculations("4.mtx");
    doCalculations("5.mtx");

    


    return 0;
}