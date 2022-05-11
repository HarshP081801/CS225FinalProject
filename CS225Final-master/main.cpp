/**
 * @file main.cpp
 */

#include <iostream>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include "avl_tree.h"
#include "graph.h"
#include <map>


using namespace std;

std::vector<std::string> file_to_vector(const std::string & filename)
{
	vector<string> output;
    ifstream input(filename);
    string str;
    while (getline(input, str))
    {
        output.push_back(str);
    }
    return output;
}

int main() 
{

    //Code for BFS_Traversal
    AVLTree<int, string> tree;
    tree.BFS_Traversal(tree, "data/sampleData/articlesTest.txt", "data/outputData/outputBFS.txt");

    //Code for our data
    Graph* graph = new Graph(true);
    graph->fillGraph("data/sampleData/articlesTest.txt", "data/sampleData/weightsTest.txt");
    graph->print("data/outputData/outputGraph.txt");
    Vertex string = "10th_century";
    graph->dijkstra(string, "data/outputData/outputDijkstra.txt");
    double d = 0.95;
    double e = 0.0001;
    graph->pageRank("data/outputData/outputPageRank.txt");

   /*  Graph* baby = new Graph(true, true);
    baby->fillGraph("data/sampleData/tinyTest.txt", "data/sampleData/tinyMatrix.txt");
    baby->pageRank("data/outputData/outputPageRank.txt"); */
}