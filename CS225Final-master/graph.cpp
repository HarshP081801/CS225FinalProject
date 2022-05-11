#include "graph.h"
#include <iostream>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <complex>
#include <iomanip>
#include <numeric>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define INF 0x3f3f3f
#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;

const Vertex Graph::InvalidVertex = "_CS225INVALIDVERTEX";
const int Graph::InvalidWeight = INT_MIN;
const string Graph:: InvalidLabel = "_CS225INVALIDLABEL";
const Edge Graph::InvalidEdge = Edge(Graph::InvalidVertex, Graph::InvalidVertex, Graph::InvalidWeight, Graph::InvalidLabel);

Graph::Graph(bool weighted) : weighted(weighted),directed(false)
{
}

Graph::Graph(bool weighted, bool directed) : weighted(weighted),directed(directed)
{
}


vector<Vertex> Graph::getAdjacent(Vertex source) const 
{
    auto lookup = adjacency_list.find(source);

    if(lookup == adjacency_list.end())
        return vector<Vertex>();

    else
    {
        vector<Vertex> vertex_list;
        unordered_map <Vertex, Edge> & map = adjacency_list[source];
        for (auto it = map.begin(); it != map.end(); it++)
        {
            vertex_list.push_back(it->first);
        }
        return vertex_list;
    }
}


Vertex Graph::getStartingVertex() const
{
    return adjacency_list.begin()->first;
}

vector<Vertex> Graph::getVertices() const
{
    vector<Vertex> ret;

    for(auto it = adjacency_list.begin(); it != adjacency_list.end(); it++)
    {
        ret.push_back(it->first);
    }

    return ret;
}

Edge Graph::getEdge(Vertex source , Vertex destination) const
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return Edge();
    Edge ret = adjacency_list[source][destination];
    return ret;
}

vector<Edge> Graph::getEdges() const
{
    if (adjacency_list.empty())
        return vector<Edge>();

    vector<Edge> ret;
    set<pair<Vertex, Vertex>> seen;

    for (auto it = adjacency_list.begin(); it != adjacency_list.end(); it++)
    {
        Vertex source = it->first;
        for (auto its = adjacency_list[source].begin(); its != adjacency_list[source].end(); its++)
        {
            Vertex destination = its->first;
            if(seen.find(make_pair(source, destination)) == seen.end())
            {
                //this pair is never added to seen
                ret.push_back(its->second);
                seen.insert(make_pair(source,destination));
                if(!directed)
                {
                    seen.insert(make_pair(destination, source));
                }
            }
        }
    }

    return ret;
}

bool Graph::vertexExists(Vertex v) const
{
    return assertVertexExists(v, "");
}

bool Graph::edgeExists(Vertex source, Vertex destination) const
{
    return assertEdgeExists(source, destination, "");
}

Edge Graph::setEdgeLabel(Vertex source, Vertex destination, string label)
{
    if (assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    Edge new_edge(source, destination, e.getWeight(), label);
    adjacency_list[source][destination] = new_edge;

    if(!directed)
    {
        Edge new_edge_reverse(destination,source, e.getWeight(), label);
        adjacency_list[destination][source] = new_edge_reverse;
    }
    return new_edge;
}


string Graph::getEdgeLabel(Vertex source, Vertex destination) const
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidLabel;
    return adjacency_list[source][destination].getLabel();
}

int Graph::getEdgeWeight(Vertex source, Vertex destination) const
{
    if (!weighted)
        error("can't get edge weights on non-weighted graphs!");

    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidWeight;
    return adjacency_list[source][destination].getWeight();
}

void Graph::insertVertex(Vertex v)
{
    // will overwrite if old stuff was there
    removeVertex(v);
    // make it empty again
    adjacency_list[v] = unordered_map<Vertex, Edge>();
}


Vertex Graph::removeVertex(Vertex v)
{

    if (adjacency_list.find(v) != adjacency_list.end())
    {
        if(!directed){
            for (auto it = adjacency_list[v].begin(); it != adjacency_list[v].end(); it++)
            {
                Vertex u = it->first;
                adjacency_list[u].erase(v); 
            }
            adjacency_list.erase(v);
            return v;
        }
        
        adjacency_list.erase(v);
        for(auto it2 = adjacency_list.begin(); it2 != adjacency_list.end(); it2++)
        {
            Vertex u = it2->first;
            if (it2->second.find(v)!=it2->second.end())
            {
                it2->second.erase(v);
            }
        }
        return v;
    }

    return InvalidVertex;
}

bool Graph::insertEdge(Vertex source, Vertex destination)
{
    if(adjacency_list.find(source)!= adjacency_list.end() 
    && adjacency_list[source].find(destination)!= adjacency_list[source].end())
    {
        //edge already exit
        return false;
    }

    if(adjacency_list.find(source)==adjacency_list.end())
    {
        adjacency_list[source] = unordered_map<Vertex, Edge>();
    }
        //source vertex exists
    adjacency_list[source][destination] = Edge(source, destination);
    if(adjacency_list.find(destination)== adjacency_list.end())
    {
        adjacency_list[destination] = unordered_map<Vertex, Edge>();
    }
    if(!directed)
    {
        adjacency_list[destination][source] = Edge(source, destination);
    }
    return true;
}

Edge Graph::removeEdge(Vertex source, Vertex destination)
{
    if(assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    adjacency_list[source].erase(destination);
    // if undirected, remove the corresponding edge
    if(!directed)
    {
        adjacency_list[destination].erase(source);
    }
    return e;
}


Edge Graph::setEdgeWeight(Vertex source, Vertex destination, int weight)
{
    if (assertEdgeExists(source, destination, __func__) == false)
        return InvalidEdge;
    Edge e = adjacency_list[source][destination];
    //std::cout << "setting weight: " << weight << std::endl;
    Edge new_edge(source, destination, weight, e.getLabel());
    adjacency_list[source][destination] = new_edge;

    if(!directed)
        {
            Edge new_edge_reverse(destination,source, weight, e.getLabel());
            adjacency_list[destination][source] = new_edge_reverse;
        }

    return new_edge;
}

bool Graph::assertVertexExists(Vertex v, string functionName) const
{
    if (adjacency_list.find(v) == adjacency_list.end())
    {
        if (functionName != "")
            error(functionName + " called on nonexistent vertices");
        return false;
    }
    return true;
}

bool Graph::assertEdgeExists(Vertex source, Vertex destination, string functionName) const
{
    if(assertVertexExists(source,functionName) == false)
        return false;
    if(adjacency_list[source].find(destination)== adjacency_list[source].end())
    {
        if (functionName != "")
            error(functionName + " called on nonexistent edge " + source + " -> " + destination);
        return false;
    }

    if(!directed)
    {
        if (assertVertexExists(destination,functionName) == false)
            return false;
        if(adjacency_list[destination].find(source)== adjacency_list[destination].end())
        {
            if (functionName != "")
                error(functionName + " called on nonexistent edge " + destination + " -> " + source);
            return false;
        }
    }
    return true;
}

bool Graph::isDirected() const
{
    return directed;
}

void Graph::clear()
{
    adjacency_list.clear();
}


/**
 * Prints a graph error and quits the program.
 * The program is exited with a segfault to provide a stack trace.
 * @param message - the error message that is printed
 */
void Graph::error(string message) const
{
    cerr << "\033[1;31m[Graph Error]\033[0m " + message << endl;
}

/**
 * Prints the graph to stdout.
 */
void Graph::print(string fileName) const
{
    ofstream output;
    output.open(fileName);
    for (auto it = adjacency_list.begin(); it != adjacency_list.end(); ++it) 
    {
        output << it->first << endl;
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) 
        {
            std::stringstream ss;
            ss << it2->first; 
            string vertexColumn = "    => " + ss.str();
            vertexColumn += " " ;
            output << std::left << std::setw(26) << vertexColumn;
            string edgeColumn = "edge label = \"" + it2->second.getLabel()+ "\"";
            output << std::left << std::setw(26) << edgeColumn;
            if (weighted)
                output << "weight = " << it2->second.getWeight();
            output << endl;
        }
        output << endl;
    }
}

void Graph::print() const
{
    for (auto it = adjacency_list.begin(); it != adjacency_list.end(); ++it) 
    {
        cout << it->first << endl;
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) 
        {
            std::stringstream ss;
            ss << it2->first; 
            string vertexColumn = "    => " + ss.str();
            vertexColumn += " " ;
            cout << std::left << std::setw(26) << vertexColumn;
            string edgeColumn = "edge label = \"" + it2->second.getLabel()+ "\"";
            cout << std::left << std::setw(26) << edgeColumn;
            if (weighted)
                cout << "weight = " << it2->second.getWeight();
            cout << endl;
        }
        cout << endl;
    }
}

void Graph::fillGraph(const string& Verticies, const string& Edges)
{
    std::ifstream wikiPages(Verticies);
    Vertex vertex;
    vector<Vertex> verticies;

    //push all the pages into the graph verticies vector and another vector to keep its order 
    while (getline(wikiPages, vertex))
    {
        this->insertVertex(vertex);
        verticies.push_back(vertex);
    }
    this->numVertices = verticies.size();
    wikiPages.close(); //close the pages file
    
    //tiMe tO dO eDgEs
    std::ifstream edgeWeights(Edges);
    string weight;
    vector<int> row; //the vector to store a given vertex's edges
    for (int v = 0; v < numVertices; v++)
    {
        getline(edgeWeights, weight); 
        pair<Vertex, int> key = make_pair(verticies[v], v); //marks the vertex's index in the matrix
        row.clear();
        for(int d = 0; d < numVertices; d++)
        {   
            char w = weight.at(d);
            if( w != '0' && w != '_')
            {
                //make the edge and set it's weight
                this->insertEdge(verticies[v], verticies[d]);
                this->setEdgeWeight(verticies[v], verticies[d], int(weight.at(d))-48);
            }

            if(w == '_')
            {
                
                row.push_back(0);
            }
            else
            {
                //pushes the int version of the char
                row.push_back(int(w)-48);
            }
           
        }
        //storing the matrix data
        edgeMatrixKey.push_back(key);
        edgeMatrix.push_back(row);

    }
    edgeWeights.close();
    
}

void Graph::dijkstra(Vertex source, string fileName)
{
    //each vertex's distance from the source vertex  
    vector<int> dist;
    //keeps track of whether a node is in the tree
    vector<bool> sptSet; 
    //shortest path tree
    vector<int> parent; 

    //Initializing vectors
    dist.resize(numVertices, 0);
    sptSet.resize(numVertices, false);
    parent.resize(numVertices, 0);

    //Finding the matrix index of the source vertex
    int src;
    for (int i = 0; i < numVertices; i++) 
    { 
        parent[0] = -1; 
        dist[i] = INT_MAX; 
        sptSet[i] = false; 

        if(edgeMatrixKey[i].first == source)
        {
            src = edgeMatrixKey[i].second;
        }
    } 
  
    
    //initializing and updating the shortest distance
    dist[src] = 0; 
    for (int count = 0; count < numVertices - 1; count++) 
    {  
        //finding the min distance to the vertex and marking that its in the parent tree
        int u = minDistance(dist, sptSet); 
        sptSet[u] = true; 

        for (int v = 0; v < numVertices; v++) 
        {
            if (!sptSet[v] && edgeMatrix[u][v] && 
                dist[u] + edgeMatrix[u][v] < dist[v]) 
            { 
                parent[v] = u; 
                dist[v] = dist[u] + edgeMatrix[u][v]; 
            }  
        }
    } 

    //Outputting the results to a .txt file
    ofstream output;
    output.open(fileName, std::ios_base::app);

    output << "Vertex->Label" << endl;
    for(int i = 0; i < numVertices; i++)
    {
        output << edgeMatrixKey[i].first << "->" << to_string(edgeMatrixKey[i].second) << endl;
    }

    output << "\t" << endl;

    output << "Vertex\t Distance" << endl; 
    for (int i = 0; i < numVertices; ++i)
    {
        output << to_string(src) << "-> " << to_string(i) << "\t\t" << to_string(dist[i]) << endl;
    }

}

int Graph::minDistance(vector<int> dist,  vector<bool> spt) 
{ 
    // Initialize min value 
    int min = INT_MAX;
    int min_index; 

    for (int v = 0; v < numVertices; v++) 
        if (spt[v] == false && dist[v] <= min) 
        {
            min = dist[v];
            min_index = v;
        } 
  
    return min_index; 
} 


void Graph::pageRank(string fileName)
{
    vector<double> oldPR(numVertices), newPR(numVertices);
    double p = 1.0 / double(numVertices);

    //Assign equal rank to every page to begin
    for (int i = 0; i < numVertices; i++) 
    {
        oldPR[i] = double(1.0 / numVertices);
        newPR[i] = oldPR[i];
    }
    
    int numIt = 0;
    double d = 0.95; //decay factor (how willing people are to keep clicking)
    bool close = false;

    while (!close)
    {
        oldPR = newPR;

        //set the new rank to 0
        newPR.clear();
        newPR.resize(numVertices, 0.0);
        double increase = 0;

        //equally split page rank between all outgoing neighbors
        for (int v = 0; v < numVertices; v++)
        {
            //find number of neighbors
            int uniqLinks = 0;
            for(int i = 0; i < numVertices; i++)
            {
                if(edgeMatrix[v][i] != 0)
                {
                    uniqLinks++;
                }
            }

            if(uniqLinks != 0) // if it has incoming neighbors
            {
                increase = (oldPR[v] / uniqLinks) * d;
                for(int n = 0; n < numVertices; n++)
                {
                    if(edgeMatrix[v][n] != 0)
                    {
                        newPR[n] = newPR[n] + increase;
                    }
                }
            }
            else //no incoming links
            {
                increase = (oldPR[v] / numVertices) * d;
                for(int n = 0; n < numVertices; n++)
                {
                    newPR[n] = newPR[n] + increase;
                }
            }

            //accounting for decay
            increase = (1.0-d)/double(numVertices);
            for(int i = 0; i < numVertices; i++)
            {
                newPR[i] = newPR[i] + increase;
            } 

            //testing whether the values have converged
            close = isClose(newPR, oldPR, 0.00001);
        }
    }
    
    //printing to output file
    ofstream output;
    output.open(fileName, std::ios_base::app);
    output << "Vertex\t Rank" << endl; 
    for (int i = 0; i < numVertices; ++i)
    {
        output << to_string(numIt) << "-> " << to_string(i) << "\t\t" << to_string(newPR[i]) << "\t\t" << endl;
    }
}


bool Graph::isClose(vector<double> newPR, vector<double> oldPR, double e)
{
    for (int i = 0; i < numVertices; i++)
    {
        double diff = newPR[i] - oldPR[i];
        if (diff < 0)
        {
            diff = diff * -1;
        }

        if (diff > e)
        {
            return false;
        }
    }

    return true;
}
