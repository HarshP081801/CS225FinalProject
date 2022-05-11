# vjh2-vfgray2-cl47-harshbp2
Final Project

Tasks:

1. print the BFS_Traversal of an AVL_Tree
2. print the Dijktstra's Algorithm
3. print the page rank algorithm


   In order to generate the files, first clear any residual data in outputData.
   Then use make && make final && ./final.

   Each algorithm/feature is outputted in txt files. Review below for each specific algorithm.


How to Use:

1. To accomplish BFS_Traversal, we implemented AVL_Trees from lab_avl and added in a function called levelorderTraversal in avltree_given.cpp and a 
 
   function BFS_Traversal in avl_tree.cpp to print the output.

   The algorithm uses a private queue of const Node pointers and proceeds with the standard level order traversal algorithm.

   To generate a BFS_Traversal, a avltree object is created in the main file and one input file and one output file is required as a parameter along 
   
   with the tree object.

   From the input file, we insert each line of the txt file into the tree as a string with a corresponding integer value given in order of insertion.

   Then the algorithm is called on the tree, outputting the result in the output file.



   *In order to change the input/output files, you must manually write the file names based on where they are stored.

   For example our outputBFS.txt file is located in data->outputData, so to output to that, set the string that output takes to "data/outputData/
    
   outputBFS.txt".
   
2. To accomplish Dijkstra's algorithm we wrote a new function for the graph class in graph.cpp
   
   Restrictions on Data for Dijkstra's to run:
   
   -Input data must be in the form of two files: 
      
      -a line separated list of vertices
      
      -a weighted distance matrix in the format:
         
         -One row per article (the "source" of the shortest paths), in the same order as the list of vertices. 
         
         -Each row contains the distances from the source to all articles, also in the order of the vertices list
         
         -The distance is represented as a single digit, with no separators between values. The distance can be no more than 9
         
         -An underscore ("_") is used to indicate that there isn't an edge between the source.
   
   To run Dijkstra's algorithm, an graph object that is weighted and directed is created in main
   
   The function fillGraph is called (in main) on the newly created graph with 2 string arguments: the first should be the name of the txt file that holds the line separated vertices and the second should be the name of the txt file that holds the distance matrix. This will input the data into a graph compatible with Dijkstra's Algorithm.
   
   Dijkstra is run by calling the dijkstra function on the graph object and inputting the name of the Vertex you wish to be the source, and the name of the output file you wish the shortest distances to be written to.
   
3. PageRank was also a new function we added to the graph class
   
   To run page rank on a data set, the data must be in the same format as dijkstra's algorithm and set up the same way (using fillGraph)
   
   Then call the pageRank function on the graph with a single argument: the name of the file you wish the page ranks to be printed to. 
   
         
   
