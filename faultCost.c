/**
 * faultCost.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE \
"\nUsage: `./faultCost [-m|-l#|-L] [-aco#O#v] [-h]`\n"

#define HELPTEXT \
"Compute fault cost or properties related to number of leaves in a\n\
spanning tree.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout\n\
in graph6 format.\n\
\n\
If -m, -l# or -L is absent the fault cost will be computed.\n\
\n\
    -a, --all\n\
            send distinct degree sequences of every computed\n\
            ml-subgraph to stderr\n\
    -c, --count-branches\n\
            count the number of branches in every computed ml-subgraph;\n\
            introduces overhead\n\
    -h, --help\n\
            print this help message\n\
    -l#, --k-leaf-guaranteed=#\n\
            send all #-leaf-guaranteed graphs to stdout\n\
    -L, --leaf-guaranteed\n\
            send all leaf-guaranteed graphs to stdout\n\
    -m, --ml-number\n\
            compute the minimum leaf number of the input graphs;\n\
            combine with -o# or -O# to send graphs with certain min\n\
            leaf numbers to stdout\n\
    -o#, --output=#\n\
            combine with no arguments or -m to send graphs with fault\n\
            cost # or ml number # to stdout respectively; can be used\n\
            multiple times to output for more values; can be combined\n\
            with -O#\n\
    -O#, --output-from=#\n\
            combine with no arguments or -m to send graphs with fault\n\
            cost at least # or ml number at least # to stdout\n\
            respectively; can be combined with -o#\n\
            \n\
    -v, --verbose\n\
            send more verbose output to stderr\n"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include "utilities/readGraph6.h"
#include "utilities/bitset.h"
#include "utilities/hamiltonicityMethods.h"

typedef struct {
  int **array;
  size_t used;
  size_t size;
} Array;

typedef struct {
  bitset **array;
  size_t used;
  size_t size;
} BSArray;

struct graph {
    int nv;
    bitset *adjacencyList;
    int ne; // NOT KEPT UP TO DATE
    Array *mlsDegrees;
    Array **vtxArrays;
    BSArray *mlsBitsets;
    BSArray **vtxBitsets;
    int mlBound[MAXBITSETSIZE]; // upper bound for ml in G-i (unless 0)
    int K1ham; // -1 UNKNOWN, 0 NO, 1 YES
    int ham[MAXBITSETSIZE+1]; // -1 NO, 0 UNKNOWN, 1 YES.
                              // Entry n -> G, Entry v -> G - v
};

struct subTree { 
    int nv;
    int ne;
    bitset treeVertices;
    bitset availableVertices;
    bitset deletedVertices;
    bitset leaves;
    bitset *adjacencyList;
};


struct options {
    int leafGuaranteed; 
    bool leafGuaranteedFlag; 
    bool allFlag;
    bool fromFlag;
    bool countFlag;
    bool verboseFlag;
    bool mlFlag;
    int output[MAXBITSETSIZE];
    int numberOfOutputIntegers;
    int outputFrom;
};

struct counters {
    long long unsigned int nBranchesFreq[MAXBITSETSIZE];
    long long unsigned int nBranchesFreqSubgraphs[MAXBITSETSIZE];
    long long unsigned int timesAddedToTree;
    long long unsigned int degSequenceAlreadyPresent;
    long long unsigned int skippedGraphs;
    long long unsigned int mlPruned;
    long long unsigned int intermediateGraphs;
    long long unsigned int spanning;
    long long unsigned int decreaseMl;
    long long unsigned int storedWithoutReason;
    long long unsigned int oneHamiltonian;
    long long unsigned int doFullHamPathCheck;
};

long long unsigned int pruneInCost = 0;
long long unsigned int callCost = 0;


//**********************************************************************
//
//                  Macros for dealing with graphs
//
//**********************************************************************

//  Add one edge. 
#define addEdge(g,i,j) {\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);\
 (g)->ne++;\
}

//  Remove one edge.
#define removeEdge(g,i,j) {\
 removeElement((g)->adjacencyList[i], j);\
 removeElement((g)->adjacencyList[j],i);\
 (g)->ne--;\
}

//**********************************************************************
//
//                     I/O 
//
//**********************************************************************

void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printBitset(bitset set) {
    forEach(el, set) {
        fprintf(stderr, "%d ", el);
    }
    fprintf(stderr, "\n");
}


void printMLSDegrees(struct graph *g) {
    fprintf(stderr, "G:\n");
    for(size_t i = 0; i < g->mlsDegrees->used; i++) {
        int *degrees = g->mlsDegrees->array[i];
        fprintf(stderr, "\t");
        for(int j = 0; j < g->nv; j++) {
            fprintf(stderr, "%d ", degrees[j]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Leaves: ");
        printBitset(g->mlsBitsets->array[i][0]);
        fprintf(stderr, "Links: ");
        printBitset(g->mlsBitsets->array[i][1]);
        fprintf(stderr, "Branches: ");
        printBitset(g->mlsBitsets->array[i][2]);
    }


    for(int k = 0; k < g->nv; k++) {
        fprintf(stderr, "G - %d:\n", k);
        for(size_t i = 0; i < g->vtxArrays[k]->used; i++) {
            int *degrees = g->vtxArrays[k]->array[i];
            fprintf(stderr, "\t");
            for(int j = 0; j < g->nv; j++) {
                fprintf(stderr, "%d ", degrees[j]);
            }
            fprintf(stderr, "\n");
        }
    }
}

void countBranchesInMlsubgraph(struct graph *g,
 struct counters *counters) {

    for(size_t i = 0; i < g->mlsDegrees->used; i++) {
        counters->nBranchesFreq[size(g->mlsBitsets->array[i][2])]++;
    }


    for(int k = 0; k < g->nv; k++) {
        for(size_t i = 0; i < g->vtxArrays[k]->used; i++) {

            counters->nBranchesFreqSubgraphs[size(g->vtxBitsets
            [k]->array[i][2])]++;

        }
    }
}

int readGraph(const char *graphString, struct graph *g,
 struct options *options, struct counters *counters) {
    g->nv = getNumberOfVertices(graphString);
    if(g->nv == -1 || g->nv > MAXBITSETSIZE) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }
    g->adjacencyList = malloc(sizeof(bitset)*g->nv);
    if(loadGraph(graphString, g->nv, g->adjacencyList) == -1) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }
    g->K1ham = -1;
    return 0;
}


void createEmptyTree(struct graph *g, struct subTree *T) {
    T->nv = g->nv; 
    T->adjacencyList = calloc(g->nv, sizeof(bitset));
    T->treeVertices = complement(EMPTY, g->nv);
    T->ne = 0;
    T->availableVertices = EMPTY;
    T->leaves = EMPTY;
}


void freeTree(struct subTree *T) {
    free(T->adjacencyList);
}

void freeGraph(struct graph *g) {
    free(g->adjacencyList);
}


//**********************************************************************
//
//                          Dynamic arrays
//
//**********************************************************************



void initArray(Array *a, size_t initialSize, int n) {
  a->array = malloc(initialSize * sizeof(int) * n);
  if(a->array == NULL) {
    fprintf(stderr, "Error: out of memory\n");
    exit(1);
  }
  a->used = 0;
  a->size = initialSize;
}

// Double arraysize when too big.
void insertArray(Array *a, int degrees[], int n) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(int) * n);
    if(a->array == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
  }
  a->array[a->used++] = degrees;
}

void insertArrayAtPos(Array *a, int element[], size_t index) {
    if (index > a->size) {
        fprintf(stderr, "Error: index does not lie in the array.\n");
        exit(1);
    }
  a->array[index] = element;
}

void freeArray(Array *a) {
    for(size_t i = 0; i < a->used; i++) {
        free(a->array[i]);
    }
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void resetArray(Array *a) {
    for(size_t i = 0; i < a->used; i++) {
        free(a->array[i]);
    }
    a->used = 0;
}


void initBitsetArray(BSArray *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(bitset) * 3);
  if(a->array == NULL) {
    fprintf(stderr, "Error: out of memory\n");
    exit(1);
  }
  a->used = 0;
  a->size = initialSize;
}

// Double arraysize when too big.
void insertBitsetArray(BSArray *a, bitset degrees[]) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(bitset) * 3);
    if(a->array == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
  }
  a->array[a->used++] = degrees;
}

void freeBitsetArray(BSArray *a) {
    for(size_t i = 0; i < a->used; i++) {
        free(a->array[i]);
    }
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void resetBSArray(BSArray *a) {
    for(size_t i = 0; i < a->used; i++) {
        free(a->array[i]);
    }
    a->used = 0;
}

//**********************************************************************
//  
//  
//                Definitions for Nauty's splay tree
// 
//  
//**********************************************************************
typedef struct SPLAYNODE {
    int* degSequence;
    struct SPLAYNODE *left, *right, *parent;
} SPLAYNODE;

SPLAYNODE *splayTreeArray[MAXBITSETSIZE + 1] = {NULL};
SPLAYNODE *root = NULL;

#define SCAN_ARGS , int n

int treeSize = 0;

#define ACTION(p)   for(int i = 0; i < n; i++) {\
                        fprintf(stderr, "%d ", p->degSequence[i]);\
                    }\
                    fprintf(stderr, "\n"); treeSize++;

#define INSERT_ARGS , int *degSequence, int n, bool *isPresent

int compareSplayNodeToGraph(SPLAYNODE* p, int *degSequence, int n) {
    return memcmp(p->degSequence, degSequence, sizeof(int) * n);
}

#define COMPARE(p) compareSplayNodeToGraph(p, degSequence, n);

#define PRESENT(p) {(*isPresent) = true;}

#define NOT_PRESENT(p) {p->degSequence = degSequence;\
 (*isPresent) = false;}

#define LOOKUP_ARGS , int *degSequence, int n

// Ignore unused variable warning from splay_insert and splay_lookup
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

    #include "utilities/splay.c" 

#pragma GCC diagnostic pop

void postOrderTraverseFree(SPLAYNODE *p) {

    if(p == NULL) {
        return;
    }

    postOrderTraverseFree(p->left);

    postOrderTraverseFree(p->right);

    free(p);
}

void resetSplayTree(SPLAYNODE **root) {
   postOrderTraverseFree(*root); 
   *root = NULL;
}

//**********************************************************************
//
//              Create degree lists for ML-subgraphs
//
//**********************************************************************

// Free the arrays afterwards!

// If we encounter hamiltonian cycle, create an array of length nv with
// all entries 2.
int* getHamCycleDegrees(int n, bitset deletedVertices) {
    int *degrees = malloc(sizeof(int) * n);
    for(int i = 0; i < n; i++) {
        if(contains(deletedVertices, i)) degrees[i] = -1;
        else degrees[i] = 2;
    }
    return degrees;
}

// If we encounter hamiltonian path, create an array of length nv with
// all entries 2, except for u,v(endpoints), which get 1. 
int* getHamPathDegrees(int n, bitset deletedVertices, int u, int v) {
    int *degrees = malloc(sizeof(int) * n);
    for(int k = 0; k < n; k++) {
        if(contains(deletedVertices, k)) degrees[k] = -1;
        else degrees[k] = 2;
    }
    degrees[u] = 1;
    degrees[v] = 1;
    return degrees;
}

// Given a tree return array of size nv, with each entry being the
// degree of that vertex.
int* getTreeDegrees(struct subTree *T) {
    int *degrees = malloc(sizeof(int) * T->nv);
    for(int k = 0; k < T->nv; k++) {
        if(contains(T->deletedVertices, k)) degrees[k] = -1;
        else degrees[k] = size(T->adjacencyList[k]);
    }
    return degrees;
}

//----------------------------------------------------------------------
//
//              Hamiltonicity
//
//----------------------------------------------------------------------

// Check if hamiltonian, avoid double checking.
bool isHam(struct graph *g, bitset deletedVertices) {

    // If v == -1 we look at G otherwise at G-v
    int v = next(deletedVertices, -1);
    v = ((v == -1) ? g->nv : v);

    if(g->ham[v]) return g->ham[v] == 1; // True if 1, false if -1

    g->ham[v] = isHamiltonian(g->adjacencyList, g->nv,
     deletedVertices, false, false);

    // Now g->ham[v] must be set.
    return g->ham[v] == 1;
}


// Check if K1-hamiltonian i.e. every G-v is ham.
bool isK1Ham(struct graph *g) {

    if(g->K1ham!= -1) return g->K1ham;

    g->K1ham = isK1Hamiltonian(g->adjacencyList, g->nv, false, false,
     -1);

    return g->K1ham;
}

// Check if there is a ham path between s and e in G - excluded.

// Storing every ham path already found and avoiding double checking
// leads to more overhead than gain. Checking for forced endpoints here
// and then pruning if wrong also leads to slightly more overhead than
// gain.
bool hasHamPath(struct graph *g, struct counters *counters, 
 bitset excludedVertices, int s, int e) {
    
    //  If start or end are excluded there cannot be a path between
    //  them.
    if(contains(excludedVertices,s) || contains(excludedVertices,e)) {
        return false;
    }

    counters->doFullHamPathCheck++;    
    
    bitset path = union(singleton(s), singleton(e));
    bitset includedVertices = complement(excludedVertices, g->nv);
    bitset remainingVertices = difference(includedVertices, path);

    //  Will return true if this path can be extended to a hamiltonian
    //  path between start and end and false otherwise..
    return canBeHamiltonian(g->adjacencyList, remainingVertices, s, e,
     size(includedVertices), 2);
}

//----------------------------------------------------------------------
//
//                          Generating trees
//
//----------------------------------------------------------------------


void recursion(struct graph *g, struct subTree *T,
 struct options *options, struct counters *counters, int *ml,
 long long unsigned int *counter,
 long long unsigned int *intermediateGraphs); 


//  Find the vertex v in set for which d_g(v) - d_T(v) is the smallest.
//  Note that forbidden edges have been removed from g here.
int getVertexWithSmallestDegree(struct graph *g, struct subTree *T,
 bitset set) {

    if(isEmpty(set)) {
        return -1;
    }

    int smallest = next(set, -1);
    int smallestDifference = 
     size(intersection(g->adjacencyList[smallest], T->treeVertices)) - 
     size(T->adjacencyList[smallest]);

    forEachAfterIndex(el, set, smallest) { 

        int difference = 
         size(intersection(g->adjacencyList[el], T->treeVertices)) -
         size(T->adjacencyList[el]);
        if(smallestDifference > difference) {
            smallest = el;
            smallestDifference = difference;
        }
    }

    return smallest;
}

//  The next edge e we consider will be incident to a vertex v in which
//  d_g(v) - d_T(v) is smallest and incident to the neighbour w of v
//  for which d_g(w) - d_T(w) is smallest such that e is not yet in
//  T. 
bool getNextEdge(struct graph *g, struct subTree *T, int edge[],
 bool *bothVerticesInTree) {

    edge[0] = getVertexWithSmallestDegree(g, T, T->availableVertices);
    if(edge[0] == -1) {
        return false;
    }

    bitset neighbours = 
     difference(intersection(g->adjacencyList[edge[0]],
      T->treeVertices), T->adjacencyList[edge[0]]);

    edge[1] = getVertexWithSmallestDegree(g, T, neighbours);
    if(edge[1] == -1) {
        return false;
    }

    // We will only remove e if the following is true.
    *bothVerticesInTree = size(T->adjacencyList[edge[1]]);

    return true;
}

// Every time we add an edge to T or remove one from g, we need to
// update which vertices still have incident edges in g not yet in T.
// These are called the availableVertices.
void updateAvailableVertices(struct graph *g, struct subTree *T,
 int edge[]) {

    for(int i = 0; i < 2; i++) {

        int degreeInT = size(T->adjacencyList[edge[i]]);
        int degreeInG = 
         size(intersection(g->adjacencyList[edge[i]], T->treeVertices));

        if(degreeInT > 0 && (degreeInG > degreeInT)) {
            add(T->availableVertices, edge[i]);
            continue;
        }
        removeElement(T->availableVertices, edge[i]);
    }
}

//  Add edge to tree, check if spanning, if not, either prune or
//  continue with recursive algorithm.
//  Creates many arrays for later checking the fault cost. If not needed
//  use simple version of this method;
void addEdgeTree(struct graph *g, struct subTree *T, struct
options *options, struct counters *counters, int *ml, int edge[], long
long unsigned int *counter, long long unsigned
int *intermediateGraphs) {

    bitset origLeaves = T->leaves;
    counters->intermediateGraphs++;

    addEdge(T, edge[0], edge[1]);
    if(size(T->adjacencyList[edge[0]]) == 1) {
        add(T->leaves, edge[0]);
    }
    else {
        removeElement(T->leaves, edge[0]);
    }
    if(size(T->adjacencyList[edge[1]]) == 1) {
        add(T->leaves, edge[1]);
    }
    else {
        removeElement(T->leaves, edge[1]);
    }

    // Can be done earlier by restricting available vertices
    // once #leaves = current ml
    if(size(T->leaves) > *ml) {
        counters->mlPruned++;

        // RESET
        removeEdge(T, edge[0], edge[1]);
        T->leaves = origLeaves;
        return;
    } 

    bitset origAvailable = T->availableVertices;
    updateAvailableVertices(g, T, edge);

    //  Check if spanning
    if(T->ne == size(T->treeVertices) - 1) {
        (*counter)++;
        counters->spanning++;

        // Store data on degree sequences in local variables for ease of
        // use. degrees is dynamic array storing distinct degree
        // sequences of ml-subgraphs of G. bitsets is a dynamic array
        // of length 3 arrays of bitsets where position i gives for the
        // ith ml-subgraph with distinct degree sequence T 1.bitset
        // containing the leaves of T, 2. bitset containing the degree
        // 2 vertices of T (links) and 3. bitset containing the
        // branches of T (deg_T(v) > 2) degSequences same as degrees
        // but stored in splay tree.
        Array *degrees;
        BSArray *bitsets;
        SPLAYNODE **degSequences;

        // ST of original graph
        if(size(T->treeVertices) == g->nv) {
           degrees = g->mlsDegrees;
           bitsets = g->mlsBitsets;
           degSequences = &splayTreeArray[g->nv];
        }
        // ST of G-v
        else {
            int v = next(T->deletedVertices, -1);
            degrees = g->vtxArrays[v]; 
            bitsets = g->vtxBitsets[v];
            degSequences = &splayTreeArray[v];
        }

        int l = size(T->leaves);

        // if found tree with less leaves, reset everything so that we
        // only have data on ml-subgraphs.
        if(*ml > l) {
            if(*ml != INT_MAX) {

                counters->decreaseMl++;
                counters->storedWithoutReason += degrees->used;
            }
            *ml = l;
            resetArray(degrees);
            resetBSArray(bitsets);
            resetSplayTree(degSequences);
        }

        // If we think it is an ml-subgraph add degree sequence data to
        // arrays if the degree sequence is distinct.
        if(*ml == l) {
            int *treeDegrees = getTreeDegrees(T);

            // Check if degree sequence was already in splay tree.
            bool isPresent;
            splay_insert(degSequences, treeDegrees, g->nv,
             &isPresent);
            counters->timesAddedToTree++;
            if(isPresent) {
                free(treeDegrees);
                counters->degSequenceAlreadyPresent++;
            }
            else { //if not insert data
                insertArray(degrees, treeDegrees, g->nv);
                bitset interiorVertices =  // non-leaves of current T
                 complement(union(T->leaves, T->deletedVertices), T->nv);
                bitset links = EMPTY;
                bitset branches = EMPTY;
                forEach(v, interiorVertices) {
                    if(size(T->adjacencyList[v]) == 2) add(links, v);
                    else add(branches, v);
                }
                bitset *degBitsets = malloc(sizeof(bitset) * 3);
                degBitsets[0] = T->leaves;
                degBitsets[1] = links;
                degBitsets[2] = branches;
                insertBitsetArray(bitsets, degBitsets);
            }
        }

        // RESET
        removeEdge(T, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        T->leaves = origLeaves;

        return;
    }

    // If not spanning continue recursion
    recursion(g, T, options, counters, ml, counter, intermediateGraphs);

    //  Reset
    removeEdge(T, edge[0], edge[1]);
    T->availableVertices = origAvailable;
    T->leaves = origLeaves;
}

//  Add edge to tree, check if spanning, if not, either prune or
//  continue with recursive algorithm. No extra storing of data.
void addEdgeTreeSimple(struct graph *g, struct subTree *T, struct
options *options, struct counters *counters, int *ml, int edge[], long
long unsigned int *counter, long long unsigned
int *intermediateGraphs) {

    bitset origLeaves = T->leaves;
    counters->intermediateGraphs++;

    addEdge(T, edge[0], edge[1]);
    if(size(T->adjacencyList[edge[0]]) == 1) {
        add(T->leaves, edge[0]);
    }
    else {
        removeElement(T->leaves, edge[0]);
    }
    if(size(T->adjacencyList[edge[1]]) == 1) {
        add(T->leaves, edge[1]);
    }
    else {
        removeElement(T->leaves, edge[1]);
    }

    if(size(T->leaves) >= *ml) {
        counters->mlPruned++;

        // RESET
        removeEdge(T, edge[0], edge[1]);
        T->leaves = origLeaves;
        return;
    } 

    bitset origAvailable = T->availableVertices;
    updateAvailableVertices(g, T, edge);

    //  Check if spanning
    if(T->ne == size(T->treeVertices) - 1) {
        (*counter)++;
        counters->spanning++;

        int l = size(T->leaves);

        if(*ml > l) {
            if(*ml != INT_MAX) {

                counters->decreaseMl++;
            }
            *ml = l;
        }

        // RESET
        removeEdge(T, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        T->leaves = origLeaves;

        return;
    }

    recursion(g, T, options, counters, ml, counter, intermediateGraphs);

    //  Reset
    removeEdge(T, edge[0], edge[1]);
    T->availableVertices = origAvailable;
    T->leaves = origLeaves;
}

//  Remove edge from graph. No changes are made to T so do not need to
//  check if it is spanning. Check if we can prune, if not, continue
//  with recursive algorithm.
void removeEdgeTree(struct graph *g, struct subTree *T, struct
options *options, struct counters *counters, int *ml, int edge[], long
long unsigned int *counter, long long unsigned
int *intermediateGraphs) {

    removeEdge(g, edge[0], edge[1]);
    bitset origAvailable = T->availableVertices;
    updateAvailableVertices(g, T, edge);

    //  Check if can prune: we can if removing edge created isolated
    //  vertex. Of course, we only need to check the endpoints of the
    //  removed edge.
    int degreeVertex1 = 
     size(intersection(g->adjacencyList[edge[0]], T->treeVertices));
    int degreeVertex2 = 
     size(intersection(g->adjacencyList[edge[1]], T->treeVertices));

    if(degreeVertex1 == 0 || degreeVertex2 == 0) {
        addEdge(g, edge[0], edge[1]);
        T->availableVertices = origAvailable;
        return;
    }

    recursion(g, T, options, counters, ml, counter, intermediateGraphs);

    //  Reset
    addEdge(g, edge[0], edge[1]);
    T->availableVertices = origAvailable;
}

// Recursion for generating all spanning trees. Compute next edge, add
// and forbid it recursively.
void recursion(struct graph *g, struct subTree *T,
 struct options *options, struct counters *counters, int *ml,
 long long unsigned int *counter,
 long long unsigned int *intermediateGraphs) {

    (*intermediateGraphs)++; // Counts number of recursive calls.

    int edge[2];
    bool bothVerticesInTree = false;

    //  Check if there is a next edge and store in edge[]
    if(getNextEdge(g, T, edge, &bothVerticesInTree)) {

        //  If edge makes cycle in T do not add it to T.
        if(!bothVerticesInTree) {
            // If computing fault cost call this version of addEdgeTree
            if(!options->leafGuaranteedFlag && !options->mlFlag) {
                addEdgeTree(g, T, options, counters, ml, edge, counter,
                 intermediateGraphs);
            }
            else { // When not computing fault cost
                addEdgeTreeSimple(g, T, options, counters, ml, edge,
                 counter, intermediateGraphs);
            }
        }

        // Forbid edge (this is done by simply removing it from g)
        removeEdgeTree(g, T, options, counters, ml, edge, counter,
         intermediateGraphs);
    }
}

// Initialize the generation of all ml-subgraph spanning trees. When
// calling we know graph is non-hamiltonian and non-traceable.
long long unsigned int genSpanningTrees(struct graph *g,
 struct options *options, struct counters *counters,
 bitset deletedVertices, int *minLeafNumber) {

    // If we have some upper bound on ml of G-v already set it.
    int ml = INT_MAX;
    forEach(v, deletedVertices) {
        if(g->mlBound[v]) ml = g->mlBound[v];
    }

    struct subTree T;
    createEmptyTree(g, &T); 
    long long unsigned int counter = 0;
    long long unsigned int intermediateGraphs = 0;

    T.deletedVertices = deletedVertices;
    T.treeVertices = difference(T.treeVertices, deletedVertices);
    if(size(T.treeVertices) == 0) {
        freeTree(&T);
        return 0;
    }

    //  Get the first allowed vertex v and denote its neighbours by
    //  u1,...uk. First we start the recursion for all trees with vu1,
    //  then for all trees without vu1 with vu2, then for all trees
    //  without vu1,vu2 and with vu3 etc. 
    int v = next(T.treeVertices, -1);
    bitset nbrs = g->adjacencyList[v]; // Make copy of original nbrs.

    // Check if starting vertex is not isolated.
    if(isEmpty(intersection(nbrs, T.treeVertices))) {
        return 0;
    }

    forEach(u, intersection(nbrs, T.treeVertices)) {

        int edge[2] = {v, u};
        addEdge(&T, v, u);
        updateAvailableVertices(g,&T,edge);
        T.leaves = union(singleton(u), singleton(v));

        if(T.ne == size(T.treeVertices) - 1) {
            freeTree(&T);
            return 1;
        }

        recursion(g, &T, options,
        counters, &ml, &counter, &intermediateGraphs);

        removeEdge(&T, v, u);
        updateAvailableVertices(g,&T,edge);
        T.leaves = EMPTY;
        
        removeEdge(g, v, u);
        updateAvailableVertices(g,&T,edge);
    } 

    //  Restore our original graph
    forEach(nbr, nbrs) {
        addEdge(g, v, nbr);
    }

    // if(options->verboseFlag) {
    //     fprintf(stderr, "Number of spanning trees found: %llu\n",
    //      counter);
    //     fprintf(stderr, "Times recursion step was performed: %llu\n",
    //      intermediateGraphs);
    // }

    if(minLeafNumber != NULL) *minLeafNumber = ml;
    freeTree(&T);
    return counter;
}

//----------------------------------------------------------------------
//
//              Methods for computing fault cost
//
//----------------------------------------------------------------------

// Assumes that degrees is degreelist of full graph
int transitionCost(int n, int degrees1[], int degrees2[], int
upperBound) {

    callCost++;

    int cost = 0;
    for(int i = 0; i < n; i++) {
        if(degrees2[i] == -1) continue;
        if(degrees1[i] != degrees2[i]) cost++;
        if(cost >= upperBound) {
            pruneInCost++;
            return cost;
        }
    }
    return cost;
}

// Find minimum transition cost between the ml-subgraphs of G-v and a
// fixed ml-subgraph of G. vDegrees is Array containing degree arrays
// for each ml-subgraph of G-v, degrees1 contains the degrees of one
// ml-subgraph of G. 
int minimumOverVtxDeleted(int n, int degrees1[], Array *vDegrees,
 int *vBestIdx) {

    if(!vDegrees->used) {
        fprintf(stderr, 
         "Graphs should be 2-connected.\n");
        exit(0);
    }

    int min = INT_MAX;

    for(size_t i=0; i < vDegrees->used; i++) {
        int cost = transitionCost(n, degrees1, vDegrees->array[i], min);
        if(cost < min) {
            min = cost;
            *vBestIdx = i;
        }
    }

    return min;
}

// Find maximum value over minimum transition costs for ml-subgraphs of
// all G-v with a fixed ml-subgraph of G. 
int maxOverVertices(int n, int degrees1[], Array *vtxArrays[],
 int *vBest, int *vBestIdx) {
    int max = -1;
    for(int v = 0; v < n; v++) {
        int vIdx = -1;
        int cost = minimumOverVtxDeleted(n, degrees1, vtxArrays[v],
         &vIdx);
        if(cost > max) {
            max = cost;
            *vBest = v;
            *vBestIdx = vIdx;
        }
    }
    return max;
}

// Same as before but we try to be faster by employing bitsets of leaves
// and links. Seems to be worth the overhead of storing and computing
// these bitsets.
int transitionCost2(struct graph *g, size_t i, bitset leaves,
 bitset links, bitset branches, int v, size_t j, bitset vLeaves, 
 bitset vLinks, bitset vBranches, int lowerBound) {

    callCost++;

    // Neem mask hier
    bitset vMaskInverse = complement(singleton(v), g->nv);

    // Leaves of G st not leaves in G-v ST
    int cost = size(difference(intersection(leaves, vMaskInverse),
     vLeaves));

    if(cost >= lowerBound) {
        pruneInCost++;
        return cost;
    }

    // Links of G st not links in G-v ST
    cost += size(difference(intersection(links, vMaskInverse), vLinks));

    if(cost >= lowerBound) {
        pruneInCost++;
        return cost;
    }

    // Branches of G st not branches in G-v ST
    cost += size(difference(intersection(branches, vMaskInverse),
     vBranches));

    if(cost >= lowerBound) {
        pruneInCost++;
        return cost;
    }

    bitset branchesBoth = intersection(branches, vBranches);
    int *degrees = g->mlsDegrees->array[i];
    int *vDegrees = g->vtxArrays[v]->array[j];

    forEach(u, branchesBoth) {
        if(degrees[u] != vDegrees[u]) cost++;
    }
    return cost;
}

int minOverVtxDeleted2(struct graph *g, size_t i, int v, int *vBestIdx)
{

    bitset leaves = g->mlsBitsets->array[i][0];
    bitset links = g->mlsBitsets->array[i][1];
    bitset branches = g->mlsBitsets->array[i][2];

    BSArray *vBitsets = g->vtxBitsets[v];

    // Maybe do not perform check if this method is bottleneck
    if(!vBitsets->used) {
        fprintf(stderr, 
         "Graphs should be 2-connected.\n");
        exit(0);
    }

    int min = INT_MAX;

    for(size_t j=0; j < vBitsets->used; j++) {
        bitset vLeaves = vBitsets->array[j][0];
        bitset vLinks = vBitsets->array[j][1];
        bitset vBranches = vBitsets->array[j][2];
        int cost = transitionCost2(g, i, leaves, links, branches, v, j, vLeaves, vLinks, vBranches, min);
        if(cost < min) {
            min = cost;
            *vBestIdx = j;
        }
    }

    return min;
}

int maxOverVertices2(struct graph *g, size_t i, int *vBest,
int *vBestIdx) {

    int max = -1;
    for(int v = 0; v < g->nv; v++) {
        int vIdx = -1;
        int cost= minOverVtxDeleted2(g, i, v, &vIdx);
        if(cost > max) {
            max = cost;
            *vBest = v;
            *vBestIdx = vIdx;
        }
    }
    return max;
}

// Compute the fault cost of the graph. Uses method 2 with bitsets.
int computeCost(struct graph *g, struct options *options) {

    int min = INT_MAX;
    int minIdx = -1;
    int vBest = -1;
    int vBestIdx = -1;

    if(!g->mlsDegrees->used) {
        fprintf(stderr,
         "Error: Call method only after generateMlSubgraphs\n");
        exit(0);
    }

    for(size_t i = 0; i < g->mlsDegrees->used; i++) {
        int v = -1;
        int vIdx = -1;
        // int cost = maxOverVertices(g->nv, g->mlsDegrees->array[i],
         // g->vtxArrays, &v, &vIdx);
        int cost = maxOverVertices2(g, i, &v, &vIdx);
        if(cost < min) {
            min = cost;
            minIdx = i;
            vBest = v;
            vBestIdx = vIdx;
        }
    }

    if(options->verboseFlag) {
        fprintf(stderr, "Degrees of best ml-subgraph\n\t");
        for(int i = 0; i < g->nv; i++) {
            fprintf(stderr, "%d ", g->mlsDegrees->array[minIdx][i]);
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "Degrees of best ml-subgraph of G - %d\n\t",
         vBest);
        for(int i = 0; i < g->nv; i++) {
            fprintf(stderr, "%d ",
             g->vtxArrays[vBest]->array[vBestIdx][i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Fault cost is %d\n", min);
    }

    return min;
}

// Generate ml-subgraphs of g and store degree sequence data
void generateMlSubgraphs(struct graph *g, struct options *options,
 struct counters *counters, bitset deletedVertices) {

    // a is dynamic array containing distinct degree sequences of
    // ml-subgraphs. b contains bitsets giving leaves, linkes and
    // branches for each of these ml-subgraphs whose deg sequence was
    // stored in a.
    Array *a;
    BSArray *b;
    if(isEmpty(deletedVertices)) {
        a = g->mlsDegrees;
        b = g->mlsBitsets;
    }
    else {
        int v = next(deletedVertices,-1);
        a = g->vtxArrays[v];
        b = g->vtxBitsets[v];
    }

    if(isHam(g, deletedVertices)) {
        int *degrees = getHamCycleDegrees(g->nv, deletedVertices);
        insertArray(a, degrees, g->nv);
        bitset *degBitsets = malloc(sizeof(bitset) * 3);
        degBitsets[0] = EMPTY;
        degBitsets[1] = complement(deletedVertices, g->nv);
        degBitsets[2] = EMPTY;
        insertBitsetArray(b, degBitsets);
        return;
    }

    // Trying to computed forced endpoints for hamiltonian paths does
    // not yield speedup

    // If we are here: not ham
    bool isTrac = false;
    for(int i = 0; i < g->nv; i++) {
        for(int j = i+1; j < g->nv; j++) {
            if(contains(g->adjacencyList[i], j) && g->nv > 3) continue;
            if(hasHamPath(g, counters, deletedVertices, i, j)) {

                int *degrees = 
                 getHamPathDegrees(g->nv, deletedVertices, i, j);
                insertArray(a, degrees, g->nv);

                bitset *degBitsets = malloc(sizeof(bitset) * 3);
                degBitsets[0] = union(singleton(i), singleton(j));
                degBitsets[1] = complement(
                 union(degBitsets[0],deletedVertices), g->nv);
                degBitsets[2] = EMPTY;
                insertBitsetArray(b, degBitsets);
                isTrac = true;
            }
        }
    }

    // If not hamiltonian or traceable, generate all spanning trees.
    if(!isTrac)
        genSpanningTrees(g, options, counters, deletedVertices, NULL); 
}

// For each vertex of g and each ml-subgraph of G, check if v is a leaf.
// This gives upperbound on ml-number of G-v
void getUpperbounds(struct graph *g) {
    int ml = size(g->mlsBitsets->array[0][0]);
    bitset toCheck = complement(EMPTY, g->nv);
    for(size_t i = 0; i < g->mlsBitsets->used; i++) {
        forEach(v, toCheck) {
            if(contains(g->mlsBitsets->array[i][0], v)) {
                removeElement(toCheck, v);
                g->mlBound[v] = ml;
            }
        }
    }
}

// Check if graphs is 1-ham, => fc = 0 otherwise, initialise dynamic
// arrays for storing data on degree sequences of G and all G - v,
// generate ml-subgraphs which populates these dynamic arrays, then
// compute fault cost based on this data.
int getFaultCost(struct graph *g, struct options *options,
 struct counters *counters) {


    // ML-subgraph of g is ham cycle and of G-v are also ham cycles.
    if(isHam(g, EMPTY) && isK1Ham(g)) {
        if(options->verboseFlag) {
            fprintf(stderr, "1-hamiltonian -> Fault cost is 0.\n");
        }
        counters->oneHamiltonian++;
        return 0;
    }

    Array mlsDegrees;
    initArray(&mlsDegrees, 100*g->nv, g->nv);
    g->mlsDegrees = &mlsDegrees;
    Array *vtxArrays[g->nv];
    for(int i = 0; i < g->nv; i++) {
        vtxArrays[i] = malloc(sizeof(Array));
        initArray(vtxArrays[i], 100*g->nv, g->nv);
    }
    g->vtxArrays = vtxArrays;


    BSArray mlsBitsets;
    initBitsetArray(&mlsBitsets, 100*g->nv);
    g->mlsBitsets = &mlsBitsets;
    BSArray *vtxBitsets[g->nv];
    for(int i = 0; i < g->nv; i++) {
        vtxBitsets[i] = malloc(sizeof(BSArray));
        initBitsetArray(vtxBitsets[i], 100*g->nv);
    }
    g->vtxBitsets = vtxBitsets;

    generateMlSubgraphs(g, options, counters, EMPTY); 
    if(!isHam(g, EMPTY)) getUpperbounds(g);
    else {
        for(int i = 0; i < g->nv; i++) {
            g->mlBound[i] = 2;
        }
    }
    for(int i = 0; i < g->nv; i++) {
        generateMlSubgraphs(g, options, counters, singleton(i));
    }

    // if(options->verboseFlag) {
    //     for(int i = 0; i < g->nv; i++) {
    //         fprintf(stderr, "%d: %zu\n",i, vtxArrays[i]->used);
    //     }
    // }

    if(options->allFlag) {
        printMLSDegrees(g);
    }

    if(options->countFlag)
        countBranchesInMlsubgraph(g, counters);

    int faultCost = computeCost(g, options);

    for(int i = 0; i < g->nv; i++) {
        freeArray(g->vtxArrays[i]);
        free(g->vtxArrays[i]);
        freeBitsetArray(g->vtxBitsets[i]);
        free(g->vtxBitsets[i]);
    }
    freeArray(g->mlsDegrees);
    freeBitsetArray(g->mlsBitsets);

    // Remove later
    return faultCost;
}


//**********************************************************************
//
//              Leaf-Guaranteed
//
//**********************************************************************

// Check if ham or traceable, otherwise, generate spanning trees
// recursively to find one with min number of leaves.
int getMinimumLeafNumber(struct graph *g, struct options *options,
 struct counters *counters, bitset deletedVertices) {

    if(isHam(g, deletedVertices)) {
        return 1;
    }

    // Trying to computed forced endpoints for hamiltonian paths does
    // not yield speedup

    // If we are here: not ham
    for(int i = 0; i < g->nv; i++) {
        for(int j = i+1; j < g->nv; j++) {
            if(contains(g->adjacencyList[i], j)) continue;
            if(hasHamPath(g, counters, deletedVertices, i, j)) {
                return 2;
            }
        }
    }

    int ml;
    genSpanningTrees(g, options, counters, deletedVertices, &ml); 

    return ml;
}

// Given a k check if g is k-leaf-guaranteed
bool isKLeafGuaranteed(struct graph *g, struct options *options,
 struct counters *counters, int K) {

    int ml = getMinimumLeafNumber(g, options, counters, EMPTY);
    if(ml != K) return false;

    for(int i = 0; i < g->nv; i++) {
        int ml2 = getMinimumLeafNumber(g, options, counters,
         singleton(i));
        if(ml2 > ml) return false;
    }
    return true;
}

// Check if g is leaf-guaranteed for any k
int isLeafGuaranteed(struct graph *g, struct options *options,
 struct counters *counters) {

    int ml = getMinimumLeafNumber(g, options, counters, EMPTY);
    if(options->verboseFlag) {
        fprintf(stderr, "ml(G) = %d\n", ml);
    }

    for(int i = 0; i < g->nv; i++) {
        int ml2 = getMinimumLeafNumber(g, options, counters,
         singleton(i));
        if(ml2 > ml) return 0;
        if(options->verboseFlag) {
            fprintf(stderr, "ml(G - %d) = %d\n", i, ml2);
        }
    }
    return ml;
}

// int getMaxDegree(struct graph *g) {
//     int maxDeg = 0;
//     for(int i = 0; i < g->nv; i++) {
//         int deg = size(g->adjacencyList[i]);
//         if(deg > maxDeg) maxDeg = deg;
//     }
//     return maxDeg;
// }


// bool attainsMaxDegUpperBound(struct graph *g, struct options *options,
//  struct counters *counters) {

//     int maxDeg = getMaxDegree(g);

//     int ml = getMinimumLeafNumber(g, options, counters, EMPTY);

//     for(int i = 0; i < g->nv; i++) {
//         int ml2 = getMinimumLeafNumber(g, options, counters,
//          singleton(i));
//         if(ml2 == ml + maxDeg) return true;
//     }
//     return false;
// }

int main(int argc, char ** argv) {
    struct counters counters = {0};
    struct options options = {0};
    options.outputFrom = INT_MAX;
    char *valueToCompute = "with fault cost";
    int opt;
    int temp;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"all", no_argument, NULL, 'a'},
            {"help", no_argument, NULL, 'h'},
            {"count-branches", no_argument, NULL, 'c'},
            {"k-leaf-guaranteed", required_argument, NULL, 'l'},
            {"leaf-guaranteed", no_argument, NULL, 'L'},
            {"ml-number", no_argument, NULL, 'm'},
            {"output", required_argument, NULL, 'o'},
            {"output-from", required_argument, NULL, 'O'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "achl:Lmo:O:v", long_options,
         &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'a':
                options.allFlag = true;
                break;
            case 'c':
                options.countFlag = true;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'l':
                options.leafGuaranteedFlag = true; 
                options.leafGuaranteed = 
                 (int) strtol(optarg, (char **)NULL, 10);
                break;
            case 'L':
                options.leafGuaranteedFlag = true; 
                valueToCompute = "are leaf-guaranteed with k =";
                break;
            case 'm':
                options.mlFlag = true;
                valueToCompute = "with ml number";
                break;
            case 'o':
                options.output[options.numberOfOutputIntegers] = 
                 (int) strtol(optarg, (char **)NULL, 10);
                 options.numberOfOutputIntegers++;
                break;
            case 'O':
                temp = (int) strtol(optarg, (char **)NULL, 10);
                if(temp < options.outputFrom) options.outputFrom = temp;
                options.fromFlag = true;
                break;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./faultCost --help for more detailed"
                 " instructions.\n");
                return 1;
        }
    }

    unsigned long long int counter = 0;
    unsigned long long int passedGraphs = 0;
    unsigned long long int frequencies[MAXBITSETSIZE] =
     { [ 0 ... MAXBITSETSIZE-1 ] = 0 };

    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {
        struct graph g = {0};
        if(readGraph(graphString, &g, &options, &counters) == 1) {
            fprintf(stderr, "Error: problem loading graph.\n");
            exit(1);
        }

        counter++;

        // CHECK GRAPH PROPERTY HERE
        if(options.verboseFlag) {
            fprintf(stderr, "_______________________\n");
            fprintf(stderr, "%s", graphString);
            printGraph(g.adjacencyList, g.nv);
        }

        int value;

        if(options.leafGuaranteedFlag) {
            if(options.leafGuaranteed != 0) {
                if(isKLeafGuaranteed(&g, &options, &counters,
                 options.leafGuaranteed)) {
                // if(attainsMaxDegUpperBound(&g, &options, &counters)) {
                    printf("%s", graphString);
                    passedGraphs++;
                }
            }
            else {
                value = isLeafGuaranteed(&g, &options, &counters);
                if(value) {
                    printf("%s", graphString);
                    passedGraphs++;
                    frequencies[value]++;
                }
            }
            freeGraph(&g);
            continue;
        }
        if(options.mlFlag) {
            value = getMinimumLeafNumber(&g, &options, &counters,
             EMPTY);
        }
        else {
            value = getFaultCost(&g, &options, &counters);
        }
        frequencies[value]++;

        bool alreadyOutput = false;
        for(int i = 0; i < options.numberOfOutputIntegers; i++) {
            if(value == options.output[i]) {
                printf("%s", graphString);
                passedGraphs++;
                alreadyOutput = true;
                break;
            }
        }
        if(!alreadyOutput && 
         value >= options.outputFrom && options.fromFlag) {
            printf("%s", graphString);
            passedGraphs++;
        }

        for(int i = 0; i < g.nv + 1; i++) {
            resetSplayTree(&splayTreeArray[i]);
        }
        freeGraph(&g);
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    free(graphString);

    for (int i = 0; i < MAXBITSETSIZE; ++i) {
        if(frequencies[i] != 0) {
            fprintf(stderr, "\n \t%16lld graphs %s %d",
             frequencies[i], valueToCompute, i);
        }
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "Attempted to add: %llu. Already present: %llu\n",
    counters.timesAddedToTree, counters.degSequenceAlreadyPresent);

    if(options.countFlag) {
        for (int i = 0; i < MAXBITSETSIZE; ++i) {
            if(counters.nBranchesFreq[i] != 0) {
                fprintf(stderr, "%d branches: %4lld graphs, ",
                 i, counters.nBranchesFreq[i]);
            }
        }
        fprintf(stderr, "\n");
        for (int i = 0; i < MAXBITSETSIZE; ++i) {
            if(counters.nBranchesFreqSubgraphs[i] != 0) {
                fprintf(stderr, "%d branches: %4lld subgraphs, ",
                 i, counters.nBranchesFreqSubgraphs[i]);
            }
        }
        fprintf(stderr, "\n");
    }

    fprintf(stderr, "intermediateGraphs: %llu, times leaf num >"
    "ml: %llu, difference: %.2f%%\n", counters.intermediateGraphs,
    counters.mlPruned, 100.0 *
    counters.mlPruned/counters.intermediateGraphs);

    fprintf(stderr, 
     "call computeCost: %llu, prune: %llu, difference: %.2f%%\n",
     callCost, pruneInCost, 100.0 * pruneInCost/callCost);

    fprintf(stderr, 
     "Generate spanning: %llu, decrease ML: %llu, difference: %.2f%%\n",
     counters.spanning, counters.decreaseMl, 
     100.0 * counters.decreaseMl/counters.spanning);

    fprintf(stderr, "Times stored and later removed: %llu\n",
     counters.storedWithoutReason);

    fprintf(stderr, "Times executed hamiltonian path check: %llu\n",
     counters.doFullHamPathCheck);

    fprintf(stderr,
     "Checked %lld graphs in %f seconds: %llu passed.\n",
     counter, time_spent, passedGraphs);

    if(counters.skippedGraphs > 0) {
        fprintf(stderr, "Warning: %lld graphs were skipped.\n",
         counters.skippedGraphs);
    }

    if(counters.oneHamiltonian != frequencies[0]) {
        fprintf(stderr, 
         "Error: this should not happen: more fc0 than 1-ham.\n");
    }

    return 0;
}

