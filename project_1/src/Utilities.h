#pragma once
#include <iostream>
#include <array>
#include <map>
#include <tuple>
#include <vector>
#include <string>
#include <cstring>
#include <time.h>
#include <math.h>
#include <fstream>


/* Graph methods*/
#define MAX_NUM_NODES 120
#define MAX_NUM_EDGES 100000
#define MAX_NUM_VERTICES 120
#define MAX_NUM_ELE_IN_Q 100000


#define mCUT_S1 0
#define mCUT_S2 1
#define NUM_CUT_SETS 2
#define scaleMXflow 10e4 // to scale the fractional flow as na intger


template <typename T>
class Matrix2D {
    std::vector<T> data;
    unsigned int num_rows, num_cols;
public:

    Matrix2D()
        : num_rows(0), num_cols(0) {
    }

    Matrix2D(unsigned int x, unsigned int y)
        : num_rows(x), num_cols(y) {
        data.resize(num_rows * num_cols);
    }

    void resizeMatrix2D(unsigned int num_rows, unsigned int num_cols)
    {
        this->num_rows = num_rows;
        this->num_cols = num_cols;

        this->data.resize(num_rows * num_cols);
    }

    T& operator()(unsigned int x, unsigned int y) {
        if (x >= num_rows || y >= num_cols)
            throw std::out_of_range("OOB access"); // Throw something more appropriate



        return data[num_cols * x + y]; // Stride-aware access
    }


    void assignMatrix2D(unsigned int num_rows, unsigned int num_cols, T value)
    {
        this->num_rows = num_rows;
        this->num_cols = num_cols;

        this->data.assign(num_rows * num_cols, value);
    }

    unsigned int getNumRows()
    {
        return this->num_rows;
    }

    unsigned int getNumCols()
    {
        return this->num_cols;
    }

    unsigned int getNumElem()
    {
        return this->num_rows * this->num_cols;
    }



};
template <typename T>
class Matrix3D {
    std::vector<T> data;
    unsigned int num_matrices, num_rows, num_cols;
public:

    Matrix3D()
        : num_matrices(0), num_rows(0), num_cols(0) {
    }

    Matrix3D(unsigned int z, unsigned int x, unsigned int y)
        :num_matrices(z), num_rows(x), num_cols(y) {
        data.resize(num_matrices * num_rows * num_cols);
    }

    void resizeMatrix3D(unsigned int num_matrices, unsigned int num_rows, unsigned int num_cols)
    {
        this->num_matrices = num_matrices;
        this->num_rows = num_rows;
        this->num_cols = num_cols;

        this->data.resize(num_matrices * num_rows * num_cols);
    }

    T& operator()(unsigned int z, unsigned int x, unsigned int y) {
        if (x >= num_rows || y >= num_cols)
            throw std::out_of_range("OOB access"); // Throw something more appropriate



        return data[(num_cols * num_rows) * z + (num_cols * x) + y]; // Stride-aware access
    }




    void assignMatrix3D(unsigned int num_matrices, unsigned int num_rows, unsigned int num_cols, T value)
    {
        this->num_matrices = num_matrices;
        this->num_rows = num_rows;
        this->num_cols = num_cols;

        this->data.assign(num_matrices * num_rows * num_cols, value);
    }

    unsigned int getNumRows()
    {
        return this->num_rows;
    }

    unsigned int getNumCols()
    {
        return this->num_cols;
    }

    unsigned int getNumElem()
    {
        return this->num_matrices * this->num_rows * this->num_cols;
    }

};



typedef struct queue {
    std::vector<unsigned short> elem;
    // unsigned short elem[MAX_NUM_ELE_IN_Q];
    unsigned short num_el;
    unsigned short next;
}queue_t;


typedef struct edge {
    unsigned short first;
    unsigned short second;
}edge_t;

typedef struct arch {
    unsigned short a1; // source node
    unsigned short a2; // dest. node
    unsigned short q;// capacity
    unsigned short flow; // flow

    unsigned short rev; // index of reverse arc

}arch_t;

typedef struct un_graph
{
    unsigned short num_nodes;
    unsigned short num_edges;
    unsigned short num_cc;
    //edge_t edge[MAX_NUM_EDGES];
    std::vector<edge_t> edge;
    // unsigned short nodes_map[MAX_NUM_NODES]; // maps vertex in the solution 
    std::vector<unsigned short> nodes_map;
    //  short conn_comp[MAX_NUM_NODES]; 	// Connected component
    std::vector<short> conn_comp;

}un_graph_t;


typedef struct di_graph {
    unsigned short num_ver;
    // unsigned short ver_label[MAX_NUM_VERTICES]; // var_label[graph_index] = original_model
    std::vector<short> ver_label;
    //  short level[MAX_NUM_VERTICES]; // level of a node
    std::vector<short> level;
    //  arch_t adj_archs[MAX_NUM_VERTICES][MAX_NUM_ARCHS_PER_V]; //  adj[i] stores the archs i->j
    Matrix2D<arch_t> adj_archs;
    // unsigned short adj_size[MAX_NUM_VERTICES]; // size of adj_archs
    std::vector<unsigned short> adj_size;
}di_graph_t;

/*////////////////*/
