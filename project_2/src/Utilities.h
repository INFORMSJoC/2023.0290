#pragma once
#include "Constant.h"
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
        //if (x >= num_rows || y >= num_cols)
        //    throw std::out_of_range("OOB access"); // Throw something more appropriate
       // return data[num_rows * y + x]; // Stride-aware access

        return data[num_cols * x + y]; // Stride-aware access
    }



    void assignMatrix2D(unsigned int num_rows, unsigned int num_cols,  T value)
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

    Matrix3D(unsigned int z,unsigned int x, unsigned int y)
        :num_matrices(z), num_rows(x), num_cols(y) {
        data.resize(num_matrices * num_rows * num_cols);
    }

    void resizeMatrix3D(unsigned int num_matrices,unsigned int num_rows, unsigned int num_cols)
    {
        this->num_matrices = num_matrices;
        this->num_rows = num_rows;
        this->num_cols = num_cols;

        this->data.resize(num_matrices * num_rows * num_cols);
    }

    T& operator()(unsigned int z, unsigned int x, unsigned int y) {
        //if (x >= num_rows || y >= num_cols)
        //    throw std::out_of_range("OOB access"); // Throw something more appropriate
       // return data[(num_matrices * x) + z + (num_rows * y) + x]; // Stride-aware access

     // return data[(num_matrices * num_rows)*y + (num_rows * x) +z]; // Stride-aware access

        //return data[(num_cols * num_rows) * z + (num_rows * x) + y]; // Stride-aware access

        return data[(num_cols * num_rows) * z + (num_cols * x) + y];

    }




    void assignMatrix3D(unsigned int num_matrices, unsigned int num_rows, unsigned int num_cols,  T value)
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



