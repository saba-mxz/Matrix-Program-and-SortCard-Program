// matlib24.cpp : This file contains the 'main' function and others to test matrix routines
//

// matlib24.cpp : A simple matrix library program.
// Template Written by Prof Richard Mitchell    6/12/23
// Adapted by:   << put your name here >>

#include <string>
#include <iostream>
using namespace std;

const int maxM = 16;	// allows for matrices up to size 4*4

struct myMat {			// structure to store a matrix
    int numRows;		// number of rows
    int numCols;		// number of columns
    int data[maxM];		// data are stored in row order
};

myMat zeroMat(int r, int c) {
    // create a matrix with r rows and c columns, filled with zeros
    myMat m;			// define matrix
    m.numRows = r;		// set size
    m.numCols = c;
    for (int ct = 0; ct < maxM; ct++) m.data[ct] = 0;	// set elems to 0
    return m;			// return matrix
}

int getIndex(myMat m, int r, int c) {
    // returm index of element in row r, col c of matrix m
    return r * m.numCols + c;
}

myMat matError(string errstr) {
    // if this is called, an error has occured ... output errstr and return 0 size myMat
    cout << "Error : " << errstr << "\n";
    return zeroMat(0, 0);
}

int intError(string errstr) {
    // if this is called, an error has occured ... output errstr and return 0
    cout << "Error : " << errstr << "\n";
    return 0;
}

myMat mFromStr(string s) {
    // create a matrix from string s
    // string of form "1,2;3,4;5,6"   ; separates rows and , separates columns ... No error check
    int ms;
    if (s.length() > 0) ms = 1; else ms = 0;
    myMat m = zeroMat(ms, ms);						// is s empty create 0*0 matrix, else start with 1*1 matrix
    int mndx = 0;									// used to index into array
    string sub = "";								// sub string - numbers between , or ; set empty
    for (int ct = 0; ct < s.length(); ct++) {		// each char in turn
        if ((s[ct] == ',') || (s[ct] == ';')) {	// if , or ; then sub contains number
            m.data[mndx++] = stoi(sub);				// convert number to integer, put in data array
            sub = "";								// clear sub string
            if (s[ct] == ';') m.numRows++;			// if found ; indicates an extra row
            else if (m.numRows == 1) m.numCols++;	// if , then (if in row 1) increase count of columns
        }
        else sub = sub + s[ct];						// add character to sub string
    }
    if (sub.length() > 0) m.data[mndx++] = stoi(sub);// add last sub string
    return m;
}

void printMat(const char* mess, myMat m) {
    // mess is string to be printed, followed by matrix m
    cout << mess << " = " << "\n";				// print message
    for (int r = 0; r < m.numRows; r++) {		// do each row
        for (int c = 0; c < m.numCols; c++)		// do each column
            cout << m.data[getIndex(m, r, c)] << "\t";	// outputing the element then tab
        cout << "\n";							// output new line at end of row
    }
    cout << "\n";								// and end of Matrix
}


myMat mGetRow(myMat m, int row) {
    // create a matrix from m, having one row
    myMat res = zeroMat(1, m.numCols);		// create a matrix of 1 row
    for (int col = 0; col < m.numCols; col++)		// for each element in row
        res.data[col] = m.data[getIndex(m, row, col)];		// copy col element to res
    return res;
}


myMat mGetCol(myMat m, int col) {
    // create a matrix from m, having one column
    myMat res = zeroMat(m.numRows, 1); // create a matrix of m.numRows rows and 1 column

    for (int row = 0; row < m.numRows; row++) {
        res.data[row] = m.data[getIndex(m, row, col)]; // copy elements from the specified column
    }

    return res;
}


myMat mSetCol(myMat m, int col, myMat v) {
    // insert v into given col in m
    if (m.numRows != v.numRows)
        return matError("Matrix/Vector should have same number of rows");
    else {
        myMat res = m;
        for (int row = 0; row < m.numRows; row++)
            res.data[getIndex(m, row, col)] = v.data[row];
        return res;
    }
}

int dotProd(myMat v1, myMat v2) {
    // Calculate the dot product of two vectors v1 and v2
    // Vectors can be either row or column vectors

    // Check if the vectors have the same dimensions
    if ((v1.numRows != 1 && v1.numCols != 1) || (v2.numRows != 1 && v2.numCols != 1) || (v1.numRows * v1.numCols != v2.numRows * v2.numCols)) {
        return intError("Vectors must have the same dimensions for dot product calculation");
    }

    int result = 0;

    // If both vectors are row vectors or column vectors
    if ((v1.numRows == 1 && v2.numRows == 1) || (v1.numCols == 1 && v2.numCols == 1)) {
        for (int i = 0; i < v1.numRows * v1.numCols; i++) {
            result += v1.data[i] * v2.data[i];
        }
    }
        // If one vector is a row vector and the other is a column vector
    else {
        // Make sure the dimensions match for row and column vectors
        if ((v1.numRows != v2.numCols) || (v1.numCols != v2.numRows)) {
            return intError("Dimensions do not match for dot product calculation");
        }

        // Perform dot product between row and column vectors
        for (int i = 0; i < v1.numCols; i++) {
            result += v1.data[i] * v2.data[i];
        }
    }

    return result;
}

void testVecs(myMat m1, myMat m3) {
    // test vector routines ... get row from m1, col from m3, do dot product of these
    cout << "Testing Vector routines" << "\n";
    printMat("m1 row 0", mGetRow(m1, 0));	// display row 0 of m1
    printMat("m3 col 1", mGetCol(m3, 1));	// display col 1 of m2
    cout << "Dot prod of these is " << dotProd(mGetRow(m1, 0), mGetCol(m3, 1)) << "\n\n";
    cout << "Dot prod of m1 row 1 and m3 row 1 " << dotProd(mGetRow(m1, 0), mGetRow(m3, 1)) << "\n\n";
}


myMat mTranspose(myMat m) {
    // Generate a new matrix which is the transpose of m
    myMat res;
    res.numRows = m.numCols;  // Number of rows in the transpose matrix is equal to the number of columns in the original matrix
    res.numCols = m.numRows;  // Number of columns in the transpose matrix is equal to the number of rows in the original matrix

    // Iterate through each element of the original matrix
    for (int r = 0; r < m.numRows; r++) {
        for (int c = 0; c < m.numCols; c++) {
            // Copy the element to the transposed matrix, but swap row and column indices
            res.data[getIndex(res, c, r)] = m.data[getIndex(m, r, c)];
        }
    }

    return res;
}


myMat mAdd(myMat m1, myMat m2) {
    // Create a new matrix whose elements are the sum of the equivalent elements in m1 and m2
    // Make sure the matrices have the same dimensions
    if (m1.numRows != m2.numRows || m1.numCols != m2.numCols) {
        return matError("Matrices must have the same dimensions for addition");
    }

    // Create a new matrix to store the result
    myMat res;
    res.numRows = m1.numRows;
    res.numCols = m1.numCols;

    // Iterate through each element of the matrices and compute the sum
    for (int i = 0; i < m1.numRows * m1.numCols; i++) {
        res.data[i] = m1.data[i] + m2.data[i];
    }

    return res;
}


myMat mScalarMultDiv(myMat m, int s, int isMult) {
    // multiply or divide all elements in m by s
    myMat res = zeroMat(0, 0);		// change arguments
    // write code to do multiply or divide by scalar
    return res;
}

myMat mMult(myMat m1, myMat m2) {
    // Check if the matrices can be multiplied
    if (m1.numCols != m2.numRows) {
        return matError("Matrices cannot be multiplied: Number of columns of first matrix must be equal to number of rows of second matrix");
    }

    // Create a new matrix to store the result
    myMat res;
    res.numRows = m1.numRows;
    res.numCols = m2.numCols;

    // Iterate through each element of the resulting matrix
    for (int r = 0; r < res.numRows; r++) {
        for (int c = 0; c < res.numCols; c++) {
            // Initialize the element in the resulting matrix to 0
            res.data[getIndex(res, r, c)] = 0;
            // Calculate the value of the element by summing the products of corresponding elements from rows of m1 and columns of m2
            for (int k = 0; k < m1.numCols; k++) {
                res.data[getIndex(res, r, c)] += m1.data[getIndex(m1, r, k)] * m2.data[getIndex(m2, k, c)];
            }
        }
    }

    return res;
}


void testMatOps(myMat m1, myMat m2, myMat m3) {
    // test matrix operations m1 + m2; m1 + m3 (not poss) m1 + m3'; m1 * m3; m3 * m2
    cout << "Testing Add, Transpose, Multiply routines" << "\n";
    printMat("m1+m2", mAdd(m1, m2));
    printMat("m1 + m3", mAdd(m1, m3));
    printMat("m1 + m3'", mAdd(m1, mTranspose(m3)));
    printMat("m1*m3", mMult(m1, m3));
    printMat("m3*m1", mMult(m3, m1));
    printMat("m1*m2", mMult(m1, m2));

    myMat m1_transpose = mTranspose(m1);
    printMat("Transpose of m1", m1_transpose);

    myMat m2_transpose = mTranspose(m2);
    printMat("Transpose of m2", m2_transpose);

    myMat m3_transpose = mTranspose(m3);
    printMat("Transpose of m3", m3_transpose);

}



myMat mSubMat(myMat A, int row, int col) {
    // Check if the specified row and column are within the bounds of the matrix
    if (row < 0 || row >= A.numRows || col < 0 || col >= A.numCols) {
        return matError("Invalid row or column index");
    }

    // Create a new matrix to store the result
    myMat res;
    res.numRows = A.numRows - 1;  // Decrease the number of rows by 1
    res.numCols = A.numCols - 1;  // Decrease the number of columns by 1

    // Iterate through each element of the resulting matrix
    int resIndex = 0;  // Index for the resulting matrix
    for (int r = 0; r < A.numRows; r++) {
        if (r == row) {
            // Skip the specified row
            continue;
        }
        for (int c = 0; c < A.numCols; c++) {
            if (c == col) {
                // Skip the specified column
                continue;
            }
            // Copy the element to the resulting matrix
            res.data[resIndex] = A.data[getIndex(A, r, c)];
            resIndex++;
        }
    }

    return res;
}


int mDet(myMat A) {
    // Check if the matrix is square
    if (A.numRows != A.numCols) {
        return intError("Matrix must be square for determinant calculation");
    }

    // Base case: If the matrix is 1x1, return the single element as the determinant
    if (A.numRows == 1) {
        return A.data[0];
    }

    // Base case: If the matrix is 2x2, return the determinant using the formula ad - bc
    if (A.numRows == 2) {
        return A.data[0] * A.data[3] - A.data[1] * A.data[2];
    }

    int det = 0;

    // Recursive case: If the matrix is larger than 2x2, compute the determinant using cofactor expansion along the first row
    for (int c = 0; c < A.numCols; c++) {
        // Calculate the cofactor for each element in the first row
        int cofactor = A.data[c] * mDet(mSubMat(A, 0, c));
        // Add or subtract the cofactor based on its position
        det += (c % 2 == 0 ? 1 : -1) * cofactor;
    }

    return det;
}



myMat mAdj(myMat A) {
    // Check if the matrix is square
    if (A.numRows != A.numCols) {
        return matError("Matrix must be square for adjoint calculation");
    }

    // Create a new matrix to store the adjoint
    myMat adj;
    adj.numRows = A.numRows;
    adj.numCols = A.numCols;

    // Iterate through each element of the matrix
    for (int r = 0; r < A.numRows; r++) {
        for (int c = 0; c < A.numCols; c++) {
            // Calculate the cofactor of the element
            int cofactor = ((r + c) % 2 == 0 ? 1 : -1) * mDet(mSubMat(A, r, c));
            // Assign the cofactor to the adjoint matrix, transposed
            adj.data[getIndex(adj, c, r)] = cofactor;
        }
    }

    return adj;
}

void testMatEqn(myMat A2, myMat b) {
    // Check if the matrix A2 is square and has dimensions 2x2
    if (A2.numRows != 2 || A2.numCols != 2) {
        cout << "Matrix A2 must be a 2x2 square matrix for solving equations using Cramer's rule" << endl;
        return;
    }

    // Calculate the determinant of A2
    int detA2 = mDet(A2);

    // Check if the determinant is zero
    if (detA2 == 0) {
        cout << "Matrix A2 is singular, cannot solve equations" << endl;
        return;
    }

    // Create a matrix to store the solution
    myMat x;
    x.numRows = 2;
    x.numCols = 1;

    // Calculate the determinants of matrices obtained by replacing columns of A2 with b
    int detA2_1 = mDet(mSetCol(A2, 0, b));
    int detA2_2 = mDet(mSetCol(A2, 1, b));

    // Calculate the components of the solution vector x using Cramer's rule
    x.data[0] = detA2_1 / detA2;
    x.data[1] = detA2_2 / detA2;

    // Print the solution matrix x
    printMat("(Cramer's Rule Answer for 2*2)\nx", x);
}


myMat solveMatrixEquation(myMat A, myMat b) {
    // Check if the matrix A is square and has dimensions 3x3
    if (A.numRows != 3 || A.numCols != 3) {
        cout << "Matrix A must be a 3x3 square matrix for solving equations using Cramer's rule" << endl;
        return matError("");
    }

    // Calculate the determinant of A
    int detA = mDet(A);

    // Check if the determinant is zero
    if (detA == 0) {
        cout << "Matrix A is singular, cannot solve equations" << endl;
        return matError("");
    }

    // Create a matrix to store the solution
    myMat x;
    x.numRows = 3;
    x.numCols = 1;

    // Calculate the determinants of matrices obtained by replacing columns of A with b
    int detA1 = mDet(mSetCol(A, 0, b));
    int detA2 = mDet(mSetCol(A, 1, b));
    int detA3 = mDet(mSetCol(A, 2, b));

    // Calculate the components of the solution vector x using Cramer's rule
    x.data[0] = detA1 / detA;
    x.data[1] = detA2 / detA;
    x.data[2] = detA3 / detA;

    return x;
}



int main()
{
    cout << "32009512\n";	// change to your student number
    myMat m1, m2, m3, A, A2, b, A3, b3;						// create  matrices

    m1 = mFromStr("9, 1, 10; 7, 8, 10");			// change numbers to those in A from Q1 on web page, as specified on the sheet
    m2 = mFromStr("10, 9, 1; 10, 8, 8");			// ditto but B
    m3 = mFromStr("8, 7; 5, 8; 8, 10");			// ditto  but C
    A = mFromStr("8,10; 9,1");
    A2 = mFromStr("10, 7; 10, 9");
    b = mFromStr("163;181");
    A3 = mFromStr("9, 1, 5; 10, 8, 10; 9, 8, 7");
    b3 = mFromStr("137; 236; 202");

    printMat("A", A);
    printMat("A2", A2);
    printMat("b", b);
    printMat("A3", A3);
    printMat("b", b3);
    printMat("m1", m1);						// display m1
    printMat("m2", m2);						// display m2
    printMat("m3", m3);						// display m3

    cout << "Determinant of matrix A: " << mDet(A) << endl;

    // Test mSubMat function
    int row = 0;
    int col = 1;
    myMat submatrix = mSubMat(A, row, col);
    cout << "Submatrix of A with row " << row << " and col " << col << " removed:" << endl;
    printMat("", submatrix);

    // Calculate the adjoint of matrix A
    myMat adj_A = mAdj(A);
    cout << "Adjoint of matrix A:" << endl;
    printMat("", adj_A);

    // Solve the equation using Cramer's rule
    testMatEqn(A2, b);

    // Call the solveMatrixEquation function with A3 and b3
    myMat solution = solveMatrixEquation(A3, b3);
    printMat("(Cramer's Rule Answer for 3*3)", solution);

    return 0;
}//
// Created by Saba Mirza on 22/10/2024.
//
