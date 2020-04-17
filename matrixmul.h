#ifndef MATRIX_MUL
#define MATRIX_MUL

// Constant(s)

#define MAX_ELEMENTS 10

// Type(s)

typedef struct MatrixStruct {
    unsigned char rows;
    unsigned char columns;
    double elements[MAX_ELEMENTS];
} Matrix;

// Function declarations

void setRowsColumns(unsigned char rows, unsigned char columns, Matrix *);
void setElement(unsigned char row, unsigned char column, double element, Matrix *);
double getElement(unsigned char row, unsigned char column, Matrix *);
void clearMatrix(Matrix *);
int multiply(Matrix *matrixA, Matrix *matrixB, Matrix *matrixC);

// Function implementations

void setRowsColumns(unsigned char rows, unsigned char columns, Matrix *matrix) {

    matrix->rows = rows;
    matrix->columns = columns;

    return;
}

void setElement(unsigned char row, unsigned char column, double element, Matrix *matrix) {

    matrix->elements[(row - 1) * matrix->columns + (column - 1)] = element;

    return;
}

double getElement(unsigned char row, unsigned char column, Matrix *matrix) {
    return matrix->elements[(row - 1) * matrix->columns + (column - 1)];
}

void clearMatrix(Matrix *matrix) {

    for (int i = 0; i < MAX_ELEMENTS; i++) {
        matrix->elements[i] = 0;
    }

    return;
}

int multiply(Matrix *matrixA, Matrix *matrixB, Matrix *matrixC) {

    // un-equal number of columns and rows yield an error value
    if (matrixA->columns != matrixB->rows) {
        return -1;
    }

    // set the result matrix's rows and columns numbers and initialize its elements to zero
    setRowsColumns(matrixA->rows, matrixB->columns, matrixC);
    clearMatrix(matrixC);

    for (int i = 1; i <= matrixC->rows; i++) {
        for (int j = 1; j <= matrixC->columns; j++) {
            for (int k = 1; k <= matrixA->columns; k++) {
                setElement(i, j, getElement(i, j, matrixC) + getElement(i, k, matrixA) * getElement(k, j, matrixB), matrixC);
            }
        }
    }

    return 0;
}

//transposed
int transposed_v2(Matrix* matrixA, Matrix* matrixB) {

    for (int i = 1; i <= matrixA->rows; i++) {
        for (int j = 1; j <= matrixA->columns; j++) {
                setElement(j, i, getElement(i, j, matrixA) , matrixB);
        }
    }

    return 0;
}
#endif