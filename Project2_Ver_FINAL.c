//**********Project 2**********//
// Name: Yuxuan Zhang
// UHiD: 2137717
// Date: 05,06,2023
// Last Update: 05,12,2023

//**********WARNING**********//
// This programme using the method of "malloc" to access memory.
// The return value of the function in this program is the address. 
// After using it, please use "free (address)" to avoid memory leaks (in the main function).

//**********HEAD FILES**********//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//**********MACRO DEFINE**********//
#define _LUP_   // Using the method of LUP decomposition. If you want to test the LU method please change the "_LUP_" to "_LU_"

//**********DECLARATION**********//
int* get_Size(char* file);
double** get_Matrix(char* file, int* size_matrix);
void write_Matrix(double** matrix, int* size_matrix, char* file_name);
void convert_P_vector2matrix(double* P_vector,double** P_matrix, int* size_matrix);
double*** LUP_decomposition(double** matrix, int* size_matrix);
double*** LU_decomposition(double** matrix, int* size_matrix);

//**********MAIN FUNCTION**********//
// Main funtion of this project
int main(int argc, char** argv){
    // Declaration of the pointer and variables we need.
    double*** result = NULL;                                            // The result of LUP result and LU result. 
                                                                        // result[0] is the address of Lower matrix. 
                                                                        // result[1] is the address of Upper matrix. 
                                                                        // result[2] is the address of P matrix.
    double** matrix = NULL;                                             // Store the address of input matrix.
    double** L_matrix = NULL;                                           // Store the address of Lower matrix.
    double** U_matrix = NULL;                                           // Store the address of Upper matrix.
    double** P_matrix = NULL;                                           // Store the address of P matrix.
    int* size_matrix = NULL;                                            // Store the address of matrix's size.
    char* text_name = NULL;                                             // The text name which will be wrote or read.  
    int i;                                                              // Variable i for counting loops' numbers.
    // Get the text message from the commender.
    if(argc <= 1){
        return 0;                                                       // If the enter number is too small, return the main function.
    }else{
        text_name = argv[1];                                            // Store the name of text from commender (argv[1]).
    }
    // Get the size of the matrix.
    size_matrix = get_Size(text_name);
    // Get the matrix from the "textname.txt" with known size.
    matrix = get_Matrix(text_name,size_matrix);
    
    // If we want to Using LU decomposition method.  define _LU_
    #ifdef _LU_
    // Using method of LU_decomposition and Return the address of LUP matrix. (P is an Identity matrix)
    result = LU_decomposition(matrix, size_matrix);
    L_matrix = result[0];                                               // Store L matrix.
    U_matrix = result[1];                                               // Store U matrix.
    P_matrix = result[2];                                               // Store P matrix.
    // Write the matrix into text.
    text_name = "L.txt";write_Matrix(L_matrix,size_matrix,text_name);   // Re-write L matrix into a text which is named "L.txt".
    text_name = "U.txt";write_Matrix(U_matrix,size_matrix,text_name);   // Re-write U matrix into a text which is named "U.txt".
    text_name = "P.txt";write_Matrix(P_matrix,size_matrix,text_name);   // Re-write P matrix into a text which is named "P.txt".
    #endif

    // If we want to Using LUP decomposition method. define _LUP_
    #ifdef _LUP_
    // Using method of LUP_decomposition and Return the address of LUP matrix
    result = LUP_decomposition(matrix, size_matrix);
    L_matrix = result[0];                                               // Store L matrix.
    U_matrix = result[1];                                               // Store U matrix.
    P_matrix = result[2];                                               // Store P matrix.
    // Write the matrix into text.
    text_name = "L.txt";write_Matrix(L_matrix,size_matrix,text_name);   // Re-write L matrix into a text which is named "L.txt".
    text_name = "U.txt";write_Matrix(U_matrix,size_matrix,text_name);   // Re-write U matrix into a text which is named "U.txt".
    text_name = "P.txt";write_Matrix(P_matrix,size_matrix,text_name);   // Re-write P matrix into a text which is named "P.txt".
    #endif

    // Free the memory of the matrixes. In case of the "memory leak" problem happen
    for (i = 0; i < *size_matrix; i++){
        free(matrix[i]);
    }
    free(matrix);
    for (i = 0; i < *size_matrix; i++){
        free(L_matrix[i]);
    }
    free(L_matrix);
    for (i = 0; i < *size_matrix; i++){
        free(U_matrix[i]);
    }
    free(U_matrix);
    for (i = 0; i < *size_matrix; i++){
        free(P_matrix[i]);
    }
    free(P_matrix);
    free(size_matrix);
    // End of the main function.
    return 0;
}

//**********GET SIZE OF MATRIX FROM TEXT**********//
// This function can realize that return a matrix size by reading of a text file.
int* get_Size(char* file){                          // char* file is the address of file name.
    // Declaration of variables
    int* size_matrix = (int*)malloc(1*sizeof(int)); // Store the address of the value of size of matrix.
    char* buff=(char*)malloc(100*sizeof(char));     // This is a buffer for the strings for the entering number. (MAX = 100)
    //
    FILE *fp = NULL;                                // Define the address of the file pointer fp (file).
    fp = fopen(file, "r");                          // Begain read the value from file.
    // 
    fscanf(fp, "%s", buff);                         // Put the value which we get from the file into buffer.
    *size_matrix = atoi(buff);                      // Convert the "string" into "int".
    // Free the memory of the variables.
    fclose(fp);
    free(buff);
    // End of file with return vlaue int size_matrix.
    return size_matrix;
}

//**********GET MATRIX FROM TEXT*********//
// This function can realize that return a matrix address by reading of a text file.
double** get_Matrix(char* file, int* size_matrix){                  // char* file is the address of file name. int* size_matrix is the address of number of size.
    // Declaration of variables.
    double** matrix = NULL;                                         // Store matrix pointer type of (double**).
    char* buff=(char*)malloc(50*sizeof(char));                      // This is a buffer for the strings for the entering floating number. (MAX = 50)
    int cols = 0;int rows = 0;int i;                                // Variables for counting loops' numbers. (cols -> column, rows -> row)
    // Start file operations
    FILE *fp = NULL;                                                // Define the address of the file pointer fp (file).
    fp = fopen(file, "r");                                          // Begain read the value from file.
    // Do not read the first value of the file.
    fscanf(fp, "%s", buff);
    // Initialize the matrix's memory by malloc.
    matrix=(double**)malloc(sizeof(double*)*(*size_matrix));        // The size of matrix's rows.
    for(i = 0;i < *size_matrix;i++){
        matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));   // The size of matrix's cols.
    }     
    // Initialize the matrix's value by two "for" loops.
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0; cols < *size_matrix; cols++){
            fscanf(fp, "%s", buff);                                 // Get the "string" from the text.
            matrix[rows][cols] = atof(buff);                        // Convert the "string" into "int".
        }
    }
    // Free the memory of the variables.
    fclose(fp);
    free(buff);
    // End of file with return vlaue double** matrix.
    return matrix;
}

//**********WRITE MATRIX INTO TEXT**********//
// This function can write a matrix into a file(named *file_name) by its address. 
void write_Matrix(double** matrix, int* size_matrix, char* file_name){  // double** matrix is the address of the matrix. int* size_matrix is the address of number of size. char* file_name is the address of file name.
    // Declaration of variables.
    int cols, rows;                                     // Variables for counting loops' numbers. (cols -> column, rows -> row)
    char* buff=(char*)malloc(100*sizeof(char));         // This is a buffer for the strings for the floating number. (MAX = 50)
    // Start file operations
    FILE *fp = NULL;                                    // Define the address of the file pointer fp (file).
    fp = fopen(file_name, "w+");                        // Begain write the value into file.
    //  First of all Write the size of matrix into the text file
    sprintf(buff,"%d",*size_matrix);                    // Get the value of size of matrix into buffer.
    fputs(buff, fp);                                    // Put the (int)size of matrix into the first line of text file.
    fputs("\n", fp);                                    // Every end of line we put "\n" to change line. 
    //
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0; cols < *size_matrix; cols++){   
            sprintf(buff,"%0.18E",matrix[rows][cols]);  // Get the value of matrix into buffer.
            fputs(buff, fp);                            // Put the (double)number of matrix into text file.
            fputs("\t", fp);                            // Every end of value we put "\t". 
        }
        fputs("\n", fp);                                // Every end of line we put "\n" to change line. 
    }
    // Free the memory of the variables.
    fclose(fp);
    free(buff);
    return;
}

//**********CONVERT VECTOR TO MATRIX**********//
// This function converts a one-dimensional vector into an N-dimensional vector. By supplementing zero and one.
void convert_P_vector2matrix(double* P_vector,double** P_matrix, int* size_matrix){ // double* P_vector is address of P vector, double** P_matrix is the address of P matrix, int* size_matrix, char* file_name is the address of file name.
    // Declaration of variables.
    int rows,cols;                          // Variables for counting loops' numbers. (cols -> column, rows -> row)
    // Go through all the cells of matrix.
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0; cols < *size_matrix; cols++){
            if(cols == (int)P_vector[rows]){
                P_matrix[cols][rows] = 1.0; // When the cols are equal to the value of P_vector's current value. (Set 1.0)
            }else{
                P_matrix[cols][rows] = 0.0; // When the cols are equal to the value of P_vector's current value. (Set 0.0)
            }
        }
    }
}

//**********LU DECOMPOSITION METHOD**********//
// This function using LU decomposition method (P is an Identity matrix)
double*** LU_decomposition(double** matrix, int* size_matrix){                      // double** matrix is the address of matrix, int* size_matrix is the address of value size of matrix
    // Declaration of variables.
    double*** LU_result = (double***)malloc(sizeof(double**)*3);                    // Requests for 3*(double**) sized addresses are stored in the LU_result.
    double** L_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the L_matrix.
    double** U_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the U_matrix.
    double** P_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the P_matrix.
    double* P_vector = (double*)malloc(sizeof(double)*(*size_matrix));              // Requests for (size_matrix)*(double) sized addresses are stored in the P_vector.
    int cols, rows, i, j;                                                           // Variables for counting loops' numbers. (cols -> column, rows -> row)
    double sum_Matrix;                                                              // Stroe the value of sum of the matrix multiply value.
    // Requests for memory for L_matrix, U_atrix and P_matrix
    for( i = 0; i < *size_matrix; i++){
        L_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the L_matrix[i].
        U_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the U_matrix[i].
        P_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the P_matrix[i].
    }
    // Initialize the Lower matrix (Identity matrix)
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0; cols < *size_matrix; cols++){
            if(rows == cols){
                L_matrix[rows][cols] = 1.0;                                         // The value of the diagonal of the matrix is 1.
            }else{  
                L_matrix[rows][cols] = 0.0;                                         // The value of the matrix in other position is 0.
            }
        }
    }
    // Initialize Upper matrix
    for(rows = 0;rows < *size_matrix; rows++){
        for(cols = 0;cols < *size_matrix; cols++){
            U_matrix[rows][cols] = matrix[rows][cols];                              // Upper matrix is the same matrix of matrix we enter in.
        }
    }
    // Initialize P vector
    for(i = 0; i < *size_matrix; i++){
        P_vector[i] = i;                                                            // P matrix is identity matrix in LU decomposition method
    }
    convert_P_vector2matrix(P_vector,P_matrix, size_matrix);                        // Using function convert_P_vector2matrix(double* ,double** , int* ) to convert the vector.
    // Beginning of LU demcomposion
    for(i = 0;i < *size_matrix;i++){
        // L_matrix calculation
        for(rows = i+1; rows< *size_matrix; rows++){
            L_matrix[rows][i] = U_matrix[rows][i] / U_matrix[i][i];                 // The value below the diagonal of the L matrix is the value of each row of the U matrix except the principal element divided by the value of the principal element.
        }
        // U_matrix calculation
        for(rows = i; rows< *size_matrix; rows++){
            for(cols = i; cols < *size_matrix; cols++){
                // Gaussian elimination
                sum_Matrix = (double)0;                                             // Initialize the value of sum_matrix. Matrix multiplication
                for(j = i; j < *size_matrix; j++){
                    if(j == rows){
                        sum_Matrix += 1.0*U_matrix[j][cols];                        // Do not change the value when j = rows.
                    }else{
                        sum_Matrix += L_matrix[rows][j]*U_matrix[j][cols]*(-1.0);   // Do elimination when j != rows.
                    }
                }
                U_matrix[rows][cols] = sum_Matrix;                                  // Calculate the value of Matrix multiplication
            }
        }
    }
    // Return the result of Lower matrix, Upper matrix, and P matrix.
    LU_result[0] = L_matrix;
    LU_result[1] = U_matrix;
    LU_result[2] = P_matrix;
    //  Free the memory of the matrixes. In case of the "memory leak" problem happ
    free(P_vector);
    // End of file with return vlaue double*** LU_result
    return LU_result;
}

//**********LUP DECOMPOSITION METHOD**********//
double*** LUP_decomposition(double** matrix, int* size_matrix){                     // double** matrix is the address of matrix, int* size_matrix is the address of value size of matrix
    // Declaration of variables.
    double*** LUP_result = (double***)malloc(sizeof(double**)*3);                   // Requests for 3*(double**) sized addresses are stored in the LUP_result.
    double** L_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the L_matrix.
    double** U_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the U_matrix.
    double** P_matrix = (double**)malloc(sizeof(double*)*(*size_matrix));           // Requests for (size_matrix)*(double*) sized addresses are stored in the P_matrix.
    double* P_vector = (double*)malloc(sizeof(double)*(*size_matrix));              // Requests for (size_matrix)*(double) sized addresses are stored in the P_vector.
    int cols, rows, i, j, max_pivot_index;                                          // Variables for counting loops' numbers. (cols -> column, rows -> row) 
                                                                                    // max_pivot_index is the position where max pivot exist.
    double pivot, buff, sum_Matrix;                                                 // pivot is the value of current pivot.
                                                                                    // buff is the buffer for the (double) value.
                                                                                    //sum_Matrix stroe the value of sum of the matrix multiply value. 
    // Requests for memory for L_matrix, U_atrix and P_matrix
    for(i = 0; i < *size_matrix; i++){
        L_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the L_matrix[i].
        U_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the U_matrix[i].
        P_matrix[i]=(double*)malloc(sizeof(double)*(*size_matrix));                 // Requests for (size_matrix)*(double) sized addresses are stored in the P_matrix[i].
    }
    // Initialize the Lower matrix (Identity matrix)
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0;cols < *size_matrix; cols++){
            if(rows == cols){
                L_matrix[rows][cols] = 1.0;                                         // The value of the diagonal of the matrix is 1.
            }else{
                L_matrix[rows][cols] = 0.0;                                         //  The value of the matrix in other position is 0.
            }
        }
    }
    // Initialize Upper matrix
    for(rows = 0; rows < *size_matrix; rows++){
        for(cols = 0;cols < *size_matrix; cols++){
            U_matrix[rows][cols] = matrix[rows][cols];                              // Upper matrix is the same matrix of matrix we enter in.
        }
    }
    // Initialize P vector
    for(i = 0; i < *size_matrix; i++){
        P_vector[i] = i;                                                            // Assume P vector is 0, 1, 2, 3, ..., n (n = size of matrix)
    }
    // Beginning of LUP demcomposion 
    for(i = 0; i < *size_matrix; i++){
        // Locate the row that stores the largest principal. Swap the row with the largest principal element with the current row.
        pivot = U_matrix[i][i];
        max_pivot_index = i;                                                        // Assume that the current value is the max store the index
        for(j = i+1;j < *size_matrix; j++){
            if(fabs(U_matrix[j][i]) > fabs(pivot)){                                 // If the current value is not bigger than U_matrix[i][j].
                pivot = U_matrix[j][i];                                             // Convert two lines index.
                max_pivot_index = j;                                                // Store the j in max_pivot_index.
            }
        }
        // The current and destination rows of the matrix swap. 
        if(max_pivot_index != i){                                                   // When the max_pivot_index is not equal to i row means we need to change the row
            // Swap P_vector by max_pivot_index
            buff = P_vector[i];                                                     // Store current value of P_vector
            P_vector[i] = P_vector[max_pivot_index];                                // Swap the target value with current position
            P_vector[max_pivot_index] = buff;                                       // Swap the current value with target position
            // Swap U_matrix by max_pivot_index
            for(j = 0;j < *size_matrix;j++){
                buff = U_matrix[max_pivot_index][j];                                // Store current value of U matrix
                U_matrix[max_pivot_index][j] = U_matrix[i][j];                      // Swap the target value with current position
                U_matrix[i][j] = buff;                                              // Swap the current value with target position
            }
            // Swap L_matrix by max_pivot_index
            for(j = 1; (i-j) >= 0; j++){
                buff = L_matrix[max_pivot_index][i-j];                              // Store current value of L matrix
                L_matrix[max_pivot_index][i-j] = L_matrix[i][i-j];                  // Swap the target value with current position
                 L_matrix[i][i-j] = buff;                                           // Swap the current value with target position
            }
        }
        // L_matrix calculation
        for(rows = i+1; rows< *size_matrix; rows++){
            L_matrix[rows][i] = U_matrix[rows][i] / U_matrix[i][i];                 // The value below the diagonal of the L matrix is the value of each row of the U matrix except the principal element divided by the value of the principal element.
        }
        // U_matrix calculation
        for(rows = i; rows < *size_matrix; rows++){
            for(cols = i;cols < *size_matrix; cols++){
                // Gaussian elimination 
                sum_Matrix = (double)0;                                             // Initialize the value of sum_matrix. Matrix multiplication
                for(j = i; j < *size_matrix; j++){
                    if(j == rows){
                        sum_Matrix += 1.0*U_matrix[j][cols];                        // Do not change the value when j = rows.
                    }else{
                        sum_Matrix += L_matrix[rows][j]*U_matrix[j][cols]*(-1.0);   // Do elimination when j != rows.
                    }
                }
                U_matrix[rows][cols] = sum_Matrix;                                  // Calculate the value of Matrix multiplication
            }
        }
    }
    //
    convert_P_vector2matrix(P_vector,P_matrix, size_matrix);                        // Using function convert_P_vector2matrix(double* ,double** , int* ) to convert the vector.
    // Return the result of Lower matrix, Upper matrix, and P matrix.
    LUP_result[0] = L_matrix;
    LUP_result[1] = U_matrix;
    LUP_result[2] = P_matrix;
    // Free the memory of the matrixes. In case of the "memory leak" problem happ
    free(P_vector);
    // End of file with return vlaue double*** LUP_result
    return LUP_result;
}
