//**********Project 1**********//
// Name: Yuxuan Zhang
// UHiD: 2137717
// Date: 04,03,2023
// Last Update: 04,11,2023

//**********WARNING**********//
// This programme using the method of "malloc" to access memory.
// The return value of the function in this program is the address. 
// After using it, please use "free (address)" to avoid memory leaks (in the main function).

//**********HEAD FILES**********//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//**********CONSTANT**********//
#define machine_epsilon 5.96e-8 // The machine error of float type
#define TRUE 1                   // The TRUE is 1

//**********DECLARATION**********//
void horner(double* C, int N, double x, double* Q);
double* newton(double* C, int N, double* B);
double* bisection(double* C, int N, double* B);

//**********MAIN FUNCTION**********//
void main(int argc, char** argv){
    int i;
    double* Boundary = (double*)malloc((2)*sizeof(double));    // Here the Boundary stores the boudary by using malloc function to get space
    double* Coefficient = (double*)malloc((argc-3)*sizeof(double));    //  Here the Coefficient stores the coeffficient by using malloc function to get space
    // If the number enter is smaller than 3 means that the coeficient is not enough we need to stop
    if(argc < 3){
        return;        
    }
    // Convert the numbers which are got from argv to Coefficient and Boundary array
    for(i = 1;i <= argc-1;i++){
        if(i <= argc-3){
            Coefficient[i-1] = atof(argv[i]);   // Write the value from argv to Coefficient without last two value
        } else{
            Boundary[i-argc+2] = atof(argv[i]); // From Boundary[0] to Boundary[1]
        }
    }
    double* result = newton(Coefficient,argc-3,Boundary);    // Using newton method first to get the zero point
    // Judge that whether the newton method is working (if it works result[2] = 1)
    if(result[2] == 1){
        printf("%lf %d",result[0],(int)result[1]);  // Print out the value of zero point xr and interations
    }
    // Judge that whether the newton method is working (if it not works result[2] = 0)
    if(result[2] == 0){
        result = bisection(Coefficient,argc-3,Boundary);    //Using method of bisection to get the value of zero points and iterations 
        printf("%lf %d",result[0],(int)result[1]);  // Print out the value of zero point xr and interations
    }
    // Delete the dynamic memory requested by the functions and main function
    free(result);
    free(Boundary);
    free(Coefficient);
    return;
}

//**********NEWTON METHOD**********//
// Netwon method to find the value of the zero point
double* newton(double* C, int N, double* B){           // C is the coefficient. N is the number of coefficient. B is the boundary around the zero point
    double* f_x = (double*)malloc(1*sizeof(double));   // The value of f(x)
    double* ff_x = (double*)malloc(1*sizeof(double));  // The value fo f'(x)
    double* x = (double*)malloc(1*sizeof(double));     // Store the value of the x
    double* error_before = (double*)malloc(1*sizeof(double));     // The value store the error where the x[n] - x[n-1]
    double* Q = (double*)malloc(N*sizeof(double));     // The array to store the value of f(X) and f'(x) (In horner method it is called q(x))
    double* error = (double*)malloc(1*sizeof(double)); // Store the value of x[n+1] - x[n]
    double* result = (double*)malloc(3*sizeof(double)); // Store the return address. where result[0] is the zero point (x), result[1] is the interation, and result[2] is a index judge if the Newton method is converge. (0 is diverge, 1 is converge)
    int Iter = 0;                       // Count the number of interations
    *x = (B[0] + B[1])/2;               // Initialize the value of x
    horner(C,N,*x,Q);                   // Using the method of horner to get the value of f(x)
    *f_x = Q[N-1];                      // Calculate the value of f(x)
    // If the value of f(x) value at zero point (x,f(x)) is bigger than the error we use the method of Newton
    while(TRUE && fabs(*f_x)>0){    
        horner(C,N,*x,Q);               // Using the method of horner to get the value of f(x)
        *f_x = Q[N-1];                  // Store the value into address of f_x
        horner(Q,N-1,*x,Q);             // Using the method of horner to get the value of f'(x)
        *ff_x = Q[N-2];                 // Store the value into address of ff_x
        *error_before = *error;         // Remember the last error for comparing the new one to judge whether the method is converge
        *error = fabs(*f_x / *ff_x);    // Store the new error (x[n]-x[n-1])
        *x = *x - *f_x / *ff_x;         // X[n+1] = X[n] - f(X[n])/f'(X[n])
        Iter += 1;                      // Everytime we using Newton method we pulse the interation
        // If the previous error(|x[n] - x[n+1]|) is smaller than present error, which means Newton method is not converge
        if(*error_before <= *error && Iter>1){
            result[1] = Iter;           // Return the Interation to result
            result[2] = 0;              // The index of result[2] is zero diverge happened
            return result;              // End of the function back address result to the main function
        }
        //If the x is out of boudary means the Newton method is not converge
        if(*x < B[0] | *x > B[1]){
            result[1] = Iter;           // Return the Interation to result
            result[2] = 0;              // The index of result[2] is zero diverge happened
            return result;              // End of the function back address result to the main function
        }
        // If the space of the two x (x[n] and x[n-1]) is diverge
        if(fabs(*error) <= machine_epsilon && Iter>1){
            break;
        }
    }
    // Store the values 
    result[0] = *x;
    result[1] = Iter;
    result[2] = 1;
    // Delete the dynamic memory requested by the function
    free(f_x);
    free(ff_x);
    free(x);
    free(Q);
    free(error);
    free(error_before);
    // Return the address of the result which include the value of x and Iteration times
    return result;  
}

//**********BISECTION METHOD**********//
// Bisection method to find the value of the zero point
double* bisection(double* C, int N, double* B){         // C is the coefficient. N is the number of coefficient. B is the boundary around the zero point     
    double* f_x = (double*)malloc(3*sizeof(double));    // f[0] means the left boundary value of f(x), f[1] means the right boundary value of f(x) and f[2] means zero point value of f(x)
    double* x = (double*)malloc(1*sizeof(double));      // Store the value of the x
    double* Q = (double*)malloc(N*sizeof(double));      // The array to store the value of f(X) and f'(x) (In horner method it is called q(x))
    double* result = (double*)malloc(3*sizeof(double)); /*Store the return address. where result[0] is the zero point (x), result[1] is the interation, and result[2] is a index judge if the Newton method is converge. (0 is diverge, 1 is converge)*/  
    int Iter = 0;                                       // Count the number of interations
    f_x[2] = 1;                                         // Intial the value at zero point in order to get into the loop. Here we don't need to calculate the true value of f_x[2] is because Newton method has calculated.
    // If the value of f(x) value at zero point (x,f(x)) is bigger than the error we use the method of Newton
    while(B[1]-B[0] > machine_epsilon){
        *x = B[0];                                      // Store the left bound value of x in order to calculate the value of f(a)
        horner(C,N,*x,Q);                               // Using the method of horner to get the value of f(x)
        f_x[0] = Q[N-1];                                // Store the value into of f_x[0] (left side)
        *x = B[1];                                      // Store the right bound value of x in order to calculate the value of f(b)
        horner(C,N,*x,Q);                               // Using the method of horner to get the value of f(x)
        f_x[1] = Q[N-1];                                // Store the value into of f_x[1] (right side)
        Iter += 1;                                      // Everytime we using Bisection method we pulse the interation
        // Judge if the value meet the condition of
        if(f_x[0]*f_x[1] < 0){
            *x = (B[0]+B[1])/2;                         // Claculate the mid value from [a,b]
            horner(C,N,*x,Q);                           // Using the method of horner to get the value of f(x)
            f_x[2] = Q[N-1];                            // Store the value into f_x[2]
        }else{
            return result;                              // End of the function back address result to the main function
        }
        // Judge if the sign of f_x[2] is not same with the value of f(a) (Left value of boundary)
        if(f_x[0]*f_x[2] < 0){
            B[1] = *x;                                  // Store new left value b = *x
        }else{
            B[0] = *x;                                  // Store new right value a = *x
        }
    }
    // Store the values 
    result[0] = *x;
    result[1] = Iter;
    free(f_x);
    free(Q);
    free(x);
    // Return the address of the result which include the value of x and Iteration times
    return result;
}

//**********HORNOR METHOD**********//
void horner(double* C, int N, double x, double* Q){     // C is the coefficient. N is the number of coefficient. x is the point of p(x). Q store the value of the horner 
    int i = 0;                                          // Calculate the coefficient of p(x) (From high to low) to get the array of q(x)
    for(i;i<=N-1;i++){                                  // The i is from 0 to N-1 (because the total number is N)
        Q[i] = Q[i-1]*x+C[i];                           // Every last reslut times x pluse the coefficient of p(x)
    }   
    // Result: the value is store in the q(x) (double* Q) Use the rule of Hornor ways to calculate the p(x) and p'(x)
}
