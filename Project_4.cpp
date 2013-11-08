#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <iostream>
using namespace std;

/* This program solves the diffusion equation using both the explicit forward Euler,
 the implicit backward Euler and the implicit Crank-Nicolson scheme.*/

void backward_Euler(double a_i, double *b_i, double c_i, double *V_j_BE_old, double *V_j_BE_new, int j);
void forward_Euler(double a_i, double *b_i, double c_i, double *V_j_FE_old, double *V_j_FE_new, int j);
void Crank_Nicolson(double a_i, double *b_i, double c_i, double *V_j_CN_old, double *V_j_CN_mid, double* V_j_BE_new_CN__new, int j, double alpha);

int main(int argc, char* argv[]) {

    /*fstream outFile;
    outFile.open("data.dat", ios::out);
    outFile.close();*/


    // Defining variables.
    int j;
    //cout <<"Give an integer: ";
    //cin >> n;

    j = atoi(argv[1]);

    double delta_t, delta_x, delta_x_2;
    delta_x= 1.0/(j-1);
    delta_t = (delta_x*delta_x)/10;
    delta_x_2 = delta_x*delta_x;
    //cout << delta_x << " " << delta_t << endl;


    // Defining the matrix A
    double a_i;
    double c_i;
    double* b_i = new double[j];
    double alpha = delta_t/delta_x_2;
    //cout << alpha << endl;


    double* u_s = new double[j];
    double* V_j_BE_new = new double[j];
    double* V_j_BE_old = new double[j];
    double* V_j_FE_new = new double[j];
    double* V_j_FE_old = new double[j];
    double* V_j_CN_new = new double[j];
    double* V_j_CN_mid = new double[j];
    double* V_j_CN_old = new double[j];

    // Set initial conditions
    V_j_BE_old[0] = 0.0;
    V_j_BE_new[0] = 0.0;

    V_j_FE_new[0] = 0.0;
    V_j_FE_old[0] = 0.0;

    V_j_CN_new[0] = 0.0;
    V_j_CN_old[0] = 0.0;

    u_s[0] = 1.0;

    double x_i;
    double backward_u, forward_u, CN_u;


    for (int i=1; i < j; i++) {
        x_i = i*delta_x;
        u_s[i] = 1 - x_i;
        V_j_BE_old[i] = x_i - 1;
        V_j_FE_old[i] = V_j_BE_old[i];
        V_j_CN_old[i] = V_j_BE_old[i];
        //cout <<  << endl;
        cout << V_j_FE_old[i] + u_s[i] << endl;


    }

    double time_step;
    time_step = atoi(argv[2]);

    // Running Backward Euler

    fstream outFile;
    outFile.open("backward_euler_data.dat", ios::out);
    for (int i=0; i < time_step; i++){


        a_i = -alpha;
        c_i = -alpha;
        for (int k=0; k < j; k++) { b_i[k] = 1 + 2*alpha;}

        backward_Euler(a_i, b_i, c_i, V_j_BE_old, V_j_BE_new, j);

        //outFile << j << endl;
        for (int l=0; l < j; l++) {


            backward_u = V_j_BE_new[l] + u_s[l];

            outFile << backward_u << " ";

        }
        outFile << endl;

    }
    outFile.close();

    // Running Forward Euler

    fstream inFile;
    inFile.open("forward_euler_data.dat", ios::out);
    for (int i=0; i < time_step; i++){

        a_i = alpha;
        c_i = alpha;
        for (int k=0; k < j; k++) { b_i[k] = 1 - 2*alpha;}

        forward_Euler(a_i, b_i, c_i, V_j_FE_old, V_j_FE_new, j);

        for (int l=0; l < j; l++) {


            forward_u = V_j_FE_new[l] + u_s[l];

            inFile << forward_u << " ";
            //cout << forward_u << endl;

        }
        inFile << endl;

    }
    inFile.close();


    // Running Crank-Nicolson

    fstream newFile;
    newFile.open("crank_nicolson_data.dat", ios::out);
    for (int i=0; i < time_step; i++){

        a_i = alpha;
        c_i = alpha;
        for (int k=0; k < j; k++) { b_i[k] = 2 + 2*alpha;}

        Crank_Nicolson(a_i, b_i, c_i, V_j_CN_old, V_j_CN_mid, V_j_CN_new, j, alpha);

        for (int l=0; l < j; l++) {


            CN_u = V_j_CN_new[l] + u_s[l];

            newFile << CN_u << " ";
            cout << V_j_CN_new[l] << endl;

        }
        newFile << endl;

        for (int k=0; k < j; k++){

            V_j_CN_old[k] = V_j_CN_mid[k];
        }

    }
    newFile.close();
}

void backward_Euler(double a_i, double *b_i, double c_i, double *V_j_BE_old, double *V_j_BE_new, int j) {

    // Forward substitution
    for (int i=1; i < j-1; i++) {

        double factor;
        factor = a_i/b_i[i];

        b_i[i+1] = b_i[i+1] - c_i*factor;
        V_j_BE_old[i+1] = V_j_BE_old[i+1] - V_j_BE_old[i]*factor;
        //cout << V_j_BE_new[i] << endl;

    }

    //The first row of the backward sub is very simple
    V_j_BE_new[j-1] = V_j_BE_new[j-1]/b_i[j-1];


    // Backward substitution
    for (int i=j-2; i > 0; i--) {

        V_j_BE_new[i] = (V_j_BE_old[i] - c_i*V_j_BE_new[i+1])/b_i[i];
        //cout << V_j_BE_new[i] << endl;

    }
    //cout << endl;
    for (int i=0; i < j; i++) {

        V_j_BE_old[i] = V_j_BE_new[i];
        //cout << V_j_BE_new[i] << endl;

    }
}

void forward_Euler(double a_i, double *b_i, double c_i, double *V_j_FE_old, double *V_j_FE_new, int j) {

    //for (int i=0; i < j; i++){
    //cout << V_j_FE_old[i] << endl;
    //}
    // First row

    //V_j_FE_new[1] = b_i[1]*V_j_FE_old[1] + c_i*V_j_FE_old[2] + a_i*V_j_FE_old[0];

    // Then for the inner values of the matrix

    for (int i=1; i < j-1; i++) {

        V_j_FE_new[i] = b_i[i]*V_j_FE_old[i] + c_i*V_j_FE_old[i+1] + a_i*V_j_FE_old[i-1];
        //cout << i << endl;

    }

    // The last row

    //V_j_FE_new[j-1] = a_i*V_j_FE_old[j-1] + b_i[j]*V_j_FE_old[j];
    V_j_FE_new[j] = 0.0;

    for (int i=0; i < j; i++) {

        V_j_FE_old[i] = V_j_FE_new[i];

        //cout << V_j_FE_new[i] << endl;

    }

}


void Crank_Nicolson(double a_i, double *b_i, double c_i, double *V_j_CN_old, double *V_j_CN_mid, double *V_j_CN_new, int j, double alpha){

    // First we do a matrix-vector multiplication, as in Forward Euler
    // To find the vector V_j_CN_mid, to use in the tridiagonal solver

    for (int i=1; i < j-1; i++) {

        V_j_CN_mid[i] = b_i[i]*V_j_CN_old[i] + c_i*V_j_CN_old[i+1] + a_i*V_j_CN_old[i-1];
        //cout << V_j_CN_new[i] << endl;

    }

    // Then we find the vector, solution, using Backward Euler

    //double a_i, c_i;
    a_i = -alpha;
    c_i = -alpha;
    for (int k=0; k < j; k++) { b_i[k] = 2 - 2*alpha;}

    backward_Euler(a_i, b_i, c_i, V_j_CN_mid, V_j_CN_new, j);

}





