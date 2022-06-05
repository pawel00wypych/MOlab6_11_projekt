#include <iostream>
#include <cmath>
#include <fstream>

#include "lib/calerf.h"

using namespace std;

//=======================================================constants==========================================================//
const int D = 1;
const double t_max = 2.0;
const double t_min = 0;

const double KMB_lambda = 0.4;
const double LAASONEN_lambda = 1.0;

const double a = 6.0 * sqrt(D * t_max);

const double x_start = -a;
const double x_end = a;
const double h = 0.1;
const double omega = 0.5; //Laasonen - SOR
typedef double** matrix;
typedef double* vector;
const int SOR_ITER = 50;

//====================================================================functions==========================================//
void KMB();
void Laasonen_Thomas();
void Laasonen_SOR();

matrix provide_Laasonen_SOR(const int r,const int c);
void SOR(matrix matrixA, vector b, vector x, const int r, int const c);
void swapVect(double *x1, double *x, const int c);
bool residuum(double **matrix, double *b, double *x, const int r, const int c);
bool estym(double *x, double *x1, const int c);


matrix provide_Laasonen_Thomas(const int r, const int c);
void Thomas(vector Lower, vector Diagonal, vector Upper, vector b, vector x, const int c);
void compute_eta(vector Lower, vector Diagonal, vector Upper, vector l, const int c);
void compute_r(vector b, vector l, const int c);
void compute_x(vector Diagonal, vector Upper, vector b, vector x, const int c);

matrix provide_KMB_solution(const int r, const int c);

double provide_delta_t(const double D, const double h, const double lambda);
matrix provide_analytical_solution(const double h, const double delta_t, const int r, const int c);
matrix init_matrix_conditions(const int r, const int c);
matrix get_errors_matrix(const int r, const int c,matrix analytical, matrix computational);
vector get_max_error(matrix errors, const int r, const int c);
vector get_x_steps(const int c);
vector get_t_steps(double delta_t, const int r);


void matrix_to_file(matrix mtrx, const int r, const int c, const char* file_name);
void vector_to_file(vector _vector, const int c, const char* file_name);
void clear(matrix analytical, matrix computational, vector max_error, matrix errors, vector x_steps, vector t_steps, const int r);


//=========================================================main===========================================================//
int main() {
    KMB();
    Laasonen_Thomas();
    Laasonen_SOR();
    return 0;
}



//================================================functions used in every method============================================//
//errors
matrix get_errors_matrix(const int r, const int c,matrix analytical, matrix computational)
{
    matrix _errors = new vector[r];
    for(int k = 0; k < r; k++)
    {
        _errors[k] = new double[c];
    }

    for(int k = 0; k < r; k++)
    {
        for(int i = 0; i < c; i++) {
            _errors[k][i] = fabs(analytical[k][i] - computational[k][i]);
        }
    }
    return _errors;
}

vector get_max_error(matrix errors, const int r, const int c)
{
    vector _max_error = new double[r];
    double current_max;
    for(int k = 0; k < r; k++)
    {
        current_max = fabs(errors[k][0]);
        for(int i = 1; i < c; i++)
        {
            current_max = std::max(current_max, fabs(errors[k][i]));
        }
        _max_error[k] = current_max;
    }
    return _max_error;
}


double provide_delta_t(const double D, const double h, const double lambda)
{
    return (lambda * h * h) / D;
    //ponieważ lambda = (D * delta_t) / h^2  - PDF str.134
}

matrix provide_analytical_solution(const double h, const double delta_t, const int r, const int c)
{
    matrix _matrix = new vector[r];

    for(int k=0; k < r; k++)
    {
        _matrix[k] = new double[c];
    }

    double t = t_min;
    double x = x_start;
    for(int k = 0; k < r; k++)
    {
        for(int i = 0; i < c; i++)
        {
            _matrix[k][i] = (calerf::ERFC_L(x / (2.0 * sqrt(D * t)))) / 2.0;

            x += h;
        }
        x = x_start;
        t += delta_t;
    }
    return _matrix;
}

matrix init_matrix_conditions(const int r, const int c)
{
    matrix _matrix = new vector[r];
    for(int i = 0; i < r; i++)
    {
        _matrix[i] = new double[c];
    }

    //warunki brzegowe
    for(int i = 0; i < r; i++)
    {
        _matrix[i][0] = 1.0;
        _matrix[i][c-1] = 0;
    }


    //warunek początkowy
    double x = -a;
    for(int i = 0; i < c; i++)
    {
        if(x < 0){
            _matrix[0][i] = 1.0;

        }else
            if(x >= 0)
            {
                _matrix[0][i] = 0;
            }
            x += h;
    }

    return _matrix;
}


//printing to files
void matrix_to_file(matrix _matrix, const int r, const int c, const char* file_name)
{
    fstream file(file_name, ios::out);
    for(int k = 0; k < r; k++)
    {
        for(int i = 0; i < c; i++)
        {
            file << _matrix[k][i] << " ";
        }
        file << endl;
    }
    file.close();
}

void vector_to_file(vector _vector, const int c, const char* file_name)
{
    fstream file(file_name, ios::out);

    for(int i = 0; i < c; i++)
    {
        file << _vector[i] << endl;
    }
    file.close();
}

//get steps
vector get_x_steps(const int c)
{
    vector _x_steps = new double[c];
    double x = x_start;
    for(int i=0; i < c; i++)
    {
        _x_steps[i] = x;
        x += h;
    }
    return _x_steps;
}

vector get_t_steps(double delta_t, const int r)
{
    vector _t_steps = new double[r];
    double t = t_min;

    for(int k=0; k < r; k++)
    {
        _t_steps[k] = t;
        t += delta_t;
    }
    return _t_steps;
}

//clear memory
void clear(matrix analytical, matrix computational, vector max_error, matrix errors, vector x_steps, vector t_steps, const int r)
{
    for(int i = 0; i < r; i++)
    {
        delete[] analytical[i];
        delete[] computational[i];
        delete[] errors[i];
    }

    delete[] analytical;
    delete[] computational;
    delete[] errors;
    delete[] max_error;
    delete[] x_steps;
    delete[] t_steps;

}


//======================================================KMB method========================================================//
void KMB()
{
    const double delta_t = provide_delta_t(D,h,KMB_lambda);
    const int r = ((t_max - t_min) / delta_t) + 2;
    const int c = ((x_end - x_start) / h);

    matrix _analytical = provide_analytical_solution(h, delta_t, r, c);
    matrix _kmb = provide_KMB_solution(r, c);

    matrix errors_matrix = get_errors_matrix(r ,c ,_analytical ,_kmb);
    vector max_error = get_max_error(errors_matrix, r, c);

    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);



    matrix_to_file(_analytical, r, c, "kmb_analytical.txt");
    matrix_to_file(_kmb, r, c, "kmb_computational.txt");
    matrix_to_file(errors_matrix, r, c, "kmb_errors_matrix.txt");

    vector_to_file(max_error, r, "kmb_max_error.txt");
    vector_to_file(t_steps, r, "kmb_t_steps.txt");
    vector_to_file(x_steps, r, "kmb_x_steps.txt");


    clear(_analytical, _kmb, max_error, errors_matrix, x_steps, t_steps, r);
}

matrix provide_KMB_solution(const int r, const int c)
{
    matrix _matrix = init_matrix_conditions(r, c);//inicjalizujemy macierz z warunkami początkowymi oraz brzegowymi

    for(int k = 1; k < r; k++)
    {
        for(int i = 1; i < c - 1; i++)
        {
            _matrix[k][i] = _matrix[k-1][i] + KMB_lambda * (_matrix[k-1][i-1] - (2.0 * _matrix[k-1][i]) + _matrix[k-1][i+1]);
            //ze wzoru str.134 PDF
            // U(i,k+1) = _matrix[k][i]
        }
    }
    return _matrix;
}


//============================================================Laasonen_Thomas method======================================//
void Laasonen_Thomas()
{
    const double delta_t = provide_delta_t(D,h,LAASONEN_lambda);
    const int r = ((t_max - t_min) / delta_t) + 2;
    const int c = ((x_end - x_start) / h);

    matrix _analytical = provide_analytical_solution(h, delta_t, r, c);
    matrix _laasonen_thomas = provide_Laasonen_Thomas(r, c);

    matrix errors_matrix = get_errors_matrix(r ,c ,_analytical ,_laasonen_thomas);
    vector max_error = get_max_error(errors_matrix, r, c);

    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);


    matrix_to_file(_analytical, r, c, "LTh_analytical.txt");
    matrix_to_file(_laasonen_thomas, r, c, "LTh_computational.txt");
    matrix_to_file(errors_matrix, r, c, "LTh_errors_matrix.txt");

    vector_to_file(max_error, r, "LTh_max_error.txt");
    vector_to_file(t_steps, r, "LTh_t_steps.txt");
    vector_to_file(x_steps, r, "LTh_x_steps.txt");


    clear(_analytical, _laasonen_thomas, max_error, errors_matrix, x_steps, t_steps, r);
};

matrix provide_Laasonen_Thomas(const int r, const int c)
{
    const double _lambda = 1.0 + 2.0 * LAASONEN_lambda;//PDF str.136 macierz trójdiagonalna
    matrix _matrix_A = init_matrix_conditions(r, c);

    vector Lower    = new double[c];
    vector Diagonal = new double[c];
    vector Upper    = new double[c];
    vector vector_b = new double[c];
    vector vector_x = new double[c];


    for (int k = 1; k < r; k++) {

        Lower[0]        = 0.0;
        Diagonal[0]     = 1.0;
        Upper[0]        = 0.0;
        vector_b[0]     = _matrix_A[k-1][0];

        for (int i = 1; i < c - 1; i++) {
            Lower[i]    = LAASONEN_lambda;
            Diagonal[i] = -_lambda;
            Upper[i]    = LAASONEN_lambda;
            vector_b[i] = -_matrix_A[k - 1][i];
        }

        Lower[c - 1]    = 0.0;
        Diagonal[c - 1] = 1.0;
        Upper[c - 1]    = 0.0;
        vector_b[c - 1] = _matrix_A[k - 1][c - 1];


        Thomas(Lower, Diagonal, Upper, vector_b, vector_x, c);

        for (int i = 1; i < c; i++) {
            _matrix_A[k][i] = vector_x[i];
        }
    }

    return _matrix_A;
}

void Thomas(vector Lower, vector Diagonal, vector Upper, vector b, vector x, const int c)
{
    vector vector_l = new double[c];
    //Thomas dla macierzy 3 - diagonalnej PDF str.62
    compute_eta(Lower, Diagonal, Upper, vector_l, c);

    compute_r(b, vector_l, c);

    compute_x(Diagonal, Upper, b, x, c);

    delete[] vector_l;
}

void compute_eta(vector Lower, vector Diagonal, vector Upper, vector l, const int c)
{
    //PDF str.63 algorytm thomasa
    for (int i = 1; i < c; i++)
    {
        l[i] = Lower[i] * (1 / Diagonal[i - 1]);
        Diagonal[i] = Diagonal[i] - l[i] * Upper[i-1];//diagonal = eta
    }
}

void compute_r(vector b, vector l, const int c) {//b = r
    for (int i = 1; i < c; i++) {
        b[i] = b[i] - l[i] * b[i - 1]; // ri = bi - li * (ri-1 / etai-1) bo l[i] = Lower[i] * (1/Diag[i-1]), więc
                                        // li * (ri-1 / etai-1) = Lower[i] * (1/Diag[i-1]) * b[i-1] a b[i-1] to r[i-1]
    }
}

void compute_x(vector Diagonal, vector Upper, vector b, vector x, const int c)
{
    x[c - 1] = (1 / Diagonal[c - 1] * b[c - 1]);
    for (int i = c - 2; i >= 0; i--) {
        x[i] = (1 / Diagonal[i]) * (b[i] - Upper[i] * x[i + 1]);// (1/eta[i]) * (r[i] - upper[i] *x[i+1])
    }
}


//=======================================================Laasonen - SOR method===========================================//
void Laasonen_SOR()
{
    //PDF str.75

    const double delta_t = provide_delta_t(D,h,LAASONEN_lambda);
    const int r = ((t_max - t_min) / delta_t) + 2;
    const int c = r;//((x_end - x_start) / h);
    const double omega = 2.0;

    matrix _analytical = provide_analytical_solution(h, delta_t, r, c);
    matrix _laasonen_sor = provide_Laasonen_SOR(r, c);

    matrix errors_matrix = get_errors_matrix(r ,c ,_analytical ,_laasonen_sor);
    vector max_error = get_max_error(errors_matrix, r, c);

    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);


    matrix_to_file(_analytical, r, c, "LSOR_analytical.txt");
    matrix_to_file(_laasonen_sor, r, c, "LSOR_computational.txt");
    matrix_to_file(errors_matrix, r, c, "LSOR_errors_matrix.txt");

    vector_to_file(max_error, r, "LSOR_max_error.txt");
    vector_to_file(t_steps, r, "LSOR_t_steps.txt");
    vector_to_file(x_steps, r, "LSOR_x_steps.txt");


    clear(_analytical, _laasonen_sor, max_error, errors_matrix, x_steps, t_steps, r);
}

matrix provide_Laasonen_SOR(const int r,const int c)
{
    const double _lambda = 1.0 + 2.0 * LAASONEN_lambda;//PDF str.136 macierz trójdiagonalna
    matrix _matrix_A = init_matrix_conditions(r, c);

    vector vector_b = new double[c];
    vector vector_x = new double[c];
    for(int i = 0; i < c; i++)
    {
        vector_b[i] = 0.0;
        vector_x[i] = 0.0;
    }

    for (int k = 1; k < r; k++) {

        _matrix_A[1][0] = 0.0;//Lower[0]  = 0.0;
        _matrix_A[0][0] = 1.0;//Diagonal[0] = 1.0;
        _matrix_A[0][1] = 0.0;//Upper[0] = 0.0;
        vector_b[0] = _matrix_A[k-1][0];

        for (int i = 1; i < c - 1; i++) {
            _matrix_A[i+1][i] = LAASONEN_lambda;// Lower[i]    = LAASONEN_lambda;
            _matrix_A[i][i] = -_lambda;//Diagonal[i] = -_lambda;
            _matrix_A[i-1][i] = LAASONEN_lambda;//Upper[i]    = LAASONEN_lambda;
            vector_b[i] = -_matrix_A[k - 1][i];
        }

        _matrix_A[1][c-1] = 0.0;//Lower[c-1] = 0.0
        _matrix_A[r-1][c-1] = 1.0;//Diagonal[c - 1] = 1.0;
        _matrix_A[r-2][c-1] = 0.0;//Upper[c - 1] = 0.0;

        vector_b[c - 1] = _matrix_A[k - 1][c - 1];


        SOR( _matrix_A, vector_b, vector_x, r, c);

        for (int i = 1; i < c; i++) {
            _matrix_A[k][i] = vector_x[i];
        }
    }

    return _matrix_A;
}

void SOR(matrix  matrix_A, vector b, vector x, const int r, int const c)
{
    vector x_plus_1 = new double[c];
    for( int i = 0; i < c; i++)
    {
        x_plus_1[i] = 0.0;
    }

    double sum = 0;
    for(int k=0; k<SOR_ITER; k++) {
        for (int i = 0; i < r; i++) {
            sum = 0;
            for(int j=0; j<c; j++)
            {
                if(i != j)
                {
                    sum = sum + matrix_A[i][j] * x[j];
                }
            }
            x_plus_1[i] = (1.0  - omega) * x[i] + (omega / matrix_A[i][i]) * (b[i] - sum);
        }

        bool est = estym(x, x_plus_1, c);
        bool res = residuum(matrix_A, b, x, r, c);
        if ( est && res )
            break;
        swapVect(x_plus_1, x, c);
    }
}

bool estym(double *x, double *x1, const int c)
{
    int counter = 0;
    for (int i = 0; i < c; i++) {
        if (fabs(x1[i] - x[i]) < 1e-16)
            counter++;
    }
    return counter == c;
}

bool residuum(double **matrix, double *b, double *x, const int r, const int c)
{
    int counter = 0;
    for (int i = 0; i < c; i++) {
        double temp = 0;
        for (int j = 0; j < r; j++) {
            temp += matrix[i][j] * x[j];
        }

        if (fabs(temp - b[i]) < 1e-16)
            counter++;
    }
    return counter == c;
}

void swapVect(double *x1, double *x, const int c)
{
    double tmp = 0;
    for(int i=0;i<c; i++)
    {
        tmp = x1[i];
        x1[i] = x[i];
        x[i] = tmp;
    }
}