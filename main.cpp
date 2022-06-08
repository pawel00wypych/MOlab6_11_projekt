#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
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
const double H = 0.1;
typedef double** matrix;
typedef double* vector;

//====================================================================functions==========================================//
void KMB();
matrix provide_KMB_solution(const int r, const int c, double h);

void Laasonen_SOR();
matrix provide_Laasonen_SOR(const int r,const int c,double h);
void SOR(matrix matrixA, vector b, vector x, const int r, int const c);
vector residuum(double **macierz, double *b, double *x, int m);
double estymator(double *xNowe, double *xPoprzednie, int m);

void Laasonen_Thomas();
matrix provide_Laasonen_Thomas(const int r, const int c, double h);
void Thomas(vector Lower, vector Diagonal, vector Upper, vector b, vector x, const int c);
void compute_eta(vector Lower, vector Diagonal, vector Upper, vector l, const int c);
void compute_r(vector b, vector l, const int c);
void compute_x(vector Diagonal, vector Upper, vector b, vector x, const int c);


double provide_delta_t(const double D, const double h, const double lambda);
matrix provide_analytical_solution(const double h, const double delta_t, const int r, const int c);
matrix init_matrix_conditions(const int r, const int c, double h);
matrix get_errors_matrix(const int r, const int c,matrix analytical, matrix computational);
vector get_max_error(matrix errors, const int r, const int c);
double norm_max(vector row, int c);
vector get_x_steps(const int c);
vector get_t_steps(double delta_t, const int r);

void task2_save_data(matrix matrix_comp,matrix matrix_analytical, vector steps, const int c, const int row, string file_name, bool log_10);
void matrix_to_file(matrix mtrx, const int r, const int c, string file_name);
void vector_to_file(vector _vector, const int c, string file_name, bool log_10);
void vector_to_file(vector _vector, vector _vector2, const int c, string file_name, bool log_10);
void vector_to_file(vector _vector, vector _vector2, vector _vector3, const int c, string file_name, bool log_10);
void clear(matrix analytical, matrix computational, vector max_error, matrix errors, vector x_steps, vector t_steps, const int r);



//=========================================================main===========================================================//

int main()
{
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
    for(int k = 0; k < r; k++)
    {
        _max_error[k] = norm_max(errors[k],c);
    }
    return _max_error;
}

double norm_max(vector row, int c)
{

    double current_max = fabs(row[0]);

    for (int i = 1; i < c; i++)
    {
        if (current_max < fabs(row[i])) {
            current_max = fabs(row[i]);
        }
    }
    return current_max;
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

matrix init_matrix_conditions(const int r, const int c, double h)
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



//========================================================printing to files================================================//
void matrix_to_file(matrix _matrix, const int r, const int c, string file_name)
{
    fstream file(file_name.c_str(), ios::out);

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

void vector_to_file(vector _vector, const int c, string file_name, bool log_10)
{
    fstream file(file_name.c_str(), ios::out);

    if(log_10) {
        for (int i = 0; i < c; i++) {
            file << _vector[i] << " " << endl;
        }
    }else
    {
        for (int i = 0; i < c; i++) {
            file << _vector[i] << " " << endl;
        }
    }

    file.close();
}

void vector_to_file(vector _vector, vector _vector2, const int c, string file_name, bool log_10)
{
    fstream file(file_name.c_str(), ios::out);

    if(log_10) {
        for (int i = 0; i < c; i++) {
            file << _vector[i] << " " << _vector2[i] << endl;
        }
    }else
    {
        for (int i = 0; i < c; i++) {
            file << _vector[i] << " " << _vector2[i] << endl;
        }
    }

    file.close();
}

void vector_to_file(vector _vector, vector _vector2, vector _vector3, const int c, string file_name, bool log_10)
{
    fstream file(file_name.c_str(), ios::out);

    if(log_10) {
        for (int i = 0; i < c; i++) {
            file << log10(_vector[i]) << " " << log10(_vector2[i]) << " " << log10(_vector3[i]) << endl;
        }
    }else
    {
        for (int i = 0; i < c; i++) {
            file << _vector[i] << " " << _vector2[i] << " " <<_vector3[i] << endl;
        }
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
        x += 0.1;
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


void task2_save_data(matrix matrix_comp ,matrix matrix_analytical , vector steps, const int c, const int row, string file_name, bool log_10) {
    vector temp_comp = new double[c];
    vector temp_analytical = new double[c];
    cout<<"row "<<row<<endl;
    for (int i = 0; i < c; i++)
    {
        temp_comp[i] = matrix_comp[row][i];
        temp_analytical[i] = matrix_analytical[row][i];
    }

    vector_to_file(steps, temp_comp, temp_analytical, c, file_name, log_10);
    delete [] temp_comp;
    delete [] temp_analytical;
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
    //zad 2

    double delta_t = provide_delta_t(D, H, KMB_lambda);
    int r = ((t_max - t_min) / delta_t);
    int c = ((x_end - x_start) / H);
    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);

    matrix _analytical = provide_analytical_solution(H, delta_t, r, c);
    matrix _kmb = provide_KMB_solution(r, c, H);
    matrix errors_matrix = get_errors_matrix(r, c, _analytical, _kmb);
    vector max_error = get_max_error(errors_matrix, r, c);

    matrix_to_file(_analytical, r, c, "kmb_analytical_data.txt");
    matrix_to_file(_kmb, r, c, "kmb_computational_data.txt");
    matrix_to_file(errors_matrix, r, c, "kmb_errors_matrix_data.txt");
    vector_to_file(t_steps, r, "kmb_t_steps.txt", false);
    vector_to_file(x_steps, r, "kmb_x_steps.txt", false);

    task2_save_data(_kmb, _analytical, x_steps, c, 20, "kmb_zad2_1.txt", false);
    task2_save_data(_kmb, _analytical, x_steps, c, int(c/3), "kmb_zad2_2.txt", false);
    task2_save_data(_kmb, _analytical, x_steps, c, int(c/2), "kmb_zad2_3.txt", false);
    //zad3
    vector_to_file(t_steps,max_error, r, "kmb_zad3_err.txt", false);


    //zad 1
    int steps = 150;
    int r1 = 0;
    vector max_t_error = new double[steps];
    vector x_steps_1 = new double[steps];
    double h1 = 0.225;

    for(int i=0;i<steps;i++){
        max_t_error[i] = 0.0;
    }
    double x = h1;
    matrix errors_matrix_1;
    vector max_error_1;

    for(int j=0;j<steps;j++) {
        x_steps_1[j] = x;
        const double delta_t = provide_delta_t(D, x, KMB_lambda);
        r1 = ((t_max - t_min) / delta_t);
        c = ((x_end - x_start) / x);

        _analytical = provide_analytical_solution(x, delta_t, r1, c);
        _kmb = provide_KMB_solution(r1, c, x);

        errors_matrix_1 = get_errors_matrix(r1, c, _analytical, _kmb);
        max_error_1 = get_max_error(errors_matrix_1, r1, c);

        max_t_error[j] = max_error_1[r1-1];

        x /= 1.009;
    }
    vector_to_file(x_steps_1,max_t_error, steps, "kmb_zad1_tmax_error.txt", true);



    delete [] max_t_error;
    delete [] max_error_1;
    delete [] x_steps_1;
    for(int i = 0; i < r1; i++)
    {
        delete[] errors_matrix_1[i];
    }

    clear(_analytical, _kmb, max_error, errors_matrix, x_steps, t_steps, r);

}

matrix provide_KMB_solution(const int r, const int c, double h)
{
    matrix _matrix = init_matrix_conditions(r, c, h);//inicjalizujemy macierz z warunkami początkowymi oraz brzegowymi

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
    //zad 2

    double delta_t = provide_delta_t(D, H, LAASONEN_lambda);
    int r = ((t_max - t_min) / delta_t);
    int c = ((x_end - x_start) / H);
    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);

    matrix _analytical = provide_analytical_solution(H, delta_t, r, c);
    matrix _laasonen_thomas = provide_Laasonen_Thomas(r, c, H);
    matrix errors_matrix = get_errors_matrix(r, c, _analytical, _laasonen_thomas);
    vector max_error = get_max_error(errors_matrix, r, c);

    matrix_to_file(_analytical, r, c, "lth_analytical_data.txt");
    matrix_to_file(_laasonen_thomas, r, c, "lth_computational_data.txt");
    matrix_to_file(errors_matrix, r, c, "lth_errors_matrix_data.txt");
    vector_to_file(t_steps, r, "lth_t_steps.txt", false);
    vector_to_file(x_steps, r, "lth_x_steps.txt", false);

    task2_save_data(_laasonen_thomas, _analytical, x_steps, c, 20, "lth_zad2_1.txt", false);
    task2_save_data(_laasonen_thomas, _analytical, x_steps, c, int(c/3), "lth_zad2_2.txt", false);
    task2_save_data(_laasonen_thomas, _analytical, x_steps, c, int(c/2), "lth_zad2_3.txt", false);
    //zad3
    vector_to_file(t_steps,max_error, r, "lth_zad3_err.txt", false);


    //zad 1
    int steps = 150;
    int r1 = 0;
    vector max_t_error = new double[steps];
    vector x_steps_1 = new double[steps];
    double h1 = 0.225;

    for(int i=0;i<steps;i++){
        max_t_error[i] = 0.0;
    }
    double x = h1;
    matrix errors_matrix_1;
    vector max_error_1;

    for(int j=0;j<steps;j++) {
        x_steps_1[j] = x;
        delta_t = provide_delta_t(D, x, LAASONEN_lambda);
        r1 = ((t_max - t_min) / delta_t);
        c = ((x_end - x_start) / x);

        _analytical = provide_analytical_solution(x, delta_t, r1, c);
        _laasonen_thomas = provide_Laasonen_Thomas(r1, c,x);

        errors_matrix_1 = get_errors_matrix(r1, c, _analytical, _laasonen_thomas);
        max_error_1 = get_max_error(errors_matrix_1, r1, c);

        max_t_error[j] = max_error_1[r1-1];

        x /= 1.009;
    }
    vector_to_file(x_steps_1,max_t_error, steps, "lth_zad1_tmax_error.txt", true);



    delete [] max_t_error;
    delete [] max_error_1;
    delete [] x_steps_1;
    for(int i = 0; i < r1; i++)
    {
        delete[] errors_matrix_1[i];
    }

    clear(_analytical, _laasonen_thomas, max_error, errors_matrix, x_steps, t_steps, r);
}

matrix provide_Laasonen_Thomas(const int r, const int c, double h)
{
    const double _lambda = 1.0 + 2.0 * LAASONEN_lambda;//PDF str.136 macierz trójdiagonalna
    matrix _matrix_A = init_matrix_conditions(r, c, h);

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

    //zad 2

    double delta_t = provide_delta_t(D, H, LAASONEN_lambda);
    int r = ((t_max - t_min) / delta_t);
    int c = r ;// ((x_end - x_start) / h);
    vector t_steps = get_t_steps(delta_t, r);
    vector x_steps = get_x_steps(c);

    matrix _analytical = provide_analytical_solution(H, delta_t, r, c);
    matrix _laasonen_sor = provide_Laasonen_SOR(r, c, H);
    matrix errors_matrix = get_errors_matrix(r, c, _analytical, _laasonen_sor);
    vector max_error = get_max_error(errors_matrix, r, c);

    matrix_to_file(_analytical, r, c, "lsor_analytical_data.txt");
    matrix_to_file(_laasonen_sor, r, c, "lsor_computational_data.txt");
    matrix_to_file(errors_matrix, r, c, "lsor_errors_matrix_data.txt");
    vector_to_file(t_steps, r, "lsor_t_steps.txt", false);
    vector_to_file(x_steps, r, "lsor_x_steps.txt", false);

    task2_save_data(_laasonen_sor, _analytical, x_steps, c, 20, "lsor_zad2_1.txt", false);
    task2_save_data(_laasonen_sor, _analytical, x_steps, c, int(c/3), "lsor_zad2_2.txt", false);
    task2_save_data(_laasonen_sor, _analytical, x_steps, c, int(c/2), "lsor_zad2_3.txt", false);
    //zad3
    vector_to_file(t_steps,max_error, r, "lsor_zad3_err.txt", false);


    //zad 1
    int steps = 130;
    int r1 = 0;
    vector max_t_error = new double[steps];
    vector x_steps_1 = new double[steps];
    double h1 = 0.35;

    for(int i=0;i<steps;i++){
        max_t_error[i] = 0.0;
    }
    double x = h1;
    matrix errors_matrix_1;
    vector max_error_1;

    for(int j=0;j<steps;j++) {
        x_steps_1[j] = x;
        delta_t = provide_delta_t(D, x, LAASONEN_lambda);
        r1 = ((t_max - t_min) / delta_t);
        c = r1;//((x_end - x_start) / x);
        cout<<"j "<<j<<endl;
        _analytical = provide_analytical_solution(x, delta_t, r1, c);
        _laasonen_sor = provide_Laasonen_SOR(r1, c, x);

        errors_matrix_1 = get_errors_matrix(r1, c, _analytical, _laasonen_sor);
        max_error_1 = get_max_error(errors_matrix_1, r1, c);

        max_t_error[j] = max_error_1[r1-1];

        x /= 1.009;
    }
    vector_to_file(x_steps_1,max_t_error, steps, "lsor_zad1_tmax_error.txt", false);



    delete [] max_t_error;
    delete [] max_error_1;
    delete [] x_steps_1;
    for(int i = 0; i < r1; i++)
    {
        delete[] errors_matrix_1[i];
    }

    clear(_analytical, _laasonen_sor, max_error, errors_matrix, x_steps, t_steps, r1);

}

matrix provide_Laasonen_SOR(const int r,const int c, double h)
{


    const double _lambda = 1.0 + 2.0 * LAASONEN_lambda;//PDF str.136 macierz trójdiagonalna
    matrix _matrix_A = init_matrix_conditions(r, c, h);

    vector vector_b = new double[c];
    vector vector_x = new double[c];
    matrix matrix_temp = new vector[c];
    for (int i = 0; i < c; i++) {
        matrix_temp[i] = new double[c];
    }

    for(int i = 0; i < c; i++)
    {
        vector_b[i] = 0.0;
        vector_x[i] = 0.0;
    }

    for (int k = 1; k < r; k++) {

        matrix_temp[0][0] = 1.0;//Diagonal[0] = 1.0;
        vector_b[0] = _matrix_A[k-1][0];

        for (int i = 1; i < c - 1; i++) {
            matrix_temp[i+1][i] = LAASONEN_lambda;// Lower[i]    = LAASONEN_lambda;
            matrix_temp[i][i] = -_lambda;//Diagonal[i] = -_lambda;
            matrix_temp[i-1][i] = LAASONEN_lambda;//Upper[i]    = LAASONEN_lambda;
            vector_b[i] = -_matrix_A[k - 1][i];
        }

        matrix_temp[r-1][c-1] = 1.0;//Diagonal[c - 1] = 1.0;
        vector_b[c - 1] = 0.0;


        SOR( matrix_temp, vector_b, vector_x, r, c);

        for (int i = 1; i < c-1; i++) {
            _matrix_A[k][i] = vector_x[i];
        }
    }

    return _matrix_A;
}

void SOR(matrix  matrix_A, vector b, vector x, const int r, int const c)
{

    double TOL = 1e-16;
    double *xnplus1 = new double[c];

    for( int i = 0; i < c; i++)
   {
       xnplus1[i] = 0.0;
   }

    double suma;
    double OMEGA = 0.5;

    for (int iter = 0; iter < 40; iter++)
    {
        for (int i = 1; i < c; i++)
        {
            suma = 0.0;
            for (int j = 1; j < c; j++)
            {
                    if (i != j) {
                        suma += (matrix_A[i][j] * x[j]);
                    }
            }
            xnplus1[i]  = (1.0 - OMEGA) * x[i] + ((OMEGA * (b[i] - suma)) / matrix_A[i][i]);
            x[i] = xnplus1[i];
        }

        if ((fabs((norm_max(residuum(matrix_A, b, x, c), c))) < TOL) &&
            (fabs((estymator(xnplus1, x, c))) < TOL))
                break;
    }
}

double estymator(vector xNowe, vector xPoprzednie, int m) {
    double max = 0.0;
    double *p = new double[m];

    for (int i = 0; i < m; i++)
        p[i] = xNowe[i] - xPoprzednie[i];

    if (fabs(p[0]) > fabs(p[1]))
        max = fabs(p[0]);
    else
        max = fabs(p[1]);

    for (int i = 0; i < m; i++) {
        if (fabs(p[i]) > max)
            max = fabs(p[i]);
    }

    delete[] p;
    return max;
}

vector residuum(matrix macierz, vector b, vector x, int c) {
    double sum = 0.0;
    double *wynik = new double[c];
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++) {
            sum += macierz[i][j] * x[j];
        }
        wynik[i] = sum - b[i];
        sum = 0.0;
    }
    return wynik;
}