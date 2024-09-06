#include <stdio.h>
#include<stdlib.h>
#include <math.h>

void Filling_of_diagonals(int n, double *low_d, double *main_d, double *upper_d, double *x, double *y, double *free_ch){
    double* h = (double*)calloc(n, sizeof(double));
    for (int i = 1; i <= n-1; i++){
        h[i] = x[i] - x[i-1];
    }
    for (int i = 0; i < n-2; i++){
        if (i != 0 && i != n-3){
            low_d[i] = h[i+1] / 6;
            main_d[i] = (h[i+1]+h[i+1]) / 3;
            upper_d[i] = h[i+2] / 6;
        }
        else if (i == 0){
            main_d[i] = (h[i+1]+h[i+2]) / 3;
            upper_d[i] = h[i+2] / 6;
        }
        else {
            low_d[i] = h[i+1] / 6;
            main_d[i] = (h[i+1]+h[i+2]) / 3;
        }
        free_ch[i] = ((y[i+2]-y[i+1])/h[i+2]) - ((y[i+1]-y[i])/h[i+1]);

    }
    free(h);
}


void Finding_gamma(int n, double *low_d, double *main_d, double *upper_d, double *free, double *gamma)
{

    double m;
    for (int i = 1; i < n-1; i++)
    {
        m = low_d[i]/main_d[i-1];
        main_d[i] = main_d[i] - m*upper_d[i-1];
        free[i] = free[i] - m*free[i-1];
    }

    gamma[n-1] = free[n-1]/main_d[n-1];

    for (int i = n - 2; i >= 0; i--)
    {
        gamma[i]=(free[i]-upper_d[i]*gamma[i+1])/main_d[i];
    }
}
double Spline(int n, const double *gamma, const double *x, const double *y, double f_x){
    double* h = (double*)calloc(n, sizeof(double));
    for (int i = 1; i <= n-1; i++){
        h[i] = x[i] - x[i-1];
    }

    for (int i = 0; i < n-1; i++){
        if (f_x >= x[i] && f_x <= x[i+1]){
            return y[i]*((x[i+1]-f_x)/h[i+1]) + y[i+1]*((f_x - x[i])/h[i+1]) + gamma[i]*((pow((x[i+1]-f_x), 3) - pow(h[i+1], 2)*(x[i+1] - f_x))/(6*h[i+1])) + gamma[i+1]*((pow((f_x - x[i]), 3) - pow(h[i+1], 2)*(f_x-x[i]))/(6*h[i+1]));
        }
    }
    free(h);



}





void Joint_spline(double * x1, double * x2, double * y1, double * y2, double * new_x, double * new_y,  int n1, int n2, double * gamma1, double * gamma2) {
    int index_1 = 0;
    int index_2 = 0;
    int index_n = 0;
    if (x1[index_1] < x2[index_2]) {

        while ((x1[index_1] < x2[index_2])){
            index_1++;
        }
        if (x1[index_1] == x2[index_2]){
            new_x[index_n] = x1[index_1];
            new_y[index_n] = fabs(y1[index_1] - y2[index_2]);
            index_2++; index_n++; index_1++;
        }
        else {
            new_x[index_n] = x2[index_2];
            new_y[index_n] = y2[index_2];
            index_2++; index_n++;
        }
    }
    else if (x1[index_1] > x2[index_2]) {
        while ((x1[index_1] > x2[index_2])){
            index_2++;
        }
        if (x1[index_1] == x2[index_2]){
            new_x[index_n] = x1[index_1];
            new_y[index_n] = fabs(y1[index_1] - y2[index_2]);
            index_2++; index_n++; index_1++;
        }
        else {
            new_x[index_n] = x1[index_1];
            new_y[index_n] = y1[index_1];
            index_1++; index_n++;
        }
    }
    else {
        new_x[index_n] = x1[index_1];
        new_y[index_n] = fabs(y1[index_1] - y2[index_2]);
        index_2++; index_n++; index_1++;
    }

    while (index_2 < n2 && index_1 < n1) {
        if (x1[index_1] < x2[index_2]) {
            new_x[index_n] = x1[index_1];
            double temp_y = Spline(n2, gamma2, x2, y2, x1[index_1]);
            new_y[index_n] = fabs(y1[index_1] - temp_y);
            index_1++; index_n++;
        }
        else if (x1[index_1] > x2[index_2]) {
            new_x[index_n] = x2[index_2];
            double temp_y = Spline(n1, gamma1, x1, y1, x2[index_2]);
            new_y[index_n] = fabs(temp_y - y2[index_2]);
            index_2++; index_n++;
        }
        else {
            new_x[index_n] = x1[index_1];
            new_y[index_n] = fabs(y1[index_1] - y2[index_2]);
            index_2++; index_n++; index_1++;
        }
    }


}

void Print_spline(int n, double * x, double *y, double * gamma){
    double st = 0.01;
    for (int i = 0; i < n-1; i++){
        double start = x[i];
        for (int j = 0; start != x[i+1]; j++){
            start += st;
            printf("%lf\n", Spline(n, gamma, x, y, start));
        }
    }
}
int Find_cross(int n, double * x, double * y, double * gamma, double * x2, double * y2, int n2, double * gamma2){
    double e = 0.001;
    int f = 0;
    for (int i = 0; i < n-1; i++) {
        double der = (Spline(n, gamma, x, y, x[i]+e) - Spline(n, gamma, x, y, x[i]))/e;
        double xn = x[i];
        while (xn <= x[i+1] && xn >= x[i]) {
            if(fabs(Spline(n, gamma, x, y, xn)) < 0.0001) {
                printf("Coordinate of crossing: %lf %lf\n", xn, Spline(n2, gamma2, x2, y2, xn));
                f = 1;
                break;

            }
            der = (Spline(n, gamma, x, y, xn+e) - Spline(n, gamma, x, y, xn))/e;
            xn = xn - (Spline(n, gamma, x, y, xn)/der);
        }


    }
    for (int i = n-1; i > 0; i--) {
        double der = (Spline(n, gamma, x, y, x[i]) - Spline(n, gamma, x, y, x[i] - e)) / e;
        double xn = x[i];
        while (xn <= x[i] && xn >= x[i - 1]) {
            if (fabs(Spline(n, gamma, x, y, xn)) < 0.0001) {
                printf("Coordinate of crossing: %lf %lf\n", xn, Spline(n2, gamma2, x2, y2, xn));
                f = 1;
                break;

            }
            der = (Spline(n, gamma, x, y, xn) - Spline(n, gamma, x, y, xn - e)) / e;
            xn = xn - (Spline(n, gamma, x, y, xn) / der);
        }
    }
    return f;
}

void Finding_min(int n, double * x, double * y, double * gamma){
    double e = 0.001, lmd = 0.001;
    double min_d = pow(10, 5);
    for (int i = 0; i < n; i ++){
        double der = (Spline(n, gamma, x, y, x[i]+e) - Spline(n, gamma, x, y, x[i]))/e;
        double xn = x[i];
        while (fabs(lmd*der) > 0.000001) {
            xn = xn - lmd*der;
            if ((Spline(n, gamma, x, y, xn) < min_d) && (Spline(n, gamma, x, y, xn) > 0.001)) {
                min_d = Spline(n, gamma, x, y, xn);
            }
            der = (Spline(n, gamma, x, y, xn+e) - Spline(n, gamma, x, y, xn))/e;

        }

    }
    printf("Min d = %lf\n", min_d);
}



int main() {
    FILE * f1 = fopen("data_file.txt", "r");
    int n1;
    fscanf(f1, "%d\n", &n1);
    double* x1 = (double*)calloc(n1, sizeof(double));
    double* y1 = (double*)calloc(n1, sizeof(double));
    for (int i = 0; i <= 2*n1; i++){
        if (i < n1 ) {
            fscanf(f1, "%lf\n", &x1[i]);
        }
        else if (i >= n1+1){
            fscanf(f1, "%lf\n",  &y1[i-(n1+1)]);
        }
    }

    FILE * f2 = fopen("data_file_2.txt", "r");
    int n2;
    fscanf(f2, "%d\n", &n2);
    double* x2 = (double*)calloc(n2, sizeof(double));
    double* y2 = (double*)calloc(n2, sizeof(double));
    for (int i = 0; i <= 2*n2; i++){
        if (i < n2) {
            fscanf(f2, "%lf\n", &x2[i]);
        }
        else if (i >= n2+1){
            fscanf(f2, "%lf\n",  &y2[i-(n2+1)]);
        }
    }

    double* low_d1 = (double*)calloc(n1-2, sizeof(double));
    double* main_d1 = (double*)calloc(n1-2, sizeof(double));
    double* upper_d1 = (double*)calloc(n1-2, sizeof(double));
    double* vector_of_free_term1 = (double*)calloc(n1-2, sizeof(double));
    double* gamma1 = (double*)calloc(n1-2, sizeof(double));
    Filling_of_diagonals(n1, low_d1, main_d1, upper_d1, x1, y1, vector_of_free_term1);
    Finding_gamma(n1-2, low_d1, main_d1, upper_d1, vector_of_free_term1, gamma1);
    free(low_d1); free(main_d1); free(upper_d1); free(vector_of_free_term1);
    printf("\n");
    double* gamma_f1 = (double*)calloc(n1, sizeof(double));
    for (int i = 0; i < n1; i++){
        if (i == 0 || i == n1-1){
            gamma_f1[i] = 0;
        }
        else {
            gamma_f1[i] = gamma1[i-1];
        }
    }
    free(gamma1);

    double* low_d2 = (double*)calloc(n2-2, sizeof(double));
    double* main_d2 = (double*)calloc(n2-2, sizeof(double));
    double* upper_d2 = (double*)calloc(n2-2, sizeof(double));
    double* vector_of_free_term2 = (double*)calloc(n2-2, sizeof(double));
    double* gamma2 = (double*)calloc(n2-2, sizeof(double));
    Filling_of_diagonals(n2, low_d2, main_d2, upper_d2, x2, y2, vector_of_free_term2);
    Finding_gamma(n2-2, low_d2, main_d2, upper_d2, vector_of_free_term2, gamma2);
    free(low_d2); free(main_d2); free(upper_d2); free(vector_of_free_term2);
    printf("\n");
    double* gamma_f2 = (double*)calloc(n1, sizeof(double));
    for (int i = 0; i < n1; i++){
        if (i == 0 || i == n1-1){
            gamma_f2[i] = 0;
        }
        else {
            gamma_f2[i] = gamma2[i-1];
        }
    }
    free(gamma2);
    //printf("%lf", Spline(n2, gamma_f2, x2, y2));


    double* new_x = (double *)calloc(n1+n2, sizeof(double));
    double* new_y = (double *)calloc(n1+n2, sizeof(double));
    Joint_spline(x1, x2, y1, y2, new_x, new_y, n1, n2, gamma_f1, gamma_f2);
    int size_of_general_data = 0;
    for (int i = 0; i < n1+n2 - 1;  i ++){
        if(new_x[i] == 0 && new_y[i] == 0){
            if(new_x[i+1] == 0 && new_y[i+1] == 0){
                break;
            }
        }
        size_of_general_data ++;
    }
    double* low_dn = (double*)calloc(size_of_general_data-2, sizeof(double));
    double* main_dn = (double*)calloc(size_of_general_data-2, sizeof(double));
    double* upper_dn = (double*)calloc(size_of_general_data-2, sizeof(double));
    double* vector_of_free_term_n = (double*)calloc(size_of_general_data-2, sizeof(double));
    double* gamma_n = (double*)calloc(size_of_general_data-2, sizeof(double));
    Filling_of_diagonals(size_of_general_data, low_dn, main_dn, upper_dn, new_x, new_y, vector_of_free_term_n);
    Finding_gamma(size_of_general_data-2, low_dn, main_dn, upper_dn, vector_of_free_term_n, gamma_n);
    free(low_dn); free(main_dn); free(upper_dn); free(vector_of_free_term_n);
    double* gamma_fn = (double*)calloc(size_of_general_data, sizeof(double));
    for (int i = 0; i < size_of_general_data; i++){
        if (i == 0 || i == size_of_general_data-1){
            gamma_fn[i] = 0;
        }
        else {
            gamma_fn[i] = gamma_n[i-1];
        }
    }
    free(gamma_n);

    int f = Find_cross(size_of_general_data, new_x, new_y, gamma_fn, x1, y1, n1, gamma_f1);
    if (f) {
        Find_cross(size_of_general_data, new_x, new_y, gamma_fn, x1, y1, n1, gamma_f1);
    }
    else {
        Finding_min(size_of_general_data, new_x, new_y, gamma_fn);
    }
    return 0;
}
