#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

void cholesky(double**, int);
double calc_norm(double**, int);
int factorial(int);
void substract_matrices(double**, double**, int, double**);
void calc_lt_inv(double**, int);
void calc_up_inv(double**, int);
void print_inv(double**, int);
void find_inv(double**, int, double**);
void substract_eye(double**, int, double**);
void print_inv(double**, int);
void mult_matr(double** , double**, int, double**);

int main(int argc, char* argv[]){
    int n;//dimension of matrix 
    double** matrix, **matrix_og;
    double** inv_expected;
    double** inv_calc;

    if(atoi(argv[1]) == 1){//otan o xristis thelei na eisagei monos tou ta dedomena, eisagei kai ton pinaka A kai to dianisma b
        printf("please enter the data of the matrix you want to test\n");
        printf("Dimension of matrix: ");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n) * sizeof(double));
        inv_calc = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_calc[i] = (double*)malloc((n) * sizeof(double));
        inv_expected  = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_expected[i] = (double*)malloc((n) * sizeof(double));
        printf("enter the elements of the matrix row by row\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                scanf("%lf", &matrix[i][j]);
                matrix_og[i][j] = matrix[i][j];
            }
        }
        printf("the matrix you entered is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("enter the elements of the expected inverted matrix row by row\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                scanf("%lf", &inv_expected[i][j]);
            }
        }
    }


    //stin periptosi tou tixaiou pinaka xrisimopoio x=(1,1,...,1)^T kai b=Ax
    else if(atoi(argv[1]) == 2){//otan theloume tixaio pinaka, gia eukoleia stis prakseis evala anotato orio to 100
        printf("a random matrix will be created\n");
        printf("please enter the dimension of the matrix: ");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n) * sizeof(double));
        inv_calc = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_calc[i] = (double*)malloc((n) * sizeof(double));
        inv_expected  = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_expected[i] = (double*)malloc((n) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][j] = rand()%100;
                matrix_og[i][j] = matrix[i][j];
            }
        } 
        
        for(int i=0; i<n; i++){
            matrix[i][n] = 0.0;
            matrix_og[i][n] = 0.0;
        }
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
    }

    else if(atoi(argv[1]) == 3){//otan theloume na diavasei apo arxeio, sto arxeio prepei na nai se grammes ksexorista ta stoixeia
        char* X;
        X = argv[2];
        printf("reading the elements of the matrix from given file %s\n", X);
        FILE* fptr;
        fptr = fopen(X,"r");
        if (fptr == NULL){
            printf("Error opening file"); 
            return 1;
        }         
        fscanf(fptr, "%d", &n);
        printf("dimension of matrix is: %d\n", n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n) * sizeof(double));
        inv_calc = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_calc[i] = (double*)malloc((n) * sizeof(double));
        inv_expected  = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_expected[i] = (double*)malloc((n) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                fscanf(fptr, "%lf", &matrix[i][j]);
                matrix_og[i][j] = matrix[i][j];
            }
        }
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                fscanf(fptr, "%lf", &inv_expected[i][j]);
            }
        }
        printf("the matrix that was found in the file is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("the expected inverse of the matrix is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", inv_expected[i][j]);
            }
            printf("\n");
        }
        fclose(fptr);
    }

    else if(atoi(argv[1]) == 4){//endeiktiki efarmogoi 1.b.3
        printf("creating matrix from 1.b.3\n choose dimension of matrix: ");
        scanf("%d", &n);
        printf("implement cholesky decomposition on hilbert matrix with n=%d\n", n);
        printf("dimension of hilbert matrix is: %d\n", n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n+1) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n+1) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][j] = (double)1.0/(i+j+1.0);
                matrix_og[i][j] = matrix[i][j];
            }
        }
        inv_calc = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_calc[i] = (double*)malloc((n) * sizeof(double));
        inv_expected  = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_expected[i] = (double*)malloc((n) * sizeof(double));
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                int a, b, c, d, e, f, g, h, k, l, m;
                a = factorial(n+j);
                b = factorial(n+i);
                c = pow(-1, i+j+2);
                d = a*b*c;
                e = factorial(i);
                f = factorial(j);
                g = pow(e*f, 2);
                h = factorial(n-i-1);
                k = factorial(n-j-1);
                l = (i+j+1);
                m = l*g*h*k;
                inv_expected[i][j] = (double)d/m;
            }
        }
        printf("the expected inverted matrix is:\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", inv_expected[i][j]);
            }
            printf("\n");
        }
    }
    else if(atoi(argv[1]) == 5){//endeiktiki efarmogi 1.b.2
        printf("creating matrix from 1.b.2\n choose dimension of matrix: ");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n) * sizeof(double));
        inv_calc = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_calc[i] = (double*)malloc((n) * sizeof(double));
        inv_expected  = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            inv_expected[i] = (double*)malloc((n) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(i!=j){
                    matrix[i][j] = 1.0;
                }
                else{
                    matrix[i][j] = n;
                }
                matrix_og[i][j] = matrix[i][j];
            }
        }
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }

        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(i!=j){
                    inv_expected[i][j] = (double)-1.0/((n-1)*(2*n-1));
                }
                else{
                    inv_expected[i][j] = (double)2/(2*n-1);
                }
            }
        }
        printf("the expected inverted matrix is:\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", inv_expected[i][j]);
            }
            printf("\n");
        }
    }
    

    printf("Calculating inverse of matrix with Cholesky method\n");
    clock_t t;
    t = clock();

    //cholesky decomposition 
    cholesky(matrix, n);

    //evresi tou antistrofou
    calc_lt_inv(matrix, n);
    find_inv(matrix, n, inv_calc);
    print_inv(inv_calc, n);



    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("the algorithm took %f seconds to execute\n", time_taken);

    //evresi condition number
    double norm = calc_norm(matrix_og,n);
    printf("the norm of given matrix is: %lf\n", norm);
    double norm_inv = calc_norm(inv_expected, n);
    printf("the norm of expected inverse matrix is: %lf\n", norm_inv);
    double cond = norm*norm_inv;
    printf("the condition number of the matrix is: %lf\n", cond);

    //ipologismos sfalmaton
    printf("Calculating errors\n");
    double error1, error2, tmp1, tmp2;
    double** diff;
    diff = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        diff[i] = (double*)malloc(1 * sizeof(double));
    substract_matrices(inv_expected, inv_calc, n, diff);
    tmp1 = calc_norm(diff, n);
    tmp2 = calc_norm(inv_expected, n);
    error1 = (double)tmp1/tmp2;
    printf("||A^-1-Acal^-1||/||A^-1|| = %e\n", error1);
    mult_matr(matrix_og, inv_calc, n, diff);
    substract_eye(diff, n, diff);
    tmp1 = calc_norm(diff, n);
    error2 = (double)tmp1/tmp2;
    printf("||AAcal^-1-I||/||A^-1|| = %e\n", error2);

    //apodesmeusi xorou
    for(int i = 0; i < n; i++)
        free(diff[i]);
    free(diff);
    for(int i = 0; i <n; i++)
        free(matrix[i]);
    free(matrix);
    for(int i = 0; i <n; i++)
        free(matrix_og[i]);
    free(matrix_og);
    for(int i = 0; i <n; i++)
        free(inv_expected[i]);
    free(inv_expected);
    for(int i = 0; i <n; i++)
        free(inv_calc[i]);
    free(inv_calc);
    return 0;
}


void cholesky(double** matrix, int n){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            if (j == i){
                for (int k = 0; k < j; k++){
                    sum += pow(matrix[j][k], 2);
                }
                matrix[j][j] = sqrt(matrix[j][j] - sum);
            }
            else{
                for (int k = 0; k < j; k++){
                    sum += (matrix[i][k] * matrix[j][k]);
                }
                matrix[i][j] = (matrix[i][j] - sum) /matrix[j][j];
                matrix[j][i] = 0.0;
            }
        }
    }
    printf("the L matrix that was calculated based on Cholesky method is: \n");//gia eksoikonomisi xorou oi pinakes L,L^T einai sto idio pinaka
    for(int q=0;q<n;q++){
        for(int l=0;l<n;l++){
            printf("%lf\t ",matrix[q][l]);
        }
        printf("\n");
    }
}
void find_inv(double** matrix, int n, double** inv_calc){
    double** temp;
    temp = (double**)malloc((n) * sizeof(double*));
    for (int i = 0; i < n; i++)
        temp[i] = (double*)malloc((n) * sizeof(double));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv_calc[i][j] = 0.0;
            temp[i][j] = matrix[j][i];

      }
    }
    for (int i = 0; i <n; ++i){
        for (int j = 0; j <n; ++j){ 
            for (int k = 0; k <n; ++k) {
                inv_calc[i][j] += (double)temp[i][k] * matrix[k][j];
            }
        }
    }
    for(int i = 0; i < n; i++)
        free(temp[i]);
    free(temp);
}


double calc_norm(double** matrix, int n){
    double norm[n], max;
    for(int i=0; i<n; i++){
        norm[i]=0.0;
    }
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<n; j++) {
            norm[i] += fabs(matrix[i][j]);
        }
    }
    max = norm[0];
    for(int i = 0; i<n; i++) {
        if (max < norm[i]) {
            max = norm[i];
        }
    }
    return max;
}

int factorial(int n){
    int i,f=1;
    for(i=1;i<=n;i++){
        f=f*i;
    }
    return f;
}

void substract_matrices(double** first, double** second, int n, double** diff){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            diff[i][j] = (double) first[i][j] - second[i][j];
        }
    }
}

//ipologizei me forward substitution ton antistrofo tou L
void calc_lt_inv(double** matrix, int n){
    double eye[n];
    double temp[n];
    for(int i=0; i<n; i++){
        eye[i] = 0.0;
    }
    for(int cnt=0; cnt<n; cnt++){
        eye[cnt-1] = 0.0;
        eye[cnt] = 1.0;
        for(int i=0; i<n; i++){
            temp[i] = 0.0;
        }
        for(int i=0;i<n;i++){       //pros ta empros
            double sum=0.0;
            if(i == 0){
                temp[i] = eye[i]/matrix[i][i];
                if(cnt<=i){//
                    matrix[i][cnt] = temp[i];
                }
            }
            else{
                for(int j=0;j<=i;j++){
                    sum+=(double)matrix[i][j]*temp[j];
                }
                if(cnt<=i){///
                    temp[i]=(double)(eye[i]-sum)/matrix[i][i];
                    matrix[i][cnt] = temp[i];
                }
            }
        }
    } 
}

void print_inv(double** matrix, int n){
    printf("\nthe calculated inverse matrix is\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void mult_matr(double** matrix, double** inv, int n, double** diff){
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
         diff[i][j] = 0.0;
      }
    }
    for (int i = 0; i <n; ++i){
        for (int j = 0; j <n; ++j){ 
            for (int k = 0; k <n; ++k) {
                diff[i][j] += (double)matrix[i][k] * inv[k][j];
            }
        }
    }
}

void substract_eye(double** first, int n, double** diff){
    double eye[n][n];
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j)
                eye[i][j] = 1;
            else    
                eye[i][j] = 0;
        }
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            diff[i][j] = (double) first[i][j] - eye[i][j];
        }
    }
}