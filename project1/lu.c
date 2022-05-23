#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

void lu_pivot(double**, int);
void solve_LUx_b(double**, int, double**);
void print_solution(double**, int);
double calc_norm(double**, int);
double calc_vector_norm(double**, int, int);
void substract_matrices(double**, double**, int, double**);
void substract_vectors(double**, double**, int, double**, int);
void mult_matrix_vector(double**, double**, int, double**);

int main(int argc, char* argv[]){
    int n;//dimension of matrix 
    double** matrix, **matrix_og;
    double** x;
    double** x_given;

    if(atoi(argv[1]) == 1){//otan o xristis thelei na eisagei monos tou ta dedomena, eisagei kai ton pinaka A kai to dianisma b
        printf("please enter the data of the matrix you want to test\n");
        printf("Dimension of matrix: ");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n+1) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n+1) * sizeof(double));
        printf("enter the elements of the matrix row by row\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                scanf("%lf", &matrix[i][j]);
                matrix_og[i][j] = matrix[i][j];
            }
        }
        printf("please enter elements of vector b\n");
        for(int k=0; k<n; k++){
            scanf("%lf", &matrix[k][n]);
            matrix_og[k][n] = matrix[k][n];
        }
        x = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x[i] = (double*)malloc(1 * sizeof(double));

        printf("please enter the expected solution\n");
        x_given = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x_given[i] = (double*)malloc(1 * sizeof(double));
        for(int i=0; i<n; i++){
                scanf("%lf", &x_given[i][0]);
        }

        printf("the matrix you entered is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("vector b you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", matrix[i][n]);
        }
        printf("the expected solution you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", x_given[i][0]);
        }
    }


    //stin periptosi tou tixaiou pinaka xrisimopoio x=(1,1,...,1)^T kai b=Ax
    else if(atoi(argv[1]) == 2){//otan theloume tixaio pinaka, gia eukoleia stis prakseis evala anotato orio to 100
        printf("a random matrix will be created\n");
        printf("please enter the dimension of the matrix: ");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n+1) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n+1) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][j] = rand()%100;
                matrix_og[i][j] = matrix[i][j];
            }
        } 
        x = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x[i] = (double*)malloc(1 * sizeof(double));
        x_given = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x_given[i] = (double*)malloc(1 * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<1; j++){
                x_given[i][j] = 1;
            }
        }
        
        for(int i=0; i<n; i++){
            matrix[i][n] = 0.0;
            matrix_og[i][n] = 0.0;
        }
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][n] += (double)matrix[i][j] * x_given[i][0];
                matrix_og[i][n] += (double)matrix_og[i][j] * x_given[i][0];
            }
        }
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("vector b you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", matrix[i][n]);
        }
        printf("the expected solution is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", x_given[i][0]);
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
            matrix[i] = (double*)malloc((n+1) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n+1) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                fscanf(fptr, "%lf", &matrix[i][j]);
                matrix_og[i][j] = matrix[i][j];
            }
        }
        x = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x[i] = (double*)malloc(1 * sizeof(double));
        x_given = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x_given[i] = (double*)malloc(1 * sizeof(double));
        for(int i=0; i<n; i++){
            fscanf(fptr, "%lf", &matrix[i][n]);
            matrix_og[i][n] = matrix[i][n];
        }
        for(int i=0; i<n; i++){
            fscanf(fptr, "%lf", &x_given[i][0]);
        }
        printf("the matrix that was found in the file is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("vector b you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", matrix[i][n]);
        }
        printf("the expected solution you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", x_given[i][0]);
        }
        fclose(fptr);
    }

    else if(atoi(argv[1]) == 4){
        n=10;
        printf("implement lu decomposition on hilbert matrix with n=%d\n", n);//endeiktiki efarmogoi 1.a.3
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
        x = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x[i] = (double*)malloc(1 * sizeof(double));
        x_given = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x_given[i] = (double*)malloc(1 * sizeof(double));
        int flag=1;
        for(int i=0; i<n; i++){
            for(int j=0; j<1; j++){
                if(flag==1){
                    x_given[i][j] = 1;
                    flag = 0;
                }
                else if(flag==0){
                    x_given[i][j] = -2;
                    flag=1;
                }
            }
        }
        for(int i=0; i<n; i++){
                matrix[i][n] =0.0;
                matrix_og[i][n] = 0.0;
        }
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][n] += (double)matrix[i][j] * x_given[i][0];
                matrix_og[i][n] += (double)matrix_og[i][j] * x_given[i][0];
            }
        }
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("vector b you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", matrix[i][n]);
        }
        printf("the expected solution you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", x_given[i][0]);
        }
    }
    else if(atoi(argv[1]) == 5){//endeiktiki efarmogi 1.a.4
        printf("creating matrix from 1.a.4\n choose dimension of matrix from 100, 500, 1000\n");
        scanf("%d", &n);
        matrix = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix[i] = (double*)malloc((n+1) * sizeof(double));
        matrix_og = (double**)malloc((n) * sizeof(double*));
        for (int i = 0; i < n; i++)
            matrix_og[i] = (double*)malloc((n+1) * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                if(i!=j){
                    matrix[i][j] = (double)1.0/(i+j+1.0);
                }
                else{
                    matrix[i][j] = n;
                }
                matrix_og[i][j] = matrix[i][j];
            }
        }
        x = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x[i] = (double*)malloc(1 * sizeof(double));
        x_given = (double**)malloc(n * sizeof(double*));
        for (int i = 0; i < n; i++)
            x_given[i] = (double*)malloc(1 * sizeof(double));
        for(int i=0; i<n; i++){
            for(int j=0; j<1; j++){
                    x_given[i][j] = 1;
            }
        }
        for(int i=0; i<n; i++){
                matrix[i][n] =0.0;
                matrix_og[i][n] =0.0;
        }
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                matrix[i][n] += (double)matrix[i][j] * x_given[i][0];
                matrix_og[i][n] += (double)matrix_og[i][j] * x_given[i][0];
            }
        }
        printf("the matrix that was created is: \n");
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
               printf("%f\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("vector b you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", matrix[i][n]);
        }
        printf("the expected solution you entered is:\n");
        for(int i=0; i<n; i++){
            printf("%lf\n", x_given[i][0]);
        }
    }
    

    printf("Calculating LU decomposition of matrix\n");
    clock_t t;
    t = clock();

    //lu decomposition 
    lu_pivot(matrix, n);

    //epilisi sistimatos
    solve_LUx_b(matrix, n, x);

    //ektiposi apotelesmatos
    print_solution(x, n);

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("the algorithm took %f seconds to execute\n", time_taken);

    //ipologismos sfalmaton
    printf("Calculating error\n");

    double error1, error2, tmp1, tmp2;
    double** tmp;
    tmp = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        tmp[i] = (double*)malloc(1 * sizeof(double));
    substract_vectors(x_given, x, n, tmp, 1);
    tmp1 = calc_vector_norm(tmp, n, 1);
    tmp2 = calc_vector_norm(x_given, n, 1);
    error1 = (double)tmp1/tmp2;
    printf("||dx||/||x|| = %e\n", error1);

    mult_matrix_vector(matrix_og, x, n, tmp);
    substract_vectors(matrix_og, tmp, n, tmp, 0);
    tmp1 = calc_vector_norm(tmp, n, 1);
    tmp2 = calc_vector_norm(matrix_og, n, 0);
    error2 = (double)tmp1/tmp2;
    printf("||r||/||b|| = %e\n", error2);

    

    //apodesmeusi xorou
    for(int i = 0; i < n; i++)
        free(x[i]);
    free(x);
    for(int i = 0; i < n; i++)
        free(tmp[i]);
    free(tmp);
    for(int i = 0; i <n; i++)
        free(matrix[i]);
    free(matrix);
    for(int i = 0; i <n; i++)
        free(matrix_og[i]);
    free(matrix_og);
    for(int i = 0; i < n; i++)
        free(x_given[i]);
    free(x_given);
    return 0;
}


void lu_pivot(double** matrix, int n){
    int p[n];
    int i, j, k, pivot_ind = 0, temp_ind;
    double pivot, *temp_row;
    double temp;
    for (j = 0; j < n; ++j) {
        pivot = 0;
        for (i = j; i < n; ++i){
            if (fabs(matrix[i][j]) > fabs(pivot)) {
                pivot = matrix[i][j];
                pivot_ind = i;
            }
        }

        //enallagi grammon
        temp_row = matrix[j];
        matrix[j] = matrix[pivot_ind];
        matrix[pivot_ind] = temp_row;

        temp_ind  = p[j];
        p[j] = p[pivot_ind];
        p[pivot_ind] = temp_ind;

        for (k = j+1; k < n; ++k) {
            temp=matrix[k][j]/=matrix[j][j];
            for(int q=j+1;q<n;q++){
                matrix[k][q] -= temp*matrix[j][q];
            }
        }
    }

    printf("the LU matrix that was calculated is: \n");
    for(int q=0;q<n;q++){
        for(int l=0;l<n;l++){
            printf("%lf\t ",matrix[q][l]);
        }
        printf("\n");
    }
}

//pros ta empros antikatastasi gia na lisoume to Lz=b opou z=Ux
//meta pros ta piso antikatastasi gia Ux=z
void solve_LUx_b(double** matrix, int n, double** x){
    int i, j;
    double sum;
    double z[n];
    for(i=0;i<n;i++){       //pros ta empros
        sum=0.0;
        if(i==0){
            z[i]=matrix[i][n];
        }
        else{
            for(j=0;j<i;j++){
                sum+=matrix[i][j]*z[j];
            }
            z[i]=(matrix[i][n]-sum);
        }
    }
    for(i=n-1;0<=i;i--){ //pros ta piso
        sum=0.0;
        if(i == n-1){
            x[i][0]=z[i]/matrix[i][i];
        }
        else{
            for(j=n-1;i<=j;j--)
                sum+=matrix[i][j]*x[j][0];
            x[i][0]=(z[i]-sum)/matrix[i][i];
        }
    }
}

//ektiponei tin lisi tou sistimatos Ax=b
void print_solution(double** x, int n){
    printf("the calculated solution to the given linear system is:\n");
    for(int i=0;i<n;i++)
       printf("%lf\n", x[i][0]);
}

//ipologizei tin norma enos pinaka
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
    printf("the norm of given matrix is: %lf\n", max);
    return max;
}

//ipologizei tin norma dianismatos
double calc_vector_norm(double** x, int n, int flag){//flag=1-> norm tou x, flag=0, norm tou b(to dianisma einai stili enos pinaka)
    double maximum;
    double temp[n];
    for(int i = 0; i <n; i++) {
        if(flag==1){
            temp[i] = fabs(x[i][0]);
        }
        else if(flag==0){
            temp[i] = fabs(x[i][n]);
        }
    }
    maximum = temp[0];
    for(int i = 1; i < n; i++){
         maximum = fmax(maximum, temp[i]);
    }
    return maximum;
}

//afairesi pinakon
void substract_matrices(double** first, double** second, int n, double** diff){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            diff[i][j] = (double) first[i][j] - second[i][j];
        }
    }
}

//afairesi dianismaton
void substract_vectors(double** first, double** second, int n, double** diff, int flag){//flag=1->x, flag=0->b(to dianisma einai stili enos pinaka)
    for(int i=0; i<n; i++){
        if(flag == 1)
            diff[i][0] = (double) first[i][0] - second[i][0];
        else if(flag == 0)
            diff[i][0] = (double) first[i][n] - second[i][0];
    }
}

//ipologizei to dianisma Ax
void mult_matrix_vector(double** matrix, double** x, int n, double** temp){
    for(int i=0; i<n; i++){
        for(int j=0; j<1; j++){
            temp[i][j] = 0.0;
            for(int k=0; k<n; k++){
                temp[i][j] +=matrix[i][k] * x[k][j];
            }
        }
    }
}


