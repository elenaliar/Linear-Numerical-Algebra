clc; clear; format long;
maxiter = 100;
A=[];

fprintf('> These are all the available options and their code:\n');
fprintf('\t1) Create a Pentadiagonal matrix with parameters a, b, c, d, e (Code: 1)\n');
fprintf('\t2) Specific Applications (Code: 2)\n');
fprintf('\t3) Read a matrix from a file (Code: 3)\n');
fprintf('\t4) Create a random strict diagonally dominant matrix (Code: 4)\n');
code = input('\n> Give the preferable Code: ');
if(code == 1)
    fprintf("please enter data of the matrix\n");
    n = input('\n> enter dimension n: ');
    a = input('\n> enter a: ');
    b = input('\n> enter b: ');
    c = input('\n> enter c: ');
    d = input('\n> enter d: ');
    e = input('\n> enter e: ');
    [A] = create_pentadiag_matr(n, a, b, c, d, e);
    fprintf("the matrix you entered is:\n");
    %disp(A);

elseif (code==2)
    fprintf("there are 3 available options:\n a=0.1, b=0.2, c=0.3, d=0.4\n a=0.4, b=0.3, c=0.2, d=0.1\n a=1.2, b=0.9, c=0.6, d=0.3\n");
    option = input('enter 1 for first option, 2 for the second and 3 for the third: ');
    n = input('\n> enter dimension n: ');
    if(option == 1)
        [A] = create_pentadiag_matr(n, 0.1, 0.2, 0.3, 0.4, 4);
    elseif (option==2)
        [A] = create_pentadiag_matr(n, 0.4, 0.3, 0.2, 0.1, 4);
    elseif(option == 3)
        [A] = create_pentadiag_matr(n, 1.2, 0.9, 0.6, 0.3, 4);
    end
    fprintf("the matrix you entered is:\n");
    %disp(A);

elseif(code == 3)
    fileID = fopen('case1i.txt','r');%edo vazoume to onoma tou arxeiou pou theloume
    formatSpec = '%f';
    B= [];
    B = fscanf(fileID,formatSpec);
    disp(B);
    n = B(1,1);
    a=B(2,1);
    b=B(3,1);
    c=B(4,1);
    d=B(5,1);
    e=B(6,1);
    [A] = create_pentadiag_matr(n, a, b, c, d, e);
    fprintf("the matrix you entered is:\n");
    %disp(A);

elseif(code==4)
    n = input('\n> enter dimension n: ');
    [A]=create_random_matrix(n);
    fprintf("the matrix that was created is:\n");
    %disp(A);
end
b = A*ones(n,1); % lisi x = (1,1,1,....,1)^T
fprintf("vector b you entered is:\n");
%disp(b);
minradius = +Inf;
fprintf("Calculating solution of linear system Ax=b with ESOR method and using x0=b\n");
for t=0.1:0.1:1.9
    for w=0.1:0.1:1.9
        radius = calculate_r(A,n, t, w);
        if((radius<minradius) || ((t==0.1) && (w==0.1)))
            minradius = radius;
            topt = t;
            wopt = w;
            %fprintf("better parameters found\n");
        end
    end
end
starttime = tic;
[x, itcount] = esor_find_sol(A, n, topt, wopt, maxiter, b);
endtime = toc;
fprintf("Optimal value of t:");
disp(topt);
fprintf("Optimal value of w:");
disp(wopt);
fprintf("itcount:");
disp(itcount);
fprintf("radius of iterative matrix:");
disp(minradius);
if(itcount < maxiter)
    fprintf("approximate solution for the given linear system was found with the aforementioned optimal values for t and w\n");
    %disp(x);
else
    fprintf("approximate solution for the given linear system was not found\n");
end
fprintf("execution time of esor method for Ax=b resolution with optimal values of t and w from radius: ");
disp(endtime);


    















