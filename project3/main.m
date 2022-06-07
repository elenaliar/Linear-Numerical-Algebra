clc; clear; format long;
A=[];
rng(1);

fprintf('> These are all the available options and their code:\n');
fprintf('\t1) Give elements of matrix row by row(Code: 1)\n');
fprintf('\t2) creating matrix from lecture(Code: 2)\n');
fprintf('\t3) Read a matrix from a file (Code: 3)\n');
fprintf('\t4) Create a random matrix (Code: 4)\n');
code = input('\n> Give the preferable Code: ');
if(code == 1)
    fprintf("please enter data of the matrix\n");
    m = input('\n> enter dimension m: ');
    n = input('\n> enter dimension n: ');
    for i=1:m
        for j=1:n
            A(i,j) = input('enter element');
        end
    end
    fprintf("the matrix you entered is:\n");
    %disp(A);

elseif (code==2)
    m = 3;
    n=3;
    A(1,1) = 0;
    A(1,2) = 1;
    A(1,3) = 1;
    A(2, 1) = 1;
    A(2, 2) = 2;
    A(2, 3) = 3;
    A(3, 1) = 1;
    A(3, 2) = 2;
    A(3, 3) = 3;
    fprintf("the matrix A is:\n");
    disp(A);

elseif(code == 3)
    fileID = fopen('sample3.txt','r');%edo vazoume to onoma tou arxeiou pou theloume
    formatSpec = '%d';
    B= [];
    B = fscanf(fileID,formatSpec);
    %disp(B);
    m = B(1,1);
    n = B(2,1);
    count=3;
    for i=1:m
        for j=1:n
            A(i,j) = B(count, 1);
            count = count+1;

        end
    end
    fprintf("the matrix you entered is:\n");
    %disp(A);

elseif(code==4)
    m = input('\n> enter dimension m: ');
    n = input('\n> enter dimension n: ');
    for i=1:m
        for j=1:n
            A(i,j) = randi([0 100]);
        end
    end
    fprintf("the matrix that was created is:\n");
    %disp(A);
end

fprintf("Calculating Householder QR decomposition\n");

tic
[Q, R] = find_qr(A, m, n);
toc
fprintf("matrix Q is:\n");
%disp(Q);
fprintf("matrix R is:\n");
%disp(R);
fprintf("Calculating erros\n");
QT = transpose(Q);
QR = Q*R;
QTQ = QT*Q;
R1 = inv(R);
AR1 = A*R1;
condition = cond(R);
fprintf("condR = ");
disp(condition);
I = eye(n,n);
er1 = norm(A-QR, "inf");
er2 = norm(QTQ-I, "inf");
er3 = norm(AR1-Q, "inf");
fprintf("||A-QR|| = ");
disp(er1);
fprintf("||QTQ-I|| = ");
disp(er2);
fprintf("||AR1-Q|| = ");
disp(er3);


