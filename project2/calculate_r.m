function [radius]=calculate_r(A, n, t, w)
    I=eye(n);
    D = diag(diag(A));
    D1 = inv(D);
    CL = -tril(A,-1);
    L = D1*CL;
    G = I-t*(inv(I-w*L))*D1*A;
    x = eig(G);
    radius = max(abs(x));
end