function [radius] = calculate_r_psd(A, n, t, w)
    I=eye(n);
    D = diag(diag(A));
    D1 = inv(D);
    CL = -tril(A,-1);
    L = D1*CL;
    CU = -tril(A,1);
    U = D1*CU;
    k = inv(I-w*L);
    m = inv(I-w*U);
    B = m*k*D1*A;
    G = I-t*B;
    x = eig(G);
    radius = max(abs(x));

end
