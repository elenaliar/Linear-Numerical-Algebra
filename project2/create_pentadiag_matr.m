function [A] = create_pentadiag_matr(n, a, b, c, d, e)%dimiourgei pentadiagonio pinaka
    for i=1:n
        A(i,i) = e;
    end
    for i=2:n
        A(i,i-1) = -b;
    end
    for i=1:n-1
        A(i, i+1) = -c;
    end
    for i=1:n-2
        A(i,i+2) = -d;
    end
    for i=3:n
        A(i,i-2) = -a;
    end
end