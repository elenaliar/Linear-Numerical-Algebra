function [x, itcount] = esor_find_sol(A, n, t, w, maxiter, b)
    e =  0.5*(10^(-6));
    itcount=0;
    solution = b;
    while(itcount< maxiter)
        temp = zeros(n,1);
               for i = 1:n % iterative method using components
                  if i > 2 && i <= n-2
                    temp(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + t*dot(-A(i,i+1:i+2)/A(i,i),solution(i+1:i+2)') + (t*b(i)/A(i,i));
                  elseif i <=2
                    temp(i) = (1-t)*solution(i) + w*dot(-A(i,1:i-1)/A(i,i), temp(1:i-1)') + (t-w)*dot(-A(i,1:i-1)/A(i,i),solution(1:i-1)') + t*dot(-A(i,i+1:i+2)/A(i,i),solution(i+1:i+2)') + (t*b(i)/A(i,i));
                  elseif i == n-1
                    temp(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + t*dot(-A(i,n)/A(i,i),solution(i+1)) + (t*b(i)/A(i,i));
                  else
                    temp(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + (t*b(i)/A(i,i));
                  end
               end
        itcount = itcount+1;
        if norm(temp-solution, +Inf) < e
            break
        else
            solution = temp;
        end
    end  
    x= solution;
end