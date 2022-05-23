function [x, itcount] = psd_find_sol(A, n, t, w, maxiter, b)
    e =  0.5*(10^(-6));
    itcount=0;
    solution = b;
    Ux = zeros(n,1);
    while(itcount< maxiter)
        temp = zeros(n,1);
        temp2 = zeros(n,1);%to proto for loop einai idio me tin esor
               for i = 1:n % epanaliptiki methodos akolouthontas to sxima tou Niehthammer
                    if i > 2 && i <= n-2
                        Ux(i) = dot(-A(i,i+1:i+2)/A(i,i),solution(i+1:i+2)'); 
                        temp2(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp2(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + t*Ux(i) + (t*b(i)/A(i,i));                    
                    elseif i <=2
                        Ux(i) = dot(-A(i,i+1:i+2)/A(i,i),solution(i+1:i+2)'); 
                        temp2(i) = (1-t)*solution(i) + w*dot(-A(i,1:i-1)/A(i,i), temp2(1:i-1)') + (t-w)*dot(-A(i,1:i-1)/A(i,i),solution(1:i-1)') + t*Ux(i) + (t*b(i)/A(i,i));
                    elseif i == n-1
                        Ux(i) = dot(-A(i,n)/A(i,i),solution(n)); 
                        temp2(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp2(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + t*Ux(i) + (t*b(i)/A(i,i));
                    else
                        Ux(i) = 0; 
                        temp2(i) = (1-t)*solution(i) + w*dot(-A(i,i-2:i-1)/A(i,i), temp2(i-2:i-1)') + (t-w)*dot(-A(i,i-2:i-1)/A(i,i),solution(i-2:i-1)') + (t*b(i)/A(i,i));
                    end
                end
               
               for i = n:-1:1
                   if i > 2 && i <= n-2
                       temp(i) = temp2(i) + w*(dot(-A(i,i+1:i+2)/A(i,i),temp(i+1:i+2)') - Ux(i));
                    elseif i <=2
                        temp(i) = temp2(i) + w*(dot(-A(i,i+1:i+2)/A(i,i),temp(i+1:i+2)') - Ux(i));
                    elseif i == n-1
                        temp(i) = temp2(i) + w*(dot(-A(i,n)/A(i,i),temp(n)) - Ux(i));
                    elseif i == n
                        temp(i) = temp2(i);
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