function [X] = rtrsyst(A,B,C, mb, nb)

M = length(A);
N = length(B);

if M <= mb && N <= nb
    
    X = sylvester(A,-B,C);
    
else
    
    if N <= M/2
        
        newM = newdims(A);
        
        if any(any(A(newM+1:end, 1:newM) ~= 0) == 1)
            fprintf("error")
        end
        
        X2 = rtrsyst(A(newM+1:end, newM+1:end), B, C(newM+1:end, :), mb, nb);
        C1 = gemm(- A(1:newM, newM +1:end), X2, C(1:newM, :));
        X1 = rtrsyst(A(1:newM, 1:newM), B, C1, mb, nb);
        X = [X1;X2];
        
    elseif M <= N/2
        
        newN = newdims(B);
        
        if any(any(B(newN+1:end, 1:newN) ~= 0) == 1)
            fprintf("error")
        end
        
        X1 = rtrsyst(A, B(1:newN, 1:newN), C(:, 1:newN), mb, nb);
        C2 = gemm(X1, B(1:newN, newN+1:end), C(:, newN+1:end));
        X2 = rtrsyst(A, B(newN+1:end, newN+1:end), C2, mb, nb);
        X = [X1, X2];
        
    else
        newM = newdims(A); newN = newdims(B);
        
        if any(any(A(newM+1:end, 1:newM) ~= 0) == 1)
            fprintf("error")
        end
        
        if any(any(B(newN+1:end, 1:newN) ~= 0) == 1)
            fprintf("error")
        end
        
        X21 = rtrsyst(A(newM+1:end, newM+1:end), B(1:newN, 1:newN), C(newM+1:end, 1:newN), mb, nb);
        C22 = gemm(X21, B(1:newN, newN+1:end), C(newM+1:end, newN+1:end));
        C11 = gemm(-A(1:newM, newM+1:end), X21, C(1:newM, 1:newN));
        X22 = rtrsyst(A(newM+1:end, newM+1:end), B(newN+1:end, newN+1:end), C22, mb, nb);
        X11 = rtrsyst(A(1:newM, 1:newM), B(1:newN, 1:newN), C11, mb, nb);
        C12 = gemm(-A(1:newM, newM+1:end), X22, C(1:newM, newN+1:end));
        C12 = gemm(X11, B(1:newN, newN+1:end), C12);
        X12 = rtrsyst(A(1:newM, 1:newM), B(newN+1:end, newN+1:end), C12, mb, nb);
        X = [X11, X12; X21, X22];     
    end
    
    
end