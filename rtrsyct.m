function [C] = rtrsyct(A,B,C, mb, nb)

M = length(A);
N = length(B);

if M <= mb && N <= nb
    
    X = sylvester(A,-B,C);
    
else
    
    if N <= M/2
        
        newM = newdims(A);
        
        C(newM+1:end, :) = rtrsyct(A(newM+1:end, newM+1:end), B, C(newM+1:end, :), mb, nb);
        C(1:newM, :) = gemm(- A(1:newM, newM +1:end), C(newM+1:end, :), C(1:newM, :));
        C(1:newM, :) = rtrsyct(A(1:newM, 1:newM), B, C(1:newM, :), mb, nb);
        
    elseif M <= N/2
        
        newN = newdims(B);
        
        C(:, 1:newN) = rtrsyct(A, B(1:newN, 1:newN), C(:, 1:newN), mb, nb);
        C(:, newN+1:end) = gemm(C(:, 1:newN), B(1:newN, newN+1:end), C(:, newN+1:end));
        C(:, newN+1:end) = rtrsyct(A, B(newN+1:end, newN+1:end), C(:, newN+1:end), mb, nb);
        
    else
        newM = newdims(A); newN = newdims(B);
        
        C(newM+1:end, 1:newN) = rtrsyct(A(newM+1:end, newM+1:end), B(1:newN, 1:newN), C(newM+1:end, 1:newN), mb, nb);
        C(newM+1:end, newN+1:end) = gemm(C(newM+1:end, 1:newN), B(1:newN, newN+1:end), C(newM+1:end, newN+1:end)); % the 2 gemms are parallelizable
        C(1:newM, 1:newN) = gemm(-A(1:newM, newM+1:end), C(newM+1:end, 1:newN), C(1:newM, 1:newN));
        C(newM+1:end, newN+1:end) = rtrsyct(A(newM+1:end, newM+1:end), B(newN+1:end, newN+1:end), C(newM+1:end, newN+1:end), mb, nb);
        C(1:newM, 1:newN) = rtrsyct(A(1:newM, 1:newM), B(1:newN, 1:newN), C(1:newM, 1:newN), mb, nb);
        C(1:newM, newN+1:end) = gemm(-A(1:newM, newM+1:end), C(newM+1:end, newN+1:end), C(1:newM, newN+1:end));
        C(1:newM, newN+1:end) = gemm(C(1:newM, 1:newN), B(1:newN, newN+1:end), C(1:newM, newN+1:end));
        C(1:newM, newN+1:end) = rtrsyct(A(1:newM, 1:newM), B(newN+1:end, newN+1:end), C(1:newM, newN+1:end), mb, nb);
    end
    
    
end