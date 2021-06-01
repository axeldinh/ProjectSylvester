function [M] = newdims(A)

i = floor(length(A)/2);

if A(i+1, i) == 0
    M = i;
else
    M = i+1;
end

