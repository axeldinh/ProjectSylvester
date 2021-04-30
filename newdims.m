function [M] = newdims(A)

flag1 = true;
flag2 = true;
i = 1;
M = 0;
% Look for the maximum square matrix only with zeros in A (in bottom left)
% Then returns the dimensions of A11.
while flag1
    
    if any(any(A(end-i+1:end, 1:i) ~= 0) == 1)
        flag1 = false;
        j = i-1;
        
        while flag2
            if any(any(A(end-i+1:end, 1:j) ~= 0) == 1)
                flag2 = false;
                i = i-1;
            else
                i = i+1;
            end
        end
        
    else
        i = i+1;
    end
end

% if M == 1
%     M = length(A);
% else
%     M = length(A)- i;
% end

M = length(A) - i;