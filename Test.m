format long

sizes = [8, 16, 32, 64, 128, 256, 512, 1024];
errors = zeros(length(sizes));
times_schur = zeros(length(sizes));
times_recon = zeros(length(sizes));
times_mysylv = zeros(length(sizes));
times_matsylv = zeros(length(sizes));

for i = 1:length(sizes)
    
    for j = 1:length(sizes)
        
        M = sizes(i);
        N = sizes(j);
        
        A = rand(M);
        B = rand(N);
        C = rand(M, N);
        
        t = tic();
        [UA, TA] = schur(A);
        [UB, TB] = schur(B);
        time = toc(t);
        times_schur(i,j) = time;
        
        t = tic();
        X1 = rtrsyst(TA, TB, UA*C*UB', 64, 64);
        time = toc(t);
        times_mysylv(i, j) = time;
        
        t = tic();
        X1 = UA*X1*UB';
        time = toc(t);
        times_recon(i,j) = time;
        
        t = tic();
        X2 = sylvester(A, -B, C);
        time = toc(t);
        times_matsylv(i, j) = time;
        
        errors(i,j) = norm(X1 - X2);
        
    end
end

for i = 1:length(sizes)
    
    figure()
    title("M = " + num2str(sizes(i)))
    loglog(sizes, times_mysylv(i, :) + times_schur(i,:) + times_recon(i,:), sizes, times_matsylv(i, :))
    xlabel("N")
    ylabel("times")
    legend("My Code", "Matlab")
end

% surf(sizes, sizes, times_mysylv + times_schur + times_recon)
% hold on
% surf(sizes, sizes, times_matsylv)

xscale('log')