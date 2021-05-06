clc
close all
clear all

format long

sizes = [8, 16, 32, 64, 128, 256, 512, 1024];
errors = zeros(length(sizes));
times_my = zeros(length(sizes));
times_mat = zeros(length(sizes));

for i = 1:length(sizes)
    
    for j = 1:length(sizes)
        
        M = sizes(i);
        N = sizes(j);
        
        A = schur(rand(M));
        B = schur(rand(N));
        
        
        C = rand(M, N);
        
        t = tic();
        X1 = rtrsyst(A, B, C, 64, 64);
        time = toc(t);
        times_my(i, j) = time;
        
        t = tic();
        X2 = sylvester(A, -B, C);
        time = toc(t);
        times_mat(i, j) = time;
        
        errors(i,j) = norm(X1 - X2);
        
    end
end

for i = 1:length(sizes)
    
    figure()
    title("M = " + num2str(sizes(i)))
    loglog(sizes, times_my(i, :), sizes, times_mat(i, :))
    xlabel("N")
    ylabel("times")
    legend("My Code", "Matlab")
end