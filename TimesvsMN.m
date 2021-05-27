clear all

save = true;

%% Computations
rng(0);

sizes = [4,8,16,32,64,128,256,512,1024,2048];

times = zeros(length(sizes));

for i = 1:length(sizes)
    
    M = sizes(i);
    A = schur(rand(M));
    
    for j = 1:length(sizes)
        
        N = sizes(j);
        
        B = schur(rand(N));
        C = rand(M,N);
        
        t = tic();
        X = rtrsyst(A, B, C, 64, 64);
        times(i,j) = toc(t);
    end
end

%% Plotting of the timings fixing N
figure()
for j = 1:length(sizes)
    loglog(sizes, times(:,j), 'DisplayName', "N = " + sizes(j))
    hold on
end
loglog(sizes, 1e-8*sizes(end,end)*sizes, '--k', 'DisplayName', "O(M)", 'linewidth', 2);
title('Computation times for rtrsyst(M,N)', 'Interpreter', 'latex')
xlabel('M', 'Interpreter', 'latex')
ylabel('Time [s]', 'Interpreter', 'latex')
legend('location', 'best', 'Interpreter', 'latex')
xlim([4, 2048])
hold off

if save
    saveas(gcf, 'figures/ExecutionTimes_MN.fig')
end

%% Surface plot of the timings
figure()
surf(sizes, sizes, times)
hold on
surf(sizes, sizes, sizes'*sizes)
xlabel('N', 'Interpreter', 'latex')
ylabel('M', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Computation times for rtrsyst(M,N)', 'Interpreter', 'latex')
set(gca,'xscale', 'log', 'yscale', 'log', 'zscale','log')
legend('Times from rtrsyst','M*N', 'location', 'best', 'Interpreter', 'latex')

if save
    saveas(gcf, 'figures/SurfExecutionTimes_MN.fig')
end
