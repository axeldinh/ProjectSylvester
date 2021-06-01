clear all
format long

save = true;
load = true;

if ~load
    %% Computations
    rng('default');

    num_runs = 10;

    sizes = [4,8,16,32,64,128,256,512,1024,2048];

    times = zeros(length(sizes));


    for i = 1:length(sizes)

        M = sizes(i);
        A = schur(rand(M));

        for j = 1:length(sizes)

            N = sizes(j);

            B = schur(rand(N));
            C = rand(M,N);

            time = 0;
            for n = 1:num_runs
                t = tic();
                X = rtrsyct(A, B, C, 64, 64);
                time = time + toc(t);
            end

            times(i,j) = time / num_runs;
            fprintf("M=%d, N=%d: %.4fs\n", M, N, times(i,j));
        end
    end

    % saving it for later use
    clear A B C X
    save results/TimesvsMN.mat
end

load results/TimesvsMN.mat

%% Plotting of the timings fixing N
figure()
for j = 1:length(sizes)
    loglog(sizes, times(:,j), 'DisplayName', "N = " + sizes(j))
    hold on
end
loglog(sizes, 1e-8*sizes(end,end)*sizes, '--k', 'DisplayName', "O(MN)", 'linewidth', 2);
title('Computation times for rtrsyct', 'Interpreter', 'latex')
xlabel('M', 'Interpreter', 'latex')
ylabel('Time [s]', 'Interpreter', 'latex')
legend('location', 'best', 'Interpreter', 'latex')
xlim([4, 2048])
hold off

if save
    saveas(gcf, 'figures/ExecutionTimes_MN.fig')
end

%% Surface plot of the timings
% figure()
% surf(sizes, sizes, times)
% hold on
% surf(sizes, sizes, sizes'*sizes)
% xlabel('N', 'Interpreter', 'latex')
% ylabel('M', 'Interpreter', 'latex')
% zlabel('Time [s]', 'Interpreter', 'latex')
% title('Computation times for rtrsyct(M,N)', 'Interpreter', 'latex')
% set(gca,'xscale', 'log', 'yscale', 'log', 'zscale','log')
% legend('Times from rtrsyct','M*N', 'location', 'best', 'Interpreter', 'latex')
% 
% if save
%     saveas(gcf, 'figures/SurfExecutionTimes_MN.fig')
% end

%% Contour plot of the timings

subplot(1,2,1)

X = sizes' * ones(1, length(sizes));
Y = ones(length(sizes),1) * sizes;
contourf(X,Y,times);
title('Contour plot of the rtrsyct times [s] vs M and N', 'Interpreter', 'latex')
xlabel('M', 'Interpreter', 'latex')
ylabel('N', 'Interpreter', 'latex')
set(gca,'xscale', 'log', 'yscale', 'log', 'ColorScale','log')
colorbar

subplot(1,2,2)
contourf(X, Y, X.*Y)
set(gca,'xscale', 'log', 'yscale', 'log', 'ColorScale','log')
xlabel('M', 'Interpreter', 'latex')
ylabel('N', 'Interpreter', 'latex')
title('M*N', 'Interpreter', 'latex')
colorbar

if save
    saveas(gcf, 'figures/ContourExecutionTimes_MN.fig')
end