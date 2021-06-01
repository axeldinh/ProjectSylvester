clear all
format long

save = true;
load = true;
num_runs = 10;

if ~load
    %% Computations
    rng('default');

    M = 2048;
    N = 2048;

    A = schur(rand(M));
    B = schur(rand(N));
    C = rand(M, N);

    sizes = 4:20:204;
    times = zeros(length(sizes),1);

    for i = 1:length(sizes)
        for j = 1:length(sizes)
            Mb = sizes(i);
            Nb = sizes(j);

            time = 0;
            for n = 1:num_runs
                t = tic();
                X = rtrsyct(A, B, C, Mb, Nb);
                time = time + toc(t);
            end
            times(i,j) = time / num_runs;
            fprintf("Mb = %d, Nb = %d: %.2fs\n", Mb, Nb, times(i,j));
        end
    end

    % Save it for later use
    clear A B C X
    save results/FindNbMb.mat
end

load results/FindNbMb.mat

%% Plotting of a surface plot of the times

% figure()
% surf(sizes, sizes, times)
% xlabel('Mb')
% ylabel('Nb')
% zlabel('Time')
% 
% if save
%     saveas(gcf, 'figures/MbNb.fig');
% end

%% Plotting of a contour plot of the times

figure()
X = sizes' * ones(1, length(sizes));
Y = ones(length(sizes), 1) * sizes;
contourf(X, Y, times);
title("Times [s] vs Mb and Nb (M=N="+num2str(M)+")"', 'Interpreter', 'latex')
xlabel('Mb', 'Interpreter', 'latex')
ylabel('Nb', 'Interpreter', 'latex')
colorbar

if save
    saveas(gcf, 'figures/contour_MbNb.fig');
end

figure()
X = sizes(3:end)' * ones(1, length(sizes)-2);
Y = ones(length(sizes)-2, 1) * sizes(3:end);
contourf(X, Y, times(3:end, 3:end));
title("Times [s] vs Mb and Nb (M=N="+num2str(M)+")", 'Interpreter', 'latex')
xlabel('Mb', 'Interpreter', 'latex')
ylabel('Nb', 'Interpreter', 'latex')
colorbar

if save
    saveas(gcf, 'figures/contour_MbNbZoom.fig');
end

%% Getting the run with the best time (might not be the wisest pick)

[~, x_id] = min(times);
[~, y_id] = min(min(times));
x_id = x_id(y_id);

best_Mb = sizes(x_id);
best_Nb = sizes(y_id);
best_time = times(x_id, y_id);

fprintf("Best time: " + num2str(best_time))
fprintf("\nAttained with: Mb=" + num2str(best_Mb) + ", Nb=" + num2str(best_Nb))
fprintf("\nTotal components in Sylvester() call: " + num2str(best_Mb*best_Nb) + "\n")
