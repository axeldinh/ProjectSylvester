clear all

format long

save = true;

M = 1024;
N = 1024;

A = schur(rand(M));
B = schur(rand(N));
C = rand(M, N);

sizes = 4:20:512;
times = zeros(length(sizes),1);

for i = 1:length(sizes)
    for j = 1:length(sizes)
        Mb = sizes(i);
        Nb = sizes(j);
        
        t = tic();
        X = rtrsyst(A, B, C, Mb, Nb);
        time = toc(t);
        times(i,j) = time;
    end
end

figure()
surf(sizes, sizes, times)
xlabel('Mb')
ylabel('Nb')
zlabel('Time')

if save
    saveas(gcf, 'figures/MbNb.fig');
end

[~, x_id] = min(times);
[~, y_id] = min(min(times));
x_id = x_id(y_id);

best_Mb = sizes(x_id);
best_Nb = sizes(y_id);
best_time = times(x_id, y_id);

fprintf("Best time: " + num2str(best_time))
fprintf("\nAttained with: Mb=" + num2str(best_Mb) + ", Nb=" + num2str(best_Nb))
fprintf("\nTotal components in Sylvester() call: " + num2str(best_Mb*best_Nb) + "\n")
