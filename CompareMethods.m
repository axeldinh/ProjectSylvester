clear all

format long

save = true;

sizes = 1:150;

steps = 1;

times_schur = zeros(1,length(sizes));
times_recon = zeros(1,length(sizes));
times_rtrsyst = zeros(1,length(sizes));
times_sylvester = zeros(1,length(sizes));
times_lyap = zeros(1,length(sizes));

    
for i = 1:length(sizes)
    
    N = sizes(i);

    A = rand(N);
    B = rand(N);
    C = rand(N, N);

    t = tic();
    for j = 1:steps
        [UA, TA] = schur(A);
        [UB, TB] = schur(B);
    end
    time = toc(t);
    times_schur(i) = time;

    t = tic();
    for j = 1:steps
        X = rtrsyst(TA, TB, UA*C*UB', 512, 512);
    end
    time = toc(t);
    times_rtrsyst(i) = time;

    t = tic();
    for j = 1:steps
        X = UA*X*UB';
    end
    time = toc(t);
    times_recon(i) = time;

    t = tic();
    for j = 1:steps
        X = sylvester(A, -B, C);
    end
    time = toc(t);
    times_sylvester(i) = time;

    t = tic();
    for j = 1:steps
        X = lyap(A, -B, -C);
    end
    time = toc(t);
    times_lyap(i) = time;
end

% Compute the fraction of time for rtrsyst
y = [times_schur', times_rtrsyst', times_recon'];
fracs = zeros(length(y(:,1)), length(y(1,:)));

for i = 1:length(y(:,1))
    for j = 1:length(y(1,:))
        fracs(i,j) = y(i,j) / sum(y(i,:));
    end
end

% Compare methods
figure()

loglog(sizes, times_rtrsyst + times_schur + times_recon, sizes, times_sylvester, sizes, times_lyap)
title("Computation Times", 'Interpreter', 'latex')
xlabel("M=N", 'Interpreter', 'latex')
ylabel("Time [s]", 'Interpreter', 'latex')
legend("rsyst", "sylvester", "lyap", 'Interpreter', 'latex', 'location', 'best')

if save
    saveas(gcf, 'CompareMethods.fig');
end

% Compare profiles for schur, rtrsyst, recon
figure()
bar(sizes, fracs, 1.0, 'stacked')
title("Computations Times for Calls needed for rsyst", 'Interpreter', 'latex')
xlabel("M=N", 'Interpreter', 'latex')
ylabel("Time [s]", 'Interpreter', 'latex')
legend("Schur Decomposition", "rtrsyst", "Reconstruction", 'Interpreter', 'latex', 'location', 'best')

if save
    saveas(gcf, 'figures/CallsTimes.fig');
end
