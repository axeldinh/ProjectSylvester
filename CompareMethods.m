clear all
format long

rng('default')

save = true;
load = true;
num_runs = 10;

if ~load
    %% Computations

    sizes = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];

    times_schur = zeros(3,length(sizes));
    times_recon = zeros(3,length(sizes));
    times_rtrsyct = zeros(3,length(sizes));
    times_sylvester = zeros(3,length(sizes));
    times_lyap = zeros(3,length(sizes));

    %% M fixed

    M = sizes(end);

    for i = 1:length(sizes)

        N = sizes(i);

        A = rand(M);
        B = rand(N);
        C = rand(M, N);

        t = tic();
        for j = 1:num_runs
            [UA, TA] = schur(A);
            [UB, TB] = schur(B);
        end
        time = toc(t);
        times_schur(1, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = rtrsyct(TA, TB, UA*C*UB', 64, 64);
        end
        time = toc(t);
        times_rtrsyct(1, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = UA*X*UB';
        end
        time = toc(t);
        times_recon(1, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = sylvester(A, -B, C);
        end
        time = toc(t);
        times_sylvester(1, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = lyap(A, -B, -C);
        end
        time = toc(t);
        times_lyap(1, i) = time / num_runs;

        fprintf("M=%d, N=%d:\n", M, N);
        fprintf("rsyct = %fs\n", times_schur(1,i) + times_rtrsyct(1,i) + times_recon(1,i));
        fprintf("\tWith Schur = %fs, rtrsyct = %fs, Reconstruction  = %fs\n", times_schur(1,i), times_rtrsyct(1,i), times_recon(1,i));
        fprintf("sylvester() = %fs\n", times_sylvester(1,i));
        fprintf("lyap() = %fs\n\n", times_lyap(1,i));

    end

    %% N fixed

    N = sizes(end);

    for i = 1:length(sizes)

        M = sizes(i);

        A = rand(M);
        B = rand(N);
        C = rand(M, N);

        t = tic();
        for j = 1:num_runs
            [UA, TA] = schur(A);
            [UB, TB] = schur(B);
        end
        time = toc(t);
        times_schur(2, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = rtrsyct(TA, TB, UA*C*UB', 64, 64);
        end
        time = toc(t);
        times_rtrsyct(2, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = UA*X*UB';
        end
        time = toc(t);
        times_recon(2, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = sylvester(A, -B, C);
        end
        time = toc(t);
        times_sylvester(2, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = lyap(A, -B, -C);
        end
        time = toc(t);
        times_lyap(2, i) = time / num_runs;

        fprintf("M=%d, N=%d:\n", M, N);
        fprintf("rsyct = %fs\n", times_schur(2, i) + times_rtrsyct(2, i) + times_recon(2,i));
        fprintf("\tWith Schur = %fs, rtrsyct = %fs, Reconstruction  = %fs\n", times_schur(2,i), times_rtrsyct(2,i), times_recon(2,i));
        fprintf("sylvester() = %fs\n", times_sylvester(2,i));
        fprintf("lyap() = %fs\n\n", times_lyap(2,i));

    end

    %% M = N

    for i = 1:length(sizes)

        M = sizes(i);
        N = M;

        A = rand(M);
        B = rand(N);
        C = rand(M, N);

        t = tic();
        for j = 1:num_runs
            [UA, TA] = schur(A);
            [UB, TB] = schur(B);
        end
        time = toc(t);
        times_schur(3, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = rtrsyct(TA, TB, UA*C*UB', 64, 64);
        end
        time = toc(t);
        times_rtrsyct(3, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = UA*X*UB';
        end
        time = toc(t);
        times_recon(3, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = sylvester(A, -B, C);
        end
        time = toc(t);
        times_sylvester(3, i) = time / num_runs;

        t = tic();
        for j = 1:num_runs
            X = lyap(A, -B, -C);
        end
        time = toc(t);
        times_lyap(3, i) = time / num_runs;

        fprintf("M=N=%d:\n", N);
        fprintf("rsyct = %fs\n", times_schur(3, i) + times_rtrsyct(3, i) + times_recon(3, i));
        fprintf("\tWith Schur = %fs, rtrsyct = %fs, Reconstruction  = %fs\n", times_schur(3, i), times_rtrsyct(3,i), times_recon(3, i));
        fprintf("sylvester() = %fs\n", times_sylvester(3, i));
        fprintf("lyap() = %fs\n\n", times_lyap(3, i));

    end

    % save it for later use
    clear A B C TA TB UA UB X
    save results/CompareMethods.mat
    
end

load results/CompareMethods.mat

%% Compute the fraction of time for rtrsyct
y1 = [times_schur(1,:)', times_rtrsyct(1,:)', times_recon(1,:)'];
fracs1 = zeros(length(y1(:,1)), length(y1(1,:)));

for i = 1:length(y1(:,1))
    for j = 1:length(y1(1,:))
        fracs1(i,j) = y1(i,j) / sum(y1(i,:));
    end
end

y2 = [times_schur(2,:)', times_rtrsyct(2,:)', times_recon(2,:)'];
fracs2 = zeros(length(y2(:,1)), length(y2(1,:)));

for i = 1:length(y2(:,1))
    for j = 1:length(y2(1,:))
        fracs2(i,j) = y2(i,j) / sum(y2(i,:));
    end
end

y3 = [times_schur(3,:)', times_rtrsyct(3,:)', times_recon(3,:)'];
fracs3 = zeros(length(y3(:,1)), length(y3(1,:)));

for i = 1:length(y3(:,1))
        for j = 1:length(y3(1,:))
            fracs3(i,j) = y3(i,j) / sum(y3(i,:));
        end
    end

%% Compare methods

for i = 1:3
    figure()

    loglog(sizes, times_rtrsyct(i,:) + times_schur(i,:) + times_recon(i,:), sizes, times_sylvester(i,:), sizes, times_lyap(i,:))
    if i == 1
        title("Computation Times (M="+num2str(sizes(end))+")", 'Interpreter', 'latex')
        xlabel("N", 'Interpreter', 'latex')
        save_name = "CompareMethods_Mfixed.fig";
    elseif i == 2
        title("Computation Times (N="+num2str(sizes(end))+")", 'Interpreter', 'latex')
        xlabel("M", 'Interpreter', 'latex')
        save_name = "CompareMethods_Nfixed.fig";
    elseif i == 3
        title("Computation Times (M=N)", 'Interpreter', 'latex')
        xlabel("M=N", 'Interpreter', 'latex')
        save_name = "CompareMethods_MequalsN.fig";
    end
    ylabel("Time [s]", 'Interpreter', 'latex')
    legend("rsyct", "sylvester", "lyap", 'Interpreter', 'latex', 'location', 'best')
    xlim([sizes(1), sizes(end)])

    if save
        saveas(gcf, "figures/" + save_name);
    end
end

%% Compare profiles for schur, rtrsyct, recon

for i = 1:3
    
    figure()
    if i == 1
        bar(log2(sizes), fracs1, 1.0, 'stacked')
        title("Computations Times for Calls needed for rsyct (M="+num2str(sizes(end))+")", 'Interpreter', 'latex')
        xlabel("N", 'Interpreter', 'latex')
        save_name = 'figures/CallsTimes_Mfixed.fig';
    elseif i == 2
        bar(log2(sizes), fracs2, 1.0, 'stacked')
        title("Computations Times for Calls needed for rsyct (N="+num2str(sizes(end))+")", 'Interpreter', 'latex')
        xlabel("M", 'Interpreter', 'latex')
        save_name = 'figures/CallsTimes_Nfixed.fig';
    elseif i == 3
        bar(log2(sizes), fracs3, 1.0, 'stacked')
        title("Computations Times for Calls needed for rsyct (M=N)", 'Interpreter', 'latex')
        xlabel("M=N", 'Interpreter', 'latex')
        save_name = 'figures/CallsTimes_MequalsN.fig';
    end
    set(gca, 'XTickLabel', {'2^2', '2^3', '2^4', '2^5', '2^6', '2^7', '2^8',...
        '2^9', '2^{10}', '2^{11}'})
    ylabel("Fraction of Computations", 'Interpreter', 'latex')
    legend("Schur Decomposition", "rtrsyct", "Reconstruction", 'Interpreter', 'latex', 'location', 'best')

    if save
        saveas(gcf, save_name);
    end
end

