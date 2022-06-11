clear
clc
close all;

N = [4 8];          % No. of Receivers = No. of Transmitters; we'll consider 4 of each and 8 of each
SNR = -20:1:20;     % SNR will range from -20 to 20 dB
iter = 1000;        % Set the number of Monte-Carlo iterations

% QPSK
[avgBER_ZF, avgBER_VBLAST, avgBER_ML] = deal(zeros(size(SNR)));   %Initialize
qpskmod = comm.QPSKModulator('BitInput',true);
qpskdemod = comm.QPSKDemodulator('BitOutput',true);
qpskquads = [0.707106781186548 + 0.707106781186548i 0.707106781186548 - 0.707106781186548i -0.707106781186548 + 0.707106781186548i -0.707106781186548 - 0.707106781186548i];

figure
tic()
for n = 1:length(N)                                             % Loop through all N values

        % Setup for ML, will be used later
        all = dec2bin(0:(2^(N(n)*2)-1))' - '0';
        qpskall = zeros(size(all)/2);
        for a = 1:length(all)
            foo = qpskmod(all(:,a));
            qpskall(:,a) = foo;
        end
        disp('Finished generating symbols!');

    for s = 1:length(SNR)                                       % Loop through all SNR values
        SNR_lin = 10^(SNR(s)/10);                               % Convert SNR to linear scale
        H = 1/sqrt(2)*(rand(N(n), N(n)) +1i*rand(N(n),N(n)));   % Generate channel H
        [iterBER_ZF, iterBER_VBLAST, iterBER_ML] = deal(zeros(size(1, iter)));
          
        for m = 1:iter
            % This is a crude "progress bar", prints status every 1000 iterations
            if mod(m,1000) == 0
                disp('Currently at #'+string(m)+' at SNR '+ string(SNR(s))+ " for QPSKs N = "+string(N(n)));
            end

            b_in = randi([0 1],N(n)*2,1);                                           % Generate random signal
            x = qpskmod(b_in);                                                      % Shape signal to QPSK symbol
            noise = sqrt(1/(2*SNR_lin)) * (randn(N(n),1) + 1j*randn(N(n),1))/N(n);  % Generate complex noise samples with noise power = 1/SNR_lin
            y = H*x+noise;                                                          % Received signal

            % Zero-Forcing (ZF)
            G_ZF = inv(H'*H)*H';
            x_ZF = G_ZF*y;
            bOut_ZF = qpskdemod(x_ZF);
            iterBER_ZF(m) = sum(bOut_ZF ~= b_in)/length(b_in);

            % V-BLAST
            G = inv(H'*H)*H';                       % Calculate pseudoinverse
            [Gj, k] = deal(zeros([1, length(G)]));
            for j = 1:length(G)
                row = G(j,:,1);
                Gj(j) = norm(row);                  % find ||(G1)j||^2
            end
            [~,k(1)] = min(Gj(:));                  % find the location of the minimum
            r = zeros(length(y), length(y)+1);
            r(:,1) = y;
            [w] = deal(zeros([length(G(:,:,1)), length(G(:,:,1))]));
            [alpha, wye] = deal(zeros(1, length(G)));
            dist = zeros(1, length(qpskquads));
            for K = 1:length(H)
                w(k(K),:) = G(k(K),:,K);                            %w_ki = (G_i)_ki
                bar = reshape(w(k(K),:), [1, length(w(k(K),:))]);
                wye(k(K)) = bar*r(:,K);                             %y_ki = w_kiˆT*r_i
                for q = 1:length(qpskquads)
                    dist(q) = abs(norm(wye(k(K))-qpskquads(q)));
                end
                [~,loc] = min(dist);
                alpha(k(K)) = qpskquads(loc);                   %aˆˆ_ki = Q(y_ki)
                r(:,K+1) = r(:,K) - alpha(k(K))*G(:,k(K),K);    %r_{i+1} = r_i - aˆˆ_ki(H)_ki
                G(:,:,K+1) = G(:,:,K);
                G(:,k(K),K+1) = 0;                              %G_{i+1} = Gˆ+_k_i-
                Gj = zeros(1, length(G));
                for j = 1:length(G(:,:,K+1))
                    if ~ismember(j, k)
                        row = G(j,:,K+1);
                        Gj(j) = norm(row);                      % find ||(G1)j||^2
                    else
                        Gj(j) = realmax;
                    end
                end
                [~,k(K+1)] = min(Gj(:));                        % find the minimum ||(G1)j||^2 for all j
            end
            temp = reshape(alpha, [length(alpha) 1]);
            bOut_VBLAST = qpskdemod(temp);
            iterBER_VBLAST(m) = sum(bOut_VBLAST ~= b_in)/length(b_in);


            % Maximum Likelihood (ML)
            lengths = zeros(1, length(qpskall));
            for i = 1:length(qpskall)
                lengths(i) = norm(y-H*qpskall(:,i))^2;
            end
            [~,xhat] = min(lengths);
            bOut_ML = qpskdemod(qpskall(:,xhat));
            iterBER_ML(m) = sum(bOut_ML ~= b_in)/length(b_in);

        end
        avgBER_ZF(s) = mean(iterBER_ZF);
        avgBER_VBLAST(s) = mean(iterBER_VBLAST);
        avgBER_ML(s) = mean(iterBER_ML);
    end
    semilogy(SNR, avgBER_ZF, 'x-', SNR, avgBER_ML, '-o', SNR, avgBER_VBLAST, '-s')
    hold on
end
toc()
xlabel('SNR')
ylabel('BER')
legend('ZF, N = 4', 'ML, N = 4', 'VBLAST, N = 4', 'ZF, N = 8', 'ML, N  = 8', 'VBLAST, N = 8');
title('Decoding using QPSK symbols ('+string(iter)+' iterations)');



%--------------------------------------------------------------------------
% BPSK
clearvars -except N SNR iter

[avgBER_ZF, avgBER_VBLAST, avgBER_ML] = deal(zeros(size(SNR)));
bpskmod = comm.BPSKModulator;
bpskdemod = comm.BPSKDemodulator;
bpskquads = [1+0i -1+0i];   

figure
tic()
for n = 1:length(N)                                             % Loop through all N values
    % Setup for ML, will be used later
        all = dec2bin(0:(2^(N(n)))-1)' - '0';
        bpskall = zeros(size(all));
        for a = 1:length(all)
            foo = bpskmod(all(:,a));
            bpskall(:,a) = foo;
        end
    for s = 1:length(SNR)                                       % Loop through all SNR values
        SNR_lin = 10^(SNR(s)/10);                               % Convert SNR to linear scale
        H = 1/sqrt(2)*(rand(N(n), N(n)) +1i*rand(N(n),N(n)));   % Generate channel H
        [iterBER_ZF, iterBER_VBLAST, iterBER_ML] = deal(zeros(size(1, iter)));

        for m = 1:iter
%             if mod(m,1000) == 0
%                 disp('Currently at #'+string(m)+' at SNR '+ string(SNR(s))+ " for BPSKs N = "+string(N(n)));
%             end

            b_in = randi([0 1],N(n),1);                                             % Generate random signal
            x = bpskmod(b_in);                                                      % Shape signal to QPSK symbol
            noise = sqrt(1/(2*SNR_lin)) * (randn(N(n),1) + 1j*randn(N(n),1))/N(n);  % Generate complex noise samples with noise power = 1/SNR_lin
            y = H*x+noise;                                                          % Received signal

            % Zero-Forcing (ZF)
            G_ZF = inv(H'*H)*H';
            x_ZF = G_ZF*y;
            bOut_ZF = bpskdemod(x_ZF);
            iterBER_ZF(m) = sum(bOut_ZF ~= b_in)/length(b_in);

            % V-BLAST
            G = inv(H'*H)*H';
            [Gj, k] = deal(zeros([1, length(G)]));
            for j = 1:length(G)
                row = G(j,:,1);
                Gj(j) = norm(row);  % find ||(G1)j||^2
            end
            [~,k(1)] = min(Gj(:));  % find the minimum ||(G1)j||^2 for all j
            r = zeros(length(y), length(y)+1);
            r(:,1) = y;
            [w] = deal(zeros([length(G(:,:,1)), length(G(:,:,1))]));
            [alpha, wye] = deal(zeros(1, length(G)));
            dist = zeros(1, length(bpskquads));
            for K = 1:length(H)
                w(k(K),:) = G(k(K),:,K);
                bar = reshape(w(k(K),:), [1, length(w(k(K),:))]);
                wye(k(K)) = bar*r(:,K);
                for q = 1:length(bpskquads)
                    dist(q) = abs(norm(wye(k(K))-bpskquads(q)));
                end
                [~,loc] = min(dist);
                alpha(k(K)) = bpskquads(loc);
                r(:,K+1) = r(:,K) - alpha(k(K))*H(:,k(K));
                G(:,:,K+1) = G(:,:,K);
                G(:,k(K),K+1) = 0;
                Gj = zeros(1, length(G));
                for j = 1:length(G(:,:,K+1))
                    if ~ismember(j, k)
                        row = G(j,:,K+1);
                        Gj(j) = norm(row); % find ||(G1)j||^2
                    else
                        Gj(j) = realmax;
                    end
                end
                [~,k(K+1)] = min(Gj(:));    % find the minimum ||(G1)j||^2 for all j
            end
            temp = reshape(alpha, [length(alpha) 1]);
            bOut_VBLAST = bpskdemod(temp);
            iterBER_VBLAST(m) = sum(bOut_VBLAST ~= b_in)/length(b_in);


            % Maximum Likelihood (ML)
            lengths = zeros(1, length(bpskall));
            for i = 1:length(bpskall)
                lengths(i) = norm(y-H*bpskall(:,i))^2;
            end
            [val,xhat] = min(lengths);
            bOut_ML = bpskdemod(bpskall(:,xhat));
            iterBER_ML(m) = sum(bOut_ML ~= b_in)/length(b_in);

        end
        avgBER_ZF(s) = mean(iterBER_ZF);
        avgBER_VBLAST(s) = mean(iterBER_VBLAST);
        avgBER_ML(s) = mean(iterBER_ML);
    end
    semilogy(SNR, avgBER_ZF, 'x-', SNR, avgBER_ML, '-o', SNR, avgBER_VBLAST, '-s')
    hold on
end
xlabel('SNR')
ylabel('BER')
legend('ZF, N = 4', 'ML, N = 4', 'VBLAST, N = 4', 'ZF, N = 8', 'ML, N  = 8', 'VBLAST, N = 8');
title('Decoding using BPSK symbols ('+string(iter)+' iterations)');
toc()