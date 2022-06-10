clear
clc
close all;
N = [4 8];          % No. of Receivers = No. of Transmitters; we'll consider 4 of each and 8 of each
SNR = -20:1:20;     % SNR will range from -20 to 20 dB
iter = 1000;        % Set the number of Monte-Carlo iterations
[avgBER_ZF, avgBER_VBLAST, avgBER_ML] = deal(zeros(size(SNR)));   %Initialize
qpskmod = comm.QPSKModulator;
qpskdemod = comm.QPSKDemodulator;
qpskquads = [0.707106781186548 + 0.707106781186548i 0.707106781186548 - 0.707106781186548i -0.707106781186548 + 0.707106781186548i -0.707106781186548 - 0.707106781186548i];

figure
tic()
for n = 1:length(N)                                             % Loop through all N values
    for s = 1:length(SNR)                                       % Loop through all SNR values
        SNR_lin = 10^(SNR(s)/10);                               % Convert SNR to linear scale
        H = 1/sqrt(2)*(rand(N(n), N(n)) +1i*rand(N(n),N(n)));   % Generate channel H
        [iterBER_ZF, iterBER_VBLAST, iterBER_ML] = deal(zeros(size(1, iter)));

        % Setup for ML, will be used later
        all = dec2bin(0:(2^(N(n)))-1)' - '0';
        qpskall = zeros(size(all));
        for a = 1:length(all)
            foo = qpskmod(all(:,a));
            qpskall(:,a) = foo;
        end
          
        for m = 1:iter
%             % This is a crude "progress bar", prints status every 1000 iterations
%             if mod(i,1000) == 0
%                 disp('Currently at #'+string(i)+' at SNR '+ string(SNR(s)));
%             end

            b_in = randi([0 1],N(n),1);                                             % Generate random signal
            x = qpskmod(b_in);                                                      % Shape signal to QPSK symbol
            noise = sqrt(1/(2*SNR_lin)) * (randn(N(n),1) + 1j*randn(N(n),1))/N(n);  % Generate complex noise samples with noise power = 1/SNR_lin
            % Maybe divide by N here?
            y = H*x+noise;                                                          % Received signal

%             % Zero-Forcing (ZF)
%             G_ZF = inv(H'*H)*H';
%             x_ZF = G_ZF*y;
%             bOut_ZF = qpskdemod(x_ZF);
%             iterBER_ZF(m) = sum(bOut_ZF ~= b_in)/length(b_in);

            % V-BLAST
            G = inv(H'*H)*H';
            % A(:,:,x) is page x
            [Gj, k] = deal(zeros([1, length(G)]));
            for j = 1:length(G)
                row = G(j,:,1);
                Gj(j) = norm(row); % find ||(G1)j||^2
            end
            [~,k(1)] = min(Gj(:));    % find the minimum ||(G1)j||^2 for all j
            r = y;
            [w, alpha] = deal(zeros([length(G(:,:,1)), length(G(:,:,1))]));
            for K = 1:length(H)
                w(k(K),:) = G(k(K),:,K);
                bar = reshape(w(k(K),:), [1, length(w(k(K),:))]);
                wye(k(K)) = bar*r(:,K); % NEED TO FIX THIS
%                 foo = wye(k(K),:);
%                 alpha(k(K),:) = qpskdemod(reshape(foo, [length(foo), 1]));
                for q = 1:length(qpskquads)
                    dist(q) = abs(norm(wye(k(K))-qpskquads(q)));
                end
                [~,loc] = min(dist(q));
                alpha(k(K)) = qpskquads(loc);
                r(:,K+1) = r(:,K) - alpha(k(K),:)*H(:,k(K));
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
            bOut_VBLAST = qpskdemod(r(:,end));
            iterBER_VBLAST(m) = sum(bOut_VBLAST ~= b_in)/length(b_in);


%             % Maximum Likelihood (ML)
%             for i = 1:length(qpskall)
%                 lengths(i) = norm(y-H*qpskall(:,i))^2;
%             end
%             [val,xhat] = min(lengths);
%             bOut_ML = qpskdemod(qpskall(:,xhat));
%             iterBER_ML(m) = sum(bOut_ML ~= b_in)/length(b_in);

        end
        avgBER_ZF(s) = mean(iterBER_ZF);
        avgBER_VBLAST(s) = mean(iterBER_VBLAST);
        avgBER_ML(s) = mean(iterBER_ML);
    end
    %semilogy(SNR, avgBER_ZF)
    semilogy(SNR, avgBER_ZF, 'x-', SNR, avgBER_ML, '-o', SNR, avgBER_VBLAST, '-s')
    hold on
end
xlabel('SNR')
ylabel('BER')
legend('ZF, N = 4', 'ML, N = 4', 'VBLAST, N = 4', 'ZF, N = 8', 'ML, N  = 8', 'VBLAST, N = 8');
toc()
