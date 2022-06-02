clear
N = [4 8];          % No. of Receivers = No. of Transmitters; we'll consider 4 of each and 8 of each
SNR = -20:1:20;     % SNR will range from -20 to 20 dB
iter = 1000;        % Set the number of Monte-Carlo iterations
[avgBER_ZF, avgBER_VBLAST, avgBER_ML] = deal(zeros(size(SNR)));   %Initialize

figure
for n = 1:length(N)                                             % Loop through all N values
    for s = 1:length(SNR)                                       % Loop through all SNR values
        SNR_lin = 10^(SNR(s)/10);                               % Convert SNR to linear scale
        H = 1/sqrt(2)*(rand(N(n), N(n)) +1i*rand(N(n),N(n)));   % Generate channel H
        [iterBER_ZF, iterBER_VBLAST, iterBER_ML] = deal(zeros(size(1, iter)));
        for i = 1:iter
%             % This is a crude "progress bar", prints status every 1000 iterations
%             if mod(i,1000) == 0
%                 disp('Currently at #'+string(i)+' at SNR '+ string(SNR(s)));
%             end

            b_in = randi([0 1],N(n)*log2(4),1);                                 % Generate random signal
            x = qammod(b_in, 4, 'InputType', 'bit', 'UnitAveragePower', true);  % Shape signal to 4-QAM
            noise = sqrt(1/(2*SNR_lin)) * (randn(N(n),1) + 1j*randn(N(n),1));   % Generate complex noise samples with noise power = 1/SNR_linear
            y = H*x+noise;                                                      % Received signal

            % Zero-Forcing (ZF)
            G_ZF = inv(H'*H)*H';
            x_ZF = G_ZF*y;
            bOut_ZF = qamdemod(x_ZF, 4,'OutputType','bit','UnitAveragePower',true);
            iterBER_ZF(i) = sum(bOut_ZF ~= b_in)/length(b_in);

            % V-BLAST
            %CODE HERE

            % Maximum Likelihood (ML)
            %CODE HERE

        end
        avgBER_ZF(s) = mean(iterBER_ZF);
    end
    semilogy(SNR, avgBER_ZF)
    hold on
end
xlabel('SNR')
ylabel('BER')
legend('N = 4', 'N = 8')