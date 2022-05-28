clear
% PROBLEM 1.1

N_1 = [16, 64];                                                             %Number of Antennae
K = [1, 2, 4];                                                              %Resolution in bits
theta = -pi/2:0.005:pi/2;                                                   %Angles that we will use to compute beamforming gains
steering = 30;
for n=1:length(N_1)
    args = mod(-pi*sin(deg2rad(steering))*(0:N_1(n)-1), -2*pi);             %Phase shifts
    W = exp(1i*args);
    figure
    for k=1:length(K)
        res = (0:(2^K(k))-1)*((-2*pi)/(2*K(k)));                            %Resolutions
        args_quantized = interp1(res, res, args, 'nearest', 'extrap');      %Quantize phase shifts
        Wk = exp(1i*args_quantized);
        gain = zeros(1,length(theta));
        for t = 1:length(theta)
            a_theta = exp(-1j*pi*(0:N_1(n)-1)'*sin(theta(t)))/sqrt(N_1(n)); %ULA array response 
            gain(t) = abs(Wk*a_theta)^2;                                    %Array gain
        end
        plot(theta, gain, 'LineWidth', 1)
        hold on
    end
    % Now do infinite resolution
    gain = zeros(1,length(theta));
    for t = 1:length(theta)
        a_theta = exp(-1j*pi*(0:N_1(n)-1)'*sin(theta(t)))/sqrt(N_1(n));     %ULA array response 
        gain(t) = abs(W*a_theta)^2;                                         %Array gain
    end
    plot(theta, gain, 'LineWidth', 1)
    legend('k=1', 'k=2', 'k=4', 'Infinite resolution')
    xlabel('Angle (radians)')
    ylabel('Array Gain')
    title('Analog Gain With '+string(N_1(n))+' Antennae');
end

%--------------------------------------------------------------------------
clear
% PROBLEM 1.2

N_1 = [16, 64];                                                             %Number of Antennae
theta = -pi/2:0.005:pi/2;                                                   %Angles that we will use to compute beamforming gains
phi = deg2rad(30);
A = 40;
beta = deg2rad(20);
alpha = 0.5842*(A-21)^0.4 + 0.07886*(A-21);


figure
for n=1:length(N_1)
    gain = zeros(1,length(theta));
    m = [-N_1(n)/2:-1 1:N_1(n)/2];
    num = sqrt(1-(m.^2/(N_1(n)/2)^2));
    w_m = (besseli(0, alpha*num))/(besseli(0, alpha));
    a_m = zeros(length(w_m), 1);                                            %Define a to be a column vector                                     
    for M = 1:length(m)                                                     %Get a_m individually because I don't know how to do it all in one line
        a_m(M) = w_m(M)*exp(-1i*(m(M)-0.5)*pi*sin(phi))*(sin(sin(beta/2)*pi*(m(M)-0.5)))/(pi*(m(M)-0.5));
    end
    for t = 1:length(theta)    
        a_theta = exp(-1j*pi*(0:N_1(n)-1)'*sin(theta(t)))/sqrt(N_1(n));
        gain(t) = abs(a_m'*a_theta)^2;
    end
    plot(theta, gain, 'LineWidth', 1)
    hold on
end
legend('N_T=16', 'N_T=64')
xlabel('Angle (radians)')
ylabel('Array Gain')
title('Digital Gain');

%--------------------------------------------------------------------------
clear

% PROBLEM 2
% Input given parameters
Nt = 8;
Nr = 8;
L = 6;

% Generate channel H_r
H_r = 1/sqrt(2)*(rand(Nt, Nr) +1i*rand(Nt,Nr));

coeff = sqrt((Nt*Nr)/L);
sig = zeros(Nt, Nr);

for i = 1:L
    % Generate random angles between a and b
    theta_i = unifrnd(-pi/2, pi/2);
    phi_i = unifrnd(-pi/2, pi/2);

    % Generate arrays as defined by (4) and (5) of the homework specs
    transmitters = 0:1:Nt-1;
    arx = arrayfun(@(L)srv(L, theta_i, Nt), transmitters);                  %srv defined at the bottom
    receivers = 0:1:Nr-1;
    atx = arrayfun(@(L)srv(L, phi_i, Nt), receivers);
    arx = transpose(arx);
    atx = transpose(atx);

    % Generate random path gain...
    alpha = 1/sqrt(2)*(rand() +1i*rand());
    % ...and multiply it linearly to the result of a_rx*a_tx'
    sig = sig + alpha*(arx*atx');
end

H_s = coeff * sig;

SNR = -10:2:20;
for i = 1:1000
    %Generate random symbols
    b_in = randi([0 1],Nt*log2(4),1);
    x = qammod(b_in, 4, 'InputType', 'bit', 'UnitAveragePower', true);
    for s = length(SNR)
    
        % Zero-forcing combining
        G_ZFHr = H_r'/(H_r'*H_r);
        G_ZFHs = H_s'/(H_s'*H_s);
        
        % LMMSE
        Eb_N0_lin = 10^(SNR(s)/10);
        SNR_lin = 2*Eb_N0_lin;
        noise = sqrt(1/(2*SNR_lin)) * (randn(Nr,1) + 1j*randn(Nr,1));
        G_LMMSEHr = inv(H_r'*H_r+(1/SNR_lin)*eye(Nr))*H_r';
        G_LMMSEHs = inv(H_s'*H_s+(1/SNR_lin)*eye(Nr))*H_s';
    
        % SVD
        [Ur,Sigmar,Vr] = svd(H_r);
        y_precodingr = H_r*Vr*x+noise;
        x_SVDr = Ur'*y_precodingr;
        [Us,Sigmas,Vs] = svd(H_s);
        y_precodings = H_s*Vs*x+noise;
        x_SVDs = Us'*y_precodings;

        c_zf_r = 0;
        c_zf_s = 0;
        for k = 1:length(H_r)
            c_zf_r = c_zf_r + log2(1 + (G_ZFHr(k)'.*x)/((G_ZFHr(k)'.*noise).^2));
        end

        for k = 1:length(H_s)
            c_zf_s = c_zf_s + log2(1 + (G_ZFHs(k)'.*x)/((G_ZFHs(k)'.*noise).^2));
        end

        c_lm_r = 0;
        c_lm_s = 0;
        for k = 1:length(H_r)
            c_lm_r = c_lm_r + log2(1 + (G_LMMSEHr(k)'.*x)/((G_LMMSEHr(k)'.*noise).^2));
        end

        for k = 1:length(H_s)
            c_lm_s = c_lm_s + log2(1 + (G_LMMSEHs(k)'.*x)/((G_LMMSEHs(k)'.*noise).^2));
        end
    end
end

%--------------------------------------------------------------------------
clear

% PROBLEM 3
% Generate random 4-QAM symbols
bits1=randi([0 1],2,1000);
s1=qammod(bits1, 4, 'InputType', 'bit', 'UnitAveragePower', true);
bits2 = zeros(2,1000);
s2=qammod(bits2, 4, 'InputType', 'bit', 'UnitAveragePower', true);

% Upsample by a factor of 4
s1u = upsample(s1, 4);
s2u = upsample(s2, 4);

% Create a rcos filter using given parameters
myfilter = rcosdesign(0.5, 8, 4);

% Apply filter to each signal
s1t = filter(myfilter, 1, s1u);
s2t = filter(myfilter, 1, s2u);

% Convert given anges to radians
phi1 = deg2rad(30);
phi2 = deg2rad(40);

% Create channels h1 and h2, with # of antennae N
N = 32;
antennae = 0:1:N-1;
h1 = arrayfun(@(L)exp(L*-1i*pi*sin(phi1)), antennae);
h2 = arrayfun(@(L)exp(L*-1i*pi*sin(phi2)), antennae);

% Apply precoding
H = [h1.' h2.']';
P = H'/(H*H');
st = [s1t.',s2t.'].';
x = P*st;

% Adjust betas for each section
beta_1 = 1;
beta_3 = -133;

% Calculate z(t) as given in equation (10)
for i = 1:N
    z(i,:) = beta_1*x(i,:) + beta_3.*x(i,:)*(abs(x(i,:)))'.^2;
end

% Graph PST of z(1)
figure
pwelch(z(1,:), [], [], [], 'mean', 'centered');
title('PSD of z_1');

theta = -pi/2:0.005:pi/2;
a_phi = zeros(length(theta),1);
g_phi = zeros(N, length(theta));
g_phi_avg = zeros(length(theta), length(theta));
for t = 1:length(theta)
  a_phi = exp((-1i*pi*(0:N-1))*sin(theta(t)));
  for i = 1:N
        g_phi(i, t) = mean(abs(a_phi'*z(i)).^2);
  end
  g_phi_avg = mean(g_phi);
end
plot(theta, g_phi_avg);

xlabel('Angle (radians)')
ylabel('Gain')
title('Angle vs Interference');




function srv = srv(num, ang, ants)
    element = exp(num*-1i*pi*sin(ang));
    srv = 1/sqrt(ants) * element;
end
