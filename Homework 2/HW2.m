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
    a = -pi/2;
    b = pi/2; 
    theta_i = (b-a).*rand() + a;
    phi_i = (b-a).*rand() + a;

    % Generate arrays as defined by (4) and (5) of the homework specs
    transmitters = 0:1:Nt-1;
    arx = arrayfun(@(L)srv(L, theta_i, Nt), transmitters);
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

% % PROBLEM 3
% % Divide 1000 symbols into a rectangular array, one symbol per column
% s1 = reshape(out.s1,[1,1000]);
% s2 = reshape(out.s2,[1,1000]);
% 
% % Upsample by a factor of 4
% s1u = upsample(s1, 4);
% s2u = upsample(s2, 4);
% 
% % Create a rcos filter using given parameters
% filter = rcosdesign(0.5, 8, 1);
% 
% % Apply filter to each signal
% s1t = upfirdn(s1u, filter);
% s2t = upfirdn(s2u, filter);
% 
% % Convert given anges to radians
% phi1 = deg2rad(30);
% phi2 = deg2rad(40);
% 
% % Create channels h1 and h2, with # of antennae N
% N = 8;
% antennae = 0:1:N-1;
% h1 = arrayfun(@(L)exp(L*-1i*pi*sin(phi1)), antennae);
% h2 = arrayfun(@(L)exp(L*-1i*pi*sin(phi2)), antennae);
% 
% % Apply precoding
% H = [h1
%     h2];
% Hherm = H';
% P = Hherm/(H*Hherm);
% st = [reshape(s1t,[4008, 1]) reshape(s2t,[4008, 1])];
% x = P*st;

function srv = srv(num, ang, ants)
    element = exp(num*-1i*pi*sin(ang));
    srv = 1/sqrt(ants) * element;
end
