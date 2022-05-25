
% % PROBLEM 2
% % Input given parameters
% Nt = 8;
% Nr = 8;
% L = 6;
% 
% % Generate channel H_r
% H_r = 1/sqrt(2)*(rand(Nt, Nr) +1i*rand(Nt,Nr));
% 
% coeff = sqrt((Nt*Nr)/L);
% sig = zeros(Nt, Nr);
% 
% for i = 1:L
%     % Generate random angles between a and b
%     theta_i = unifrd(-pi/2, pi/2);
%     phi_i = unifrd(-pi/2, pi/2);
% 
%     % Generate arrays as defined by (4) and (5) of the homework specs
%     transmitters = 0:1:Nt-1;
%     arx = arrayfun(@(L)srv(L, theta_i, Nt), transmitters);
%     receivers = 0:1:Nr-1;
%     atx = arrayfun(@(L)srv(L, phi_i, Nt), receivers);
%     arx = transpose(arx);
%     atx = transpose(atx);
% 
%     % Generate random path gain...
%     alpha = 1/sqrt(2)*(rand() +1i*rand());
%     % ...and multiply it linearly to the result of a_rx*a_tx'
%     sig = sig + alpha*(arx*atx');
% end
% 
% H_s = coeff * sig;

% PROBLEM 3
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
N = 8;
antennae = 0:1:N-1;
h1 = arrayfun(@(L)exp(L*-1i*pi*sin(phi1)), antennae);
h2 = arrayfun(@(L)exp(L*-1i*pi*sin(phi2)), antennae);

% Apply precoding
H = [h1.' h2.'].';
Hherm = H';
P = Hherm/(H*Hherm);
st = [s1t.',s2t.'].';
x = P*st;

%TODO: Fill out Z
beta_1 = 1;
beta_3 = 0;
% z = beta_1*x + beta_3*x.*(abs(x)).^2;

for i = 1:8
    z(i,:) = beta_1*x(i,:) + beta_3.*x(i,:)*(abs(x(i,:)))'.^2;
end

pwelch(z(1,:), [], [], [], 'mean', 'centered');

%TODO: get G
% a_phi = exp(-j.*pi.*((0:N1-1)-1).*sin(phi_radian(i))).';
% g_phi(1,i) =mean(abs(a_phi'*zn).^2);


function srv = srv(num, ang, ants)
    element = exp(num*-1i*pi*sin(ang));
    srv = 1/sqrt(ants) * element;
end
