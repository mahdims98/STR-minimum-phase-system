%% discretization of the model
clear;
sysC  = zpk([-30, -0.6], [-1.3,-2.5,1.1], 0.2);
BW = bandwidth(sysC); 

disc_ratio = 10;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs;

sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v');

% bodeplot(sysC);
% hold on;
% bodeplot(sysD);

B = numD;
A = denD;

%% desired system
overshoot = 10;
settling_time = 3;

zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% sys_2d_desired  = tf([wn^2], [1, 2*zeta*wn, wn^2]);
% damp(sys_desired)
z1 = -15;
z2 = -25;
p3 = -20;

k2 = -p3/(z1*z2);

G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

BmPrime_main = num_desD; 
Am = den_desD;
%%
if B(1) == 0
    B = B(2:end);
end
%% direct with zero cancellation 

% input properties
num_samples = 400;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;


% must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % initial conditions
B_estimated = [0.01, 0.035,0.2];

assert(length(A_estimated)==deg_desA && length(B_estimated)==deg_desB, "initial condistions are not true")

syms q;

% initial B
Bplus = [B_estimated/B_estimated(1)];
Bminus = B_estimated(1);
BmPrime = BmPrime_main/B_estimated(1);



% initial R S T and Ao
R_solved = [1,-0.1,-0.1]; 
S_solved = [1,0.1,0.1];
T = [1.1,0,0];
Ao = [1];

% initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;


theta_real = [A(2:end).'; B.'];

% noise and its degree
noise_poly = [1,0,0].';
deg_noise = length(noise_poly);
noise_variance = 0.01; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);


% RLS initial parameters
rls_solver = RLSClass(100 * eye(total_parameters), 0.1 * ones([total_parameters,1]));

for i = skip_instances:length(uc)
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';     

    %only for simulation. not involved in controller calculations directly
    y(i) = phi_t.' * theta_real + noise_t;

    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';
    theta_hat_new = rls_solver.update_RLS(y(i), phi_t);

    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
%     if B_estimated(1) == 0
%         B_estimated(1) = 1;
%     end
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];

    
    [R_solved,S_solved,T] = solve_diophantin(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);
%     disp(R_solved)
end

figure()
plot(t,y)
hold on;
plot(t, uc, "--r")
% figure()
% plot(t, u)

% sys_new = tf([conv(B,T)],conv(A,T));
