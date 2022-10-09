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
main_title = "2-3-with-noise-with-disturb-with-integ";

% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
input_noise_variance = 0;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
uc(1) = -1;
uc = uc + input_noise;



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

% R S and T which were calculated in the last part
R_real = [1,-0.3223,-0.4214];
S_real = [4.9308,-5.0634,1.5007];
T_real = [0.8293,0.2863,-0.002];


% initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;


theta_real = [A(2:end).'; B.'];

% noise and its degree
noise_poly = [1,0.4,0.2].';
deg_noise = length(noise_poly);
noise_variance = 0.001; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

% disturbance
v = [zeros([1,ceil(num_samples/2)]), 1*ones([1,ceil(num_samples/2)])];

%integrator
X = [1, 0];

% for plotting
u_toPlot = zeros([num_samples, length(uc)]);
R_solved_toPlot = zeros([num_samples, length(R_solved)]);
S_solved_toPlot = zeros([num_samples, length(S_solved)]);
T_toPlot = zeros([num_samples, length(T)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
theta_epsilon_zero = [1,0,0].';

els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);


for i = skip_instances:length(uc)
    Y = -((1+X(end)) * polyval(R_solved,1))/polyval(B_estimated, 1);
    

    R_estimated_withI = conv(X, R_solved) + [0,conv(Y, B_estimated)];
    S_estimated_withI = conv(X, S_solved) - conv(Y, A_estimated);

    phi_t_real = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';  
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u_toPlot(i-1:-1:i-deg_desB)].'; 
    %only for simulation. not involved in controller calculations directly
    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly;
    y(i) = phi_t_real.' * theta_real + noise_t;
    
    u(i) = T * [uc(i-1:-1:i-(length(T)))].' + S_estimated_withI * [-y(i:-1:i-(length(S_estimated_withI)-1))].' ...
        - R_estimated_withI(2:end) * [u_toPlot(i-1:-1:i-(length(R_estimated_withI)-1))].' + v(i); % it was better to add v in the phi_t.then it was not necessary to declare u_toPlot
    u_toPlot(i) = u(i) - v(i); % pure input
    theta_hat_new = els_solver.update_ELS(y(i), phi_t);


    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    A_estimated = [1, A_estimated];
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];

    [R_solved,S_solved,T] = solve_diophantin(Am, A_estimated, Ao, Bminus, Bplus, BmPrime);

    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
    T_toPlot(i, :) = T;
%     disp(R_solved)
end

%% plotters

figure()
subplot(2,1,1);
plot(t,y, 'DisplayName','Real')
hold on;
plot(t, uc, "--r", 'DisplayName','desired')
xlabel("sample number");
title("output");
legend('Location','best');
subplot(2,1,2);
plot(t, u_toPlot, 'DisplayName','control input')
xlabel("sample number");
title("control input");
saveas(gcf,'images/q2/' + main_title+ "-y-u" + '.jpeg')
close all

% R 
total_plot_rows = ceil(length(R_solved)/2) + ceil(length(S_solved)/2)+ ceil(length(T)/2) - 1;
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_solved)) % first parameter of R is 1
    subplot(total_plot_rows,2,i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_solved_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i-1));
    legend('Location','best');
    xlabel("sample number")

end
end_of_r_plot = i-1;

% S
for i = 1:(length(S_solved)) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_solved_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end

end_of_S_plot = i + end_of_r_plot;
% T
for i = 1:(length(T)) % first parameter of R is 1
    subplot(total_plot_rows,2,i+end_of_S_plot);
    title_text = "T_%d";
    plot(1:num_samples, T_toPlot(:,i), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,'images/q2/' + main_title+ "-T" + '.jpeg')
close all

% Theta
f5 = figure();
for i = 1:length(theta_hat_toPlot(1,:))
    title_text = "θ_%d";
    subplot(ceil(length(theta_hat_toPlot(1,:)))/2,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,'images/q2/' + main_title+ "-theta" + '.jpeg')
close all