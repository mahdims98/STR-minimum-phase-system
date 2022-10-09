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
if BmPrime_main(1) == 0
    BmPrime_main = BmPrime_main(2:end);
end
%% direct with zero cancellation 
syms q;

main_title = "2-4-with-noise";
% input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
input_noise_variance = 0.00;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
uc(1) = -1;
uc = uc + input_noise;


% must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % initial conditions
B_estimated = [0.01, 0.035,0.2];

assert(length(A_estimated)==deg_desA && length(B_estimated)==deg_desB, "initial condistions are not true")

% desired B
Bplus = [BmPrime_main/BmPrime_main(1)];
Bminus = BmPrime_main(1);
Bm = conv(Bplus, Bminus);
BmPrime = BmPrime_main/BmPrime_main(1);


% initial R S T and Ao
R_estimated = [1,-0.1,-0.1]; 
S_estimated = [1,0.1,0.1];
T_estimated = [1,0,0];
Ao = [1];

% R S and T which were calculated in the last part
R_real = [1,-0.3223,-0.4214];
S_real = [4.9308,-5.0634,1.5007];
T_real = [0.8293,0.2863,-0.002];


% noise and its degree
noise_poly = [1,0.4,0.1].';
deg_noise = length(noise_poly);
noise_variance = 0.001; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

d0 = deg_desA - deg_desB;
[deg_R, deg_Ao] = find_degrees(Am, A, Ao, Bminus, Bplus, BmPrime);
deg_S = deg_R;
deg_T = deg_R; 



f_filter = conv(Am, Ao);

assert(f_filter(1)==1,"AmAo is not monic")

% initial conditions for y
skip_instances = max(deg_desA, deg_desB);

y = [];
y(1:skip_instances) = 0;
ym(1:skip_instances) = 0;
u(1:skip_instances) = 0;
uf(1:skip_instances) = 0;
yf(1:skip_instances) = 0;
ucf(1:skip_instances) = 0;

total_parameters = deg_R+1 + deg_S+1 + deg_T+1;


R_estimated_toPlot = zeros([num_samples, length(R_estimated)]);
S_estimated_toPlot = zeros([num_samples, length(S_estimated)]);
T_toPlot = zeros([num_samples, length(T_estimated)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

theta_real = [A(2:end).'; B.'];
theta_desired = [Am(2:end).'; Bm.'];


theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters+ length(theta_epsilon_zero)), 0.01 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:length(uc)
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].'; % assumed deg_desB = deg_B
    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly; %only for simulation. not involved in controller calculations directly
    
    y(i) = phi_t.' * theta_real+noise_t;
    u(i) = T_estimated * [uc(i:-1:i-(length(T_estimated)-1))].' + S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].';
    u(i) = u(i)./R_estimated(1);

    uf(i) = u(i) - f_filter(2:end) * uf(i-1:-1:i-(length(f_filter)-1)).';
    yf(i) = y(i) - f_filter(2:end) * yf(i-1:-1:i-(length(f_filter)-1)).';
    ucf(i) = uc(i) - f_filter(2:end) * ucf(i-1:-1:i-(length(f_filter)-1)).';

    phi_d0_filtered = [uf(i-d0:-1:i-d0-deg_R), yf(i-d0:-1:i-d0-deg_S), -ucf(i-d0:-1:i-d0-deg_T)].';
    
    phi_t_m = [-y(i-1:-1:i-(deg_desA - 1)), uc(i-1:-1:i-deg_desB)].';
    ym(i) = phi_t_m.' * theta_desired;
    error = y(i) - ym(i);

    theta_hat_new = els_solver.update_ELS(error, phi_d0_filtered);
   
    R_estimated = theta_hat_new(1:deg_R+1).';
    S_estimated = theta_hat_new(deg_R+2:deg_R + deg_S+2).';
    T_estimated = theta_hat_new(deg_R + deg_S+3:deg_R + deg_S+deg_T+3).';


    R_estimated_toPlot(i, :) = R_estimated;
    S_estimated_toPlot(i, :) = S_estimated;
    T_toPlot(i, :) = T_estimated;
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
plot(t, u, 'DisplayName','control input')
xlabel("sample number");
title("control input");
saveas(gcf,'images/q2/' + main_title+ "-y-u" + '.jpeg')
close all

% R 
total_plot_rows = ceil(length(R_estimated_toPlot(1,:))/2) + ceil(length(S_estimated_toPlot(1,:))/2)+ ceil(length(T_toPlot(1,:))/2) - 2;
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_estimated_toPlot(:,i)./R_estimated_toPlot(:,1), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i-1));
    legend('Location','best');
    xlabel("sample number")

end
end_of_r_plot = i-1;

% S
for i = 1:(length(S_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_estimated_toPlot(:,i)./R_estimated_toPlot(:,1), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end

end_of_S_plot = i + end_of_r_plot;
% T
for i = 1:(length(T_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i+end_of_S_plot);
    title_text = "T_%d";
    plot(1:num_samples, T_toPlot(:,i)./R_estimated_toPlot(:,1), 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end
saveas(gcf,'images/q2/' + main_title+ "-T" + '.jpeg')
close all

% Theta
% f5 = figure();
% for i = 1:length(theta_hat_toPlot(1,:))
%     title_text = "θ_%d";
%     subplot(ceil(length(theta_hat_toPlot(1,1:length(theta_real)))/2),2,i);
%     plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Predicted')
% %     hold on;
% %     plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Real')
%     title(sprintf(title_text, i));
%     legend('Location','best');
%     xlabel("sample number")
% end
% saveas(gcf,'images/q2/' + main_title+ "-theta" + '.jpeg')
% close all

