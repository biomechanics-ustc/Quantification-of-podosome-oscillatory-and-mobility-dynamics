%% === Initialization & Parameter Setup ===
clear; clc; close all;

phi = cos(pi/4); x0 = 2900; alpha = 1.5; gamma = 0.2;
ks = 1000000; kc = 30; Vp0 = 80; Vd = 50;
F0 = 1000; Fp0 = 20000; tfinal = 600;

tspan = linspace(0, tfinal, 6000); 

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
dt_noise = 5.0; 

% 30 nm Group parameters
Kf1 = 5.4; tau1 = 24; sp_1 = 5; sm_1 = 5;

% 150 nm Group parameters
Kf2 = 7.2; tau2 = 18; sp_2 = 5; sm_2 = 5;

%% === Group 1: 30 nm Simulation ===
fprintf('Running 30 nm simulation...\n');
rcs = ks/(kc+ks); rfs = 1;
L1ss1 = (Vp0-Vd)*Fp0*rcs/ks/Vp0; 
Fmss1 = Kf1/(Kf1-gamma)*(F0+(alpha-gamma/Kf1/rfs)*ks*L1ss1/phi); 
y0_1 = [L1ss1; Fmss1-2000]; 

t_n1 = 0 : dt_noise : tfinal;
np1 = sp_1 * randn(size(t_n1)); nm1 = sm_1 * randn(size(t_n1));
spara1 = [Kf1, alpha, ks, tau1, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs];

[t1, y1] = ode45(@(t,y) odefcn_local(t, y, spara1, t_n1, np1, nm1), tspan, y0_1, options);
L_core1 = (((ks*y1(:,1)/phi/rfs - y1(:,2))/Kf1 + x0)*phi + y1(:,1)) / 1000;

%% === Group 2: 150 nm Simulation ===
fprintf('Running 150 nm simulation...\n');
L1ss2 = (Vp0-Vd)*Fp0*rcs/ks/Vp0; 
Fmss2 = Kf2/(Kf2-gamma)*(F0+(alpha-gamma/Kf2/rfs)*ks*L1ss2/phi); 
y0_2 = [L1ss2; Fmss2-1200]; 

t_n2 = 0 : dt_noise : tfinal;
np2 = sp_2 * randn(size(t_n2)); nm2 = sm_2 * randn(size(t_n2));
spara2 = [Kf2, alpha, ks, tau2, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs];

[t2, y2] = ode45(@(t,y) odefcn_local(t, y, spara2, t_n2, np2, nm2), tspan, y0_2, options);
L_core2 = (((ks*y2(:,1)/phi/rfs - y2(:,2))/Kf2 + x0)*phi + y2(:,1)) / 1000;

%% === Visualization: Core Length Dynamics ===
figure('Name', 'CoreLength_Comparison', 'Position', [100, 100, 1000, 750], 'Color', 'w');
set(gcf, 'GraphicsSmoothing', 'on'); 
hold on;

c_dark_blue = [0.00, 0.20, 0.60];   
c_light_blue = [0.30, 0.75, 0.93];  

plot(t2, L_core2, 'Color', c_light_blue, 'LineWidth', 5, 'DisplayName', '150 nm Group');
plot(t1, L_core1, 'Color', c_dark_blue, 'LineWidth', 5, 'DisplayName', '30 nm Group');

ax = gca;
set(ax, 'LineWidth', 3.5, 'FontSize', 28, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'bold', 'TickDir', 'in', 'TickLength', [0.03, 0.03], 'Box', 'on');
set(gcf, 'Renderer', 'painters'); 

xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 36, 'FontWeight', 'bold');
ylabel('Core Length (\mum)', 'FontName', 'Times New Roman', 'FontSize', 36, 'FontWeight', 'bold');

lgd = legend('show', 'Location', 'best');
set(lgd, 'FontName', 'Times New Roman', 'FontSize', 26, 'FontWeight', 'bold');

xlim([0, tfinal]);
y_all = [L_core1; L_core2];
ylim([min(y_all) - (max(y_all)-min(y_all))*0.15, max(y_all) + (max(y_all)-min(y_all))*0.15]);

grid off;

%% === Feature Extraction & Output ===
[p1, a1] = extract_features(t1, L_core1, tfinal);
[p2, a2] = extract_features(t2, L_core2, tfinal);

fprintf('\n================ Simulation Report ================ \n');
fprintf('30 nm Group  | Period: %.2f s | Amplitude: %.4f um\n', p1, a1);
fprintf('150 nm Group | Period: %.2f s | Amplitude: %.4f um\n', p2, a2);
fprintf('=================================================== \n');

%% === Helper Function: ODE System ===
function dydt = odefcn_local(t, y, spara, t_n, n_p, n_m)
    Kf = spara(1); alpha = spara(2); k = spara(3); tau = spara(4); 
    Vp0 = spara(5); Vd = spara(6); F0 = spara(7); Fp0 = spara(8); 
    phi = spara(9); gamma = spara(10); rcs = spara(11); rfs = spara(12);
    
    chi_p = interp1(t_n, n_p, t, 'linear');
    chi_m = interp1(t_n, n_m, t, 'linear');
    
    dydt = zeros(2,1);
    dFmdt = (F0 + (alpha - gamma/Kf/rfs)*k*y(1)/phi - (1 - gamma/Kf)*y(2) + chi_m)/tau;
    dydt(1) = 1/(k/Kf/rfs + 1) * (Vp0*(1 - k*y(1)/Fp0/rcs) - Vd + chi_p + phi/Kf*dFmdt);
    dydt(2) = dFmdt; 
end

%% === Helper Function: Feature Extraction ===
function [period, amp] = extract_features(t, L, tfinal)
    steady_idx = t > (tfinal * 0.4);
    ts = t(steady_idx); Ls = L(steady_idx);
    
    L_filt = smoothdata(Ls, 'sgolay', 'SmoothingFactor', 0.1);
    [pks, locs_p] = findpeaks(L_filt, ts, 'MinPeakDistance', 5);
    [vls, ~] = findpeaks(-L_filt, ts, 'MinPeakDistance', 5);
    
    if length(locs_p) >= 2
        period = mean(diff(locs_p));
        amp = (mean(pks) - mean(-vls)) / 2;
    else
        period = NaN; amp = NaN;
    end
end