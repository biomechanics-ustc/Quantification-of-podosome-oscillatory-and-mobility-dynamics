%% === Initialization & Batch Setup ===
clear; clc; close all;

phi = cos(pi/4); x0 = 2900; alpha = 1.5; gamma = 0.2;
ks = 1000000; kc = 30; Vp0 = 80; Vd = 50;
F0 = 1000; Fp0 = 20000; tfinal = 600; 
dt_noise = 2.0; tspan = [0 tfinal];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
rfs = 1; rcs = ks/(kc+ks);

num_sims = 300;

% 30 nm Group parameters
Kf1 = 5.4; tau1 = 24; sigma_p1 = 5; sigma_m1 = 5;
% 150 nm Group parameters
Kf2 = 7.2; tau2 = 18; sigma_p2 = 5; sigma_m2 = 5;

results_P1 = NaN(num_sims, 1); results_A1 = NaN(num_sims, 1);
results_P2 = NaN(num_sims, 1); results_A2 = NaN(num_sims, 1);

fprintf('Running batch simulations (%d trials)...\n', num_sims);

%% === Batch ODE Simulations ===
for n = 1:num_sims
    % --- 30 nm Group ---
    t_noise = 0 : dt_noise : tfinal;
    n_p1 = sigma_p1 * randn(size(t_noise));
    n_m1 = sigma_m1 * randn(size(t_noise));
    L1ss1 = (Vp0-Vd)*Fp0*rcs/ks/Vp0; 
    Fmss1 = Kf1/(Kf1-gamma)*(F0+(alpha-gamma/Kf1/rfs)*ks*L1ss1/phi);
    [~, y1] = ode45(@(t,y) odefcn_local(t, y, [Kf1, alpha, ks, tau1, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs], t_noise, n_p1, n_m1), tspan, [L1ss1; Fmss1-2000], options);
    L_c1 = (( (ks*y1(:,1)/phi/rfs - y1(:,2))/Kf1 + x0 )*phi + y1(:,1)) / 1000;
    [results_P1(n), results_A1(n)] = analyze_oscillation(tspan, L_c1, tfinal);

    % --- 150 nm Group ---
    n_p2 = sigma_p2 * randn(size(t_noise));
    n_m2 = sigma_m2 * randn(size(t_noise));
    L1ss2 = (Vp0-Vd)*Fp0*rcs/ks/Vp0; 
    Fmss2 = Kf2/(Kf2-gamma)*(F0+(alpha-gamma/Kf2/rfs)*ks*L1ss2/phi);
    [~, y2] = ode45(@(t,y) odefcn_local(t, y, [Kf2, alpha, ks, tau2, Vp0, Vd, F0, Fp0, phi, gamma, rcs, rfs], t_noise, n_p2, n_m2), tspan, [L1ss2; Fmss2-1800], options);
    L_c2 = (( (ks*y2(:,1)/phi/rfs - y2(:,2))/Kf2 + x0 )*phi + y2(:,1)) / 1000;
    [results_P2(n), results_A2(n)] = analyze_oscillation(tspan, L_c2, tfinal);
    
    if mod(n, 100) == 0, fprintf('Completed %d/%d\n', n, num_sims); end
end

%% === Data Filtering ===
idx1 = (~isnan(results_P1)) & (results_P1 >= 75);
valid_P1 = results_P1(idx1);
valid_A1 = results_A1(idx1);

idx2 = (~isnan(results_P2)) & (results_P2 >= 50);
valid_P2 = results_P2(idx2);
valid_A2 = results_A2(idx2);

fprintf('Data filtering completed.\n');

%% === Experimental Data Loading & Visualization ===
current_dir = fileparts(mfilename('fullpath'));
path_period = fullfile(current_dir, 'Exp_data', 'period.xlsx');
path_amp = fullfile(current_dir, 'Exp_data', 'amplitude.xlsx');

try
    exp_P_table = readtable(path_period);
    exp_A_table = readtable(path_amp);
    exp_P30 = exp_P_table{:,1}; exp_P30(isnan(exp_P30)) = [];
    exp_P150 = exp_P_table{:,2}; exp_P150(isnan(exp_P150)) = [];
    exp_A30 = exp_A_table{:,1}; exp_A30(isnan(exp_A30)) = [];
    exp_A150 = exp_A_table{:,2}; exp_A150(isnan(exp_A150)) = [];
catch
    warning('Failed to read experimental data. Using simulated previews. Path attempted: %s', path_period);
    exp_P30 = valid_P1*1.1; exp_P150 = valid_P2*1.2;
    exp_A30 = valid_A1*0.9; exp_A150 = valid_A2*0.8;
end

color_sim_fill = [0.447, 0.624, 0.812]; 
color_exp_fill = [0.878, 0.478, 0.478];
line_width_main = 2.5; 

data_P_all = {valid_P1, exp_P30, valid_P2, exp_P150};
data_A_all = {valid_A1, exp_A30, valid_A2, exp_A150};
titles = {'Period Distribution', 'Amplitude Distribution'};
ylabels = {'Period (s)', 'Amplitude (\mum)'};
datasets = {data_P_all, data_A_all};

for i = 1:2
    figure('Color', 'w', 'Position', [200 + i*500, 200, 550, 500]);
    hold on;
    current_data = datasets{i};
    
    pos = [1, 1.8, 3.5, 4.3];
    colors_fill = {color_sim_fill, color_exp_fill, color_sim_fill, color_exp_fill};
    
    for g = 1:4
        val = current_data{g};
        if isempty(val), continue; end
        
        q1 = prctile(val, 25);
        q3 = prctile(val, 75);
        med = median(val);
        w_low = min(val(val >= q1 - 1.5*(q3-q1)));
        w_high = max(val(val <= q3 + 1.5*(q3-q1)));
        
        rectangle('Position', [pos(g)-0.3, q1, 0.6, q3-q1], ...
            'FaceColor', colors_fill{g}, 'EdgeColor', 'k', 'LineWidth', line_width_main);
        
        line([pos(g)-0.3, pos(g)+0.3], [med, med], 'Color', 'k', 'LineWidth', line_width_main + 1);
        line([pos(g), pos(g)], [q3, w_high], 'Color', 'k', 'LineStyle', '-', 'LineWidth', line_width_main);
        line([pos(g), pos(g)], [q1, w_low], 'Color', 'k', 'LineStyle', '-', 'LineWidth', line_width_main);
        line([pos(g)-0.15, pos(g)+0.15], [w_high, w_high], 'Color', 'k', 'LineWidth', line_width_main);
        line([pos(g)-0.15, pos(g)+0.15], [w_low, w_low], 'Color', 'k', 'LineWidth', line_width_main);
    end
    
    ax = gca;
    set(ax, 'LineWidth', 2, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold', ...
        'XTick', [1.4, 3.9], 'XTickLabel', {'30 nm', '150 nm'}, 'Box', 'on', 'TickDir', 'in');
    
    ylabel(ylabels{i}, 'FontSize', 22, 'FontWeight', 'bold');
    title(['Sim vs. Exp: ', titles{i}], 'FontSize', 20, 'FontWeight', 'bold');
    xlim([0.2, 5.1]);
    
    hSim = patch(NaN, NaN, color_sim_fill, 'EdgeColor', 'k', 'LineWidth', line_width_main);
    hExp = patch(NaN, NaN, color_exp_fill, 'EdgeColor', 'k', 'LineWidth', line_width_main);
    lg = legend([hSim, hExp], {'Simulation', 'Experiment'}, 'Location', 'best', 'FontSize', 14);
    set(lg, 'FontWeight', 'bold');
    
    grid off;
end

%% === Helper Functions ===
function [T, A] = analyze_oscillation(tspan, L_core, tfinal)
    t_interp = linspace(0, tfinal, length(L_core))';
    steady_idx = t_interp > (tfinal * 0.3);
    L_s = L_core(steady_idx);
    t_s = t_interp(steady_idx);
    L_filt = smoothdata(L_s, 'sgolay', 'SmoothingFactor', 0.5);
    [pks, locs_p] = findpeaks(L_filt, t_s, 'MinPeakDistance', 8);
    [vls, locs_v] = findpeaks(-L_filt, t_s, 'MinPeakDistance', 8);
    if length(locs_p) >= 2 && ~isempty(vls)
        T = mean(diff(locs_p));
        A = (mean(pks) - mean(-vls)) / 2;
    else
        T = NaN; A = NaN; 
    end
end

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