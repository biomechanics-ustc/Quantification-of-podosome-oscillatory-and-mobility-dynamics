%% === Initialization & Parameter Setup ===
clear; clc; close all;

current_dir = fileparts(mfilename('fullpath'));
path_msd = fullfile(current_dir, 'Exp_data', 'MSD.xlsx');

try
    exp_alpha_raw = readmatrix(path_msd);
    exp_30_alpha = exp_alpha_raw(:,1); exp_30_alpha(isnan(exp_30_alpha)) = [];
    exp_150_alpha = exp_alpha_raw(:,2); exp_150_alpha(isnan(exp_150_alpha)) = [];
catch
    error('Failed to read MSD data. Please check the file path.');
end

R = 0.8; c = 100; k_spring = 0.85; asymptote_um = 60.0/1000; 
dt = 1.0; T = 600; t = 0:dt:T; N_steps = length(t);
xr0 = 0.25 + 0.2 * sin(2*pi*(1/300.0)*t);
N_trials = 100; 

color_sim = [0.10, 0.40, 0.80]; 
color_exp = [0.85, 0.15, 0.15]; 

sim_alpha_all = cell(1,2);

%% === Monte Carlo Physical Simulations ===
fprintf('Running statistical simulations (%d trials/group)...\n', N_trials);

for g = 1:2
    if g == 1 % 30 nm 
        num_springs = 150; ron = 0.15; roff = 0.1; Fb = 2.0; tau_c = 20;
    else      % 150 nm 
        num_springs = 50; ron = 0.12; roff = 0.12; Fb = 0.8; tau_c = 10;
    end
    
    alpha_list = zeros(N_trials, 1);
    for trial = 1:N_trials
        pos = zeros(2, N_steps);
        theta = linspace(0, 2*pi, num_springs+1); theta = theta(1:end-1);
        anchor_points = [R * cos(theta); R * sin(theta)];
        spring_connected = rand(1, num_springs) < (ron/(ron+roff));
        
        current_brown_force = [0; 0]; exp_factor = exp(-dt/tau_c);
        noise_prefactor = Fb * sqrt(1 - exp_factor^2); 

        for i = 1:N_steps-1        
            % O-U Brownian force
            white_noise = randn(2, 1);
            current_brown_force = current_brown_force * exp_factor + noise_prefactor * white_noise;
            
            % Markov transitions
            spring_connected = (rand(1, num_springs) < (1-exp(-ron*dt)) & ~spring_connected) | ...
                               (~(rand(1, num_springs) < (1-exp(-roff*dt))) & spring_connected);
            
            % Elastic force
            F_spring_net = [0; 0];
            idx = find(spring_connected);
            if ~isempty(idx)
                r_vecs = pos(:,i) - anchor_points(:,idx);
                r_mags = sqrt(sum(r_vecs.^2, 1));
                valid = r_mags > 0;
                if any(valid)
                    forces = -k_spring * (r_mags(valid) - (R + xr0(i))) .* (r_vecs(:,valid) ./ r_mags(valid));
                    F_spring_net = sum(forces, 2);
                end
            end
            
            % Velocity saturation
            v_raw = (F_spring_net + current_brown_force) / c;
            v_mag = norm(v_raw);
            if v_mag > 0
                v_soft = asymptote_um * tanh(v_mag / asymptote_um);
                pos(:, i+1) = pos(:, i) + (v_raw * (v_soft / v_mag)) * dt; 
            else
                pos(:, i+1) = pos(:, i);
            end
        end
        
        % Time-averaged MSD calculation
        max_tau = floor(N_steps / 3); 
        tau_idx = unique(round(logspace(0, log10(max_tau), 25))); 
        msd = zeros(length(tau_idx), 1);
        tau_t = tau_idx' * dt;
        for k = 1:length(tau_idx)
            tau = tau_idx(k);
            msd(k) = mean(sum((pos(:, 1+tau:end) - pos(:, 1:end-tau)).^2, 1));
        end
        
        % Power-law fitting for diffusion exponent alpha
        p = polyfit(log(tau_t(tau_t>5)), log(msd(tau_t>5)), 1);
        alpha_list(trial) = p(1);
    end
    sim_alpha_all{g} = alpha_list;
end

%% === Visualization: Grouped Boxplot & Jittered Scatter ===
figure('Name', 'Grouped_Boxplot_Final', 'Color', 'w', 'Position', [100, 100, 950, 800]);
hold on;

dx = 0.18; 
pos_list = [1-dx, 1+dx, 2-dx, 2+dx];
data_cell = {sim_alpha_all{1}, exp_30_alpha, sim_alpha_all{2}, exp_150_alpha};
plot_colors = {color_sim, color_exp, color_sim, color_exp};

% Jittered scatter
for k = 1:4
    curr_data = data_cell{k};
    n_pts = length(curr_data);
    x_jitter = pos_list(k) + (rand(n_pts, 1) - 0.5) * 0.12; 
    
    scatter(x_jitter, curr_data, 70, ...
        'MarkerEdgeColor', plot_colors{k} * 0.8, ... 
        'MarkerFaceColor', plot_colors{k}, ...
        'MarkerFaceAlpha', 0.45, ...  
        'LineWidth', 0.8);
end

% Boxplots
for k = 1:4
    h = boxplot(data_cell{k}, 'Positions', pos_list(k), 'Widths', 0.25, ...
        'Symbol', '', 'Colors', 'k'); 
    
    set(h, 'LineWidth', 2.5); 
    
    box_obj = findobj(h, 'Tag', 'Box');
    patch(get(box_obj, 'XData'), get(box_obj, 'YData'), plot_colors{k}, ...
        'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 2.5);
    
    median_obj = findobj(h, 'Tag', 'Median');
    set(median_obj, 'Color', 'k', 'LineWidth', 3.5);
end

% Axes formatting
set(gca, 'LineWidth', 4.0, 'FontSize', 32, 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold', 'TickDir', 'in', 'Box', 'on', ...
    'XTick', [1, 2], 'XTickLabel', {'30 nm', '150 nm'});

ylabel('Diffusion Exponent \alpha', 'FontSize', 42, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
title('Sim vs. Exp: \alpha Distribution', 'FontSize', 36, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

h_sim_leg = patch(NaN, NaN, color_sim, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 2);
h_exp_leg = patch(NaN, NaN, color_exp, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 2);
legend([h_sim_leg, h_exp_leg], {'Simulation', 'Experiment'}, ...
    'Location', 'NorthWest', 'FontSize', 30, 'Box', 'off');

ylim([0.2, 1.5]);
xlim([0.5, 2.5]);

grid off;
hold off;

fprintf('Plotting completed successfully.\n');