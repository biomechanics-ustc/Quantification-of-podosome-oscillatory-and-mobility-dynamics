%% === Clutch Dynamics Simulation ===
clear; clc; close all;

T = 30; dt = 0.1; 
t = 0:dt:T; N_steps = length(t);

R = 0.8; c = 100; k_spring = 0.85;
asymptote_um = 60.0/1000; 

%% --- 1. 50 Springs (150 nm) ---
fprintf('Running Simulation 1 (50 springs - 150nm)...\n');
num_springs_1 = 50; tau_corr_1 = 10; ron_1 = 0.12; roff_1 = 0.12; F_brown_mag_1 = 0.8; 

theta_1 = linspace(0, 2*pi, num_springs_1+1); theta_1 = theta_1(1:end-1);
theta_1_deg = rad2deg(theta_1); 
anchor_points_1 = [R * cos(theta_1); R * sin(theta_1)];
xr0 = -0.25 + 0.2 * sin(2*pi*(1/300.0)*t);

clutch_state_50 = false(num_springs_1, N_steps);
current_pos = [0; 0]; 
spring_connected = rand(1, num_springs_1) < (ron_1/(ron_1+roff_1));
current_brown_force = [0; 0];

exp_factor = exp(-dt/tau_corr_1);
noise_prefactor = F_brown_mag_1 * sqrt(1 - exp_factor^2);

for i = 1:N_steps-1        
    % O-U Brownian force
    white_noise = randn(2, 1);
    current_brown_force = current_brown_force * exp_factor + noise_prefactor * white_noise;
    
    % Markov transitions
    spring_connected = (rand(1, num_springs_1) < (1-exp(-ron_1*dt)) & ~spring_connected) | ...
                       (~(rand(1, num_springs_1) < (1-exp(-roff_1*dt))) & spring_connected);
    clutch_state_50(:, i) = spring_connected; 
    
    % Elastic force
    F_spring_net = [0; 0];
    idx = find(spring_connected);
    if ~isempty(idx)
        r_vecs = current_pos - anchor_points_1(:,idx);
        r_mags = sqrt(sum(r_vecs.^2, 1));
        valid = r_mags > 0;
        if any(valid)
            f_scalar = -k_spring * (r_mags(valid) - (R + xr0(i)));
            forces = f_scalar .* (r_vecs(:,valid) ./ r_mags(valid));
            F_spring_net = sum(forces, 2);
        end
    end
    
    % Velocity saturation
    v_raw = (F_spring_net + current_brown_force) / c;
    v_mag = norm(v_raw);
    if v_mag > 0
        v_soft = asymptote_um * tanh(v_mag / asymptote_um);
        current_vel = v_raw * (v_soft / v_mag);
    else
        current_vel = [0;0];
    end
    current_pos = current_pos + current_vel * dt;
end
clutch_state_50(:, N_steps) = clutch_state_50(:, N_steps-1);

%% --- 2. 150 Springs (30 nm) ---
fprintf('Running Simulation 2 (150 springs - 30nm)...\n');
num_springs_2 = 150; tau_corr_2 = 20; ron_2 = 0.15; roff_2 = 0.1; F_brown_mag_2 = 2.2; 

theta_2 = linspace(0, 2*pi, num_springs_2+1); theta_2 = theta_2(1:end-1);
theta_2_deg = rad2deg(theta_2); 
anchor_points_2 = [R * cos(theta_2); R * sin(theta_2)];

clutch_state_150 = false(num_springs_2, N_steps);
current_pos = [0; 0]; 
spring_connected = rand(1, num_springs_2) < (ron_2/(ron_2+roff_2));
current_brown_force = [0; 0];

exp_factor = exp(-dt/tau_corr_2);
noise_prefactor = F_brown_mag_2 * sqrt(1 - exp_factor^2);

for i = 1:N_steps-1        
    % O-U Brownian force
    white_noise = randn(2, 1);
    current_brown_force = current_brown_force * exp_factor + noise_prefactor * white_noise;
    
    % Markov transitions
    spring_connected = (rand(1, num_springs_2) < (1-exp(-ron_2*dt)) & ~spring_connected) | ...
                       (~(rand(1, num_springs_2) < (1-exp(-roff_2*dt))) & spring_connected);
    clutch_state_150(:, i) = spring_connected; 
    
    % Elastic force
    F_spring_net = [0; 0];
    idx = find(spring_connected);
    if ~isempty(idx)
        r_vecs = current_pos - anchor_points_2(:,idx);
        r_mags = sqrt(sum(r_vecs.^2, 1));
        valid = r_mags > 0;
        if any(valid)
            f_scalar = -k_spring * (r_mags(valid) - (R + xr0(i)));
            forces = f_scalar .* (r_vecs(:,valid) ./ r_mags(valid));
            F_spring_net = sum(forces, 2);
        end
    end
    
    % Velocity saturation
    v_raw = (F_spring_net + current_brown_force) / c;
    v_mag = norm(v_raw);
    if v_mag > 0
        v_soft = asymptote_um * tanh(v_mag / asymptote_um);
        current_vel = v_raw * (v_soft / v_mag);
    else
        current_vel = [0;0];
    end
    current_pos = current_pos + current_vel * dt;
end
clutch_state_150(:, N_steps) = clutch_state_150(:, N_steps-1);

fprintf('All simulations completed. Generating plots...\n');

%% --- 3. Visualization ---
cmap = [0.85, 0.85, 0.85; 0.00, 0.45, 0.74]; 

% Figure 1: 30 nm (Downsampled)
fig1 = figure('Name', 'Clutch States - 30 nm', 'Position', [100, 150, 600, 550], 'Color', 'w');
ax1 = axes('Parent', fig1);

sample_idx = 1:3:num_springs_2;
clutch_sampled_150 = clutch_state_150(sample_idx, :);
theta_sampled_150 = theta_2_deg(sample_idx);

imagesc(ax1, [0, T], [theta_sampled_150(1), theta_sampled_150(end)], clutch_sampled_150);
colormap(ax1, cmap); caxis(ax1, [0 1]);

set(ax1, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 2.0, 'Box', 'on', 'FontSize', 16, 'FontName', 'Times New Roman');
axis(ax1, 'square'); 
xlim(ax1, [0, T]); ylim(ax1, [0, 360]); yticks(ax1, 0:90:360); 

xlabel(ax1, 'Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel(ax1, 'Clutch Position (degree)', 'FontSize', 18, 'FontWeight', 'bold');
title(ax1, '30 nm', 'FontSize', 20, 'FontWeight', 'bold');

hold(ax1, 'on'); 
h1_1 = patch(ax1, [NaN NaN], [NaN NaN], cmap(1,:), 'EdgeColor', 'k');
h1_2 = patch(ax1, [NaN NaN], [NaN NaN], cmap(2,:), 'EdgeColor', 'k');
legend(ax1, [h1_1, h1_2], {'Disconnected', 'Connected'}, 'Location', 'southeast', 'FontSize', 12, 'FontName', 'Times New Roman');

% Figure 2: 150 nm
fig2 = figure('Name', 'Clutch States - 150 nm', 'Position', [720, 150, 600, 550], 'Color', 'w');
ax2 = axes('Parent', fig2);

imagesc(ax2, [0, T], [theta_1_deg(1), theta_1_deg(end)], clutch_state_50);
colormap(ax2, cmap); caxis(ax2, [0 1]);

set(ax2, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 2.0, 'Box', 'on', 'FontSize', 16, 'FontName', 'Times New Roman');
axis(ax2, 'square'); 
xlim(ax2, [0, T]); ylim(ax2, [0, 360]); yticks(ax2, 0:90:360); 

xlabel(ax2, 'Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel(ax2, 'Clutch Position (degree)', 'FontSize', 18, 'FontWeight', 'bold');
title(ax2, '150 nm', 'FontSize', 20, 'FontWeight', 'bold');

hold(ax2, 'on'); 
h2_1 = patch(ax2, [NaN NaN], [NaN NaN], cmap(1,:), 'EdgeColor', 'k');
h2_2 = patch(ax2, [NaN NaN], [NaN NaN], cmap(2,:), 'EdgeColor', 'k');
legend(ax2, [h2_1, h2_2], {'Disconnected', 'Connected'}, 'Location', 'southeast', 'FontSize', 12, 'FontName', 'Times New Roman');

fprintf('Plotting completed successfully.\n');