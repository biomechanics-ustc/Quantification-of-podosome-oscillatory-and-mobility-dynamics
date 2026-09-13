%% === Initialization & Data Loading ===
clear; clc; close all;

current_dir = fileparts(mfilename('fullpath'));
path_speed = fullfile(current_dir, 'Exp_data', 'speed.xlsx');

try
    exp_raw = readmatrix(path_speed);
    exp_data{1} = exp_raw(:, 1); % 30 nm Group
    exp_data{2} = exp_raw(:, 2); % 150 nm Group
    exp_data{1}(isnan(exp_data{1})) = [];
    exp_data{2}(isnan(exp_data{2})) = [];
catch
    error('Failed to read Excel file. Path attempted: %s', path_speed);
end
c = 100; k_spring = 0.85;
asymptote_um = 60.0/1000; 
dt = 1.0; T = 600; 
t = 0:dt:T; N_steps = length(t);
xr0 = 0.25 + 0.2 * sin(2*pi*(1/300.0)*t);

N_trials = 100; 

% 30 nm Group
params(1).R = 0.8;      
params(1).num_s = 150; params(1).ron = 0.15; params(1).roff = 0.1;  
params(1).Fb = 2.2; params(1).tau_c = 20;
params(1).color = [0.00, 0.45, 0.74]; 
params(1).name = '30 nm Group';

% 150 nm Group
params(2).R = 0.8;      
params(2).num_s = 50;  params(2).ron = 0.12; params(2).roff = 0.12; 
params(2).Fb = 0.8; params(2).tau_c = 10;
params(2).color = [0.85, 0.33, 0.10]; 
params(2).name = '150 nm Group';

speed_data_sim = cell(1, 2);  

%% === Physical Simulation Loop ===
fprintf('Running simulations (%d trials/group)...\n', N_trials);

for g = 1:2
    R = params(g).R; 
    num_springs = params(g).num_s;
    ron = params(g).ron; roff = params(g).roff;
    F_brown_mag = params(g).Fb;
    tau_corr = params(g).tau_c;
    
    all_speeds_group = []; 
    for trial = 1:N_trials
        pos = zeros(2, N_steps);
        vel_mags = zeros(1, N_steps-1);
        theta = linspace(0, 2*pi, num_springs+1); theta = theta(1:end-1);
        anchor_points = [R * cos(theta); R * sin(theta)];
        spring_connected = rand(1, num_springs) < (ron/(ron+roff));
        current_brown_force = [0; 0]; exp_factor = exp(-dt/tau_corr);
        noise_prefactor = F_brown_mag * sqrt(1 - exp_factor^2); 

        for i = 1:N_steps-1        
            % O-U Brownian force update
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
                v_final = v_raw * (v_soft / v_mag);
                vel_mags(i) = v_soft;
                pos(:, i+1) = pos(:, i) + v_final * dt; 
            else
                pos(:, i+1) = pos(:, i);
            end
        end
        all_speeds_group = [all_speeds_group, vel_mags];
    end
    speed_data_sim{g} = all_speeds_group * 1000; % Convert to nm/s
end

%% === Visualization: Speed Distribution Comparison ===
for g = 1:2
    figure('Name', ['SpeedComp_', params(g).name], 'Color', 'w', 'Position', [100+ (g-1)*950, 100, 1000, 850]);
    hold on;
    
    sim_vals = speed_data_sim{g};
    exp_vals = exp_data{g};
    
    sim_color = params(g).color;
    exp_color = params(g).color * 0.5; 
    
    edges = linspace(0, max([max(sim_vals), max(exp_vals)]), 40); 
    
    % Simulation histogram
    h_sim = histogram(sim_vals, 'BinEdges', edges, 'Normalization', 'pdf', ...
        'FaceColor', sim_color, 'EdgeColor', sim_color * 0.7, ...
        'FaceAlpha', 0.45, 'LineWidth', 1.5, 'DisplayName', 'Simulation');
    
    % Experiment histogram
    h_exp = histogram(exp_vals, 'BinEdges', edges, 'Normalization', 'pdf', ...
        'DisplayStyle', 'stairs', 'EdgeColor', exp_color, ...
        'LineWidth', 5, 'DisplayName', 'Experiment');
    
    % Mean value indicators
    yl = ylim; 
    m_sim = mean(sim_vals); m_exp = mean(exp_vals);
    
    line([m_sim m_sim], [0 yl(2)*0.98], 'Color', sim_color, 'LineStyle', ':', 'LineWidth', 2);
    plot(m_sim, yl(2)*0.98, 'v', 'MarkerSize', 15, 'MarkerFaceColor', sim_color, 'MarkerEdgeColor', 'w');
    
    line([m_exp m_exp], [0 yl(2)*0.98], 'Color', exp_color, 'LineStyle', ':', 'LineWidth', 2);
    plot(m_exp, yl(2)*0.98, 'v', 'MarkerSize', 15, 'MarkerFaceColor', exp_color, 'MarkerEdgeColor', 'w');
    
    % Statistical annotations
    text_str = { ...
        ['$\mu_{\mathrm{Sim}} = ', num2str(m_sim, '%.1f'), '$ nm/s'], ...
        ['$\mu_{\mathrm{Exp}} = ', num2str(m_exp, '%.1f'), '$ nm/s'] ...
    };
    annotation('textbox', [0.55, 0.65, 0.35, 0.15], 'String', text_str, ...
        'Interpreter', 'latex', 'FontSize', 34, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'Color', 'k', 'HorizontalAlignment', 'left');

    % Axes formatting
    set(gca, 'LineWidth', 3, 'FontSize', 32, 'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', 'TickDir', 'in', 'Box', 'on');
    xlabel('Speed (nm/s)', 'FontSize', 44, 'FontWeight', 'bold');
    ylabel('Probability Density', 'FontSize', 44, 'FontWeight', 'bold');
    title(params(g).name, 'FontSize', 40, 'FontWeight', 'bold');
    
    lgd = legend([h_sim, h_exp], {'Simulation', 'Experiment'}, 'Location', 'northeast', ...
        'FontSize', 30, 'Box', 'off');
    
    xlim([0, edges(end)*1.05]); 
    ylim([0, yl(2)*1.15]); 
    grid off; hold off;
end

fprintf('\n--- Speed Statistics Summary ---\n');
for g = 1:2
    fprintf('%s: Sim Mean = %.1f nm/s | Exp Mean = %.1f nm/s\n', ...
        params(g).name, mean(speed_data_sim{g}), mean(exp_data{g}));
end