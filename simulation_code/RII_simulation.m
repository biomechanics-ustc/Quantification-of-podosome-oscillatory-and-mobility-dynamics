%% === Initialization & Parameter Setup ===
c = 100; k_spring = 0.85;
asymptote_um = 60.0/1000; 
dt = 1.0; T = 600; 
t = 0:dt:T; N_steps = length(t);
xr0 = 0.25 + 0.2 * sin(2*pi*(1/300.0)*t);

N_trials = 20; 

% 30 nm Group 
params(1).num_s = 150; params(1).ron = 0.15; params(1).roff = 0.1;  
params(1).Fb = 2.2; params(1).tau_c = 20; params(1).R = 0.8;
params(1).color = [0.10, 0.40, 0.80]; 
params(1).name = '30 nm Group';

% 150 nm Group
params(2).num_s = 50;  params(2).ron = 0.12; params(2).roff = 0.12; 
params(2).Fb = 0.8; params(2).tau_c = 10; params(2).R = 0.8;
params(2).color = [0.85, 0.15, 0.15];
params(2).name = '150 nm Group';

integrity_data = cell(1, 2);  

%% === Physical Simulation Loop ===
fprintf('Running simulations (N_trials = %d)...\n', N_trials);

for g = 1:2
    num_springs = params(g).num_s;
    ron = params(g).ron; roff = params(g).roff;
    F_brown_mag = params(g).Fb;
    tau_corr = params(g).tau_c;
    R = params(g).R;
    all_integrity_samples = []; 
    
    for trial = 1:N_trials
        pos = zeros(2, N_steps);
        spring_connected = rand(1, num_springs) < (ron/(ron+roff));
        current_brown_force = [0; 0]; exp_factor = exp(-dt/tau_corr);
        noise_prefactor = F_brown_mag * sqrt(1 - exp_factor^2); 
        trial_integrity = zeros(1, N_steps);

        for i = 1:N_steps-1        
            % Markov transitions
            spring_connected = (rand(1, num_springs) < (1-exp(-ron*dt)) & ~spring_connected) | ...
                               (~(rand(1, num_springs) < (1-exp(-roff*dt))) & spring_connected);
            trial_integrity(i) = sum(spring_connected) / num_springs;

            % Elastic force
            F_spring_net = [0; 0];
            idx = find(spring_connected);
            if ~isempty(idx)
                theta = linspace(0, 2*pi, num_springs+1); theta = theta(1:end-1);
                anchor_points = [R * cos(theta); R * sin(theta)];
                r_vecs = pos(:,i) - anchor_points(:,idx);
                r_mags = sqrt(sum(r_vecs.^2, 1));
                if any(r_mags > 0)
                    valid = r_mags > 0;
                    forces = -k_spring * (r_mags(valid) - (R + xr0(i))) .* (r_vecs(:,valid) ./ r_mags(valid));
                    F_spring_net = sum(forces, 2);
                end
            end
            
            % O-U Brownian force update
            current_brown_force = current_brown_force * exp_factor + noise_prefactor * randn(2, 1);
            v_raw = (F_spring_net + current_brown_force) / c;
            v_mag = norm(v_raw);
            
            % Velocity saturation
            if v_mag > 0
                v_soft = asymptote_um * tanh(v_mag / asymptote_um);
                pos(:, i+1) = pos(:, i) + v_raw * (v_soft / v_mag) * dt;
            else
                pos(:, i+1) = pos(:, i);
            end
        end
        steady_idx = round(N_steps * 0.2) : N_steps;
        all_integrity_samples = [all_integrity_samples, trial_integrity(steady_idx)];
    end
    integrity_data{g} = all_integrity_samples;
end

%% === Visualization: Enhanced Boxplot with Jittered Scatter ===
figure('Color', 'w', 'Position', [100, 100, 950, 800]); 
hold on;

line_width_main = 2.5; 
box_width = 0.45;
pos_x = [1.2, 2.8]; 
h_for_legend = gobjects(1, 2); 

for g = 1:2
    val = integrity_data{g};
    
    % Jittered scatter
    sample_rate = 120; 
    val_sampled = val(1:sample_rate:end);
    x_jitter = pos_x(g) + (rand(size(val_sampled)) - 0.5) * 0.25;
    
    scatter(x_jitter, val_sampled, 70, ...
        'MarkerEdgeColor', params(g).color * 0.8, ... 
        'MarkerFaceColor', params(g).color, ...
        'MarkerFaceAlpha', 0.45, ...  
        'LineWidth', 0.8); 
    
    % Boxplot statistics
    q1 = prctile(val, 25);
    q3 = prctile(val, 75);
    med = median(val);
    iqr_val = q3 - q1;
    w_low = min(val(val >= q1 - 1.5*iqr_val));
    w_high = max(val(val <= q3 + 1.5*iqr_val));
    
    rectangle('Position', [pos_x(g)-box_width/2, q1, box_width, q3-q1], ...
        'FaceColor', [params(g).color, 0.6], 'EdgeColor', 'k', 'LineWidth', line_width_main);
    
    h_for_legend(g) = patch(NaN, NaN, params(g).color, 'FaceAlpha', 0.6, ...
        'EdgeColor', 'k', 'LineWidth', 2, 'DisplayName', params(g).name);
    
    line([pos_x(g)-box_width/2, pos_x(g)+box_width/2], [med, med], ...
        'Color', 'k', 'LineWidth', 3.5);
    
    line([pos_x(g), pos_x(g)], [q3, w_high], 'Color', 'k', 'LineWidth', line_width_main);
    line([pos_x(g), pos_x(g)], [q1, w_low], 'Color', 'k', 'LineWidth', line_width_main);
    line([pos_x(g)-0.12, pos_x(g)+0.12], [w_high, w_high], 'Color', 'k', 'LineWidth', line_width_main);
    line([pos_x(g)-0.12, pos_x(g)+0.12], [w_low, w_low], 'Color', 'k', 'LineWidth', line_width_main);
end

% Axes formatting
ax = gca;
set(ax, 'LineWidth', 4.0, 'FontSize', 32, 'FontName', 'Times New Roman', ...
    'FontWeight', 'bold', 'TickDir', 'in', 'Box', 'on', ...
    'XTick', [1.2, 2.8], 'XTickLabel', {'30 nm', '150 nm'});

ylabel('Bound Linker Ratio', 'FontSize', 42, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
title('Ring Integrity', 'FontSize', 36, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

lgd = legend(h_for_legend, {'30 nm Group', '150 nm Group'}, ...
    'Location', 'NorthEast', 'FontSize', 30, 'Box', 'off');

xlim([0.5, 3.5]);
all_vals = [integrity_data{1}, integrity_data{2}];
ylim([min(all_vals)*0.85, max(all_vals)*1.15]);

grid off; hold off;

fprintf('Simulation and plotting completed successfully.\n');