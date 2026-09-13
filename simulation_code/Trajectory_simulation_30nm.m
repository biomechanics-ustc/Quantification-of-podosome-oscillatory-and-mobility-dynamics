%% === Initialization & Physical Simulation ===
clear; clc;

R = 0.8; c = 100; k_spring = 0.85;
asymptote_um = 60.0/1000; 
tau_corr = 20; num_springs = 150; ron = 0.15; roff = 0.1;
F_brown_mag = 2.2; dt = 1.0; T = 600; t = 0:dt:T; N_steps = length(t);

num_instances = 4;      
spacing = 1.2;          
pos_all = zeros(2, N_steps, num_instances); 

theta = linspace(0, 2*pi, num_springs+1); theta = theta(1:end-1);
anchor_points_base = [R * cos(theta); R * sin(theta)];
xr0 = 0.25 + 0.2 * sin(2*pi*(1/300.0)*t);

for inst = 1:num_instances
    current_pos = [0; 0]; 
    spring_connected = rand(1, num_springs) < (ron/(ron+roff));
    current_brown_force = [0; 0];
    exp_factor = exp(-dt/tau_corr);
    noise_prefactor = F_brown_mag * sqrt(1 - exp_factor^2);
    
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
            r_vecs = current_pos - anchor_points_base(:,idx);
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
            current_vel = v_raw * (v_soft / v_mag);
        else
            current_vel = [0;0];
        end
        
        current_pos = current_pos + current_vel * dt;
        pos_all(:, i+1, inst) = current_pos;
    end
end

%% === Visualization: Migration Paths ===
figure('Name', 'Podosome_Migration_Paths', 'Position', [100, 300, 900, 380], 'Color', 'w');
hold on;

title('Simulated Migration Paths of Multiple Podosomes', ...
    'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

set(gca, 'LineWidth', 2.2, 'Box', 'on', ...
         'XTick', [], 'YTick', [], ...
         'XColor', 'k', 'YColor', 'k');

for inst = 1:num_instances
    x_offset = (inst-1) * spacing;
    
    this_x = pos_all(1,:,inst) + x_offset;
    this_y = pos_all(2,:,inst);
    surface([this_x; this_x], [this_y; this_y], zeros(2, N_steps), [t; t], ...
            'FaceColor', 'no', 'EdgeColor', 'interp', 'LineWidth', 2.5);
    
    scatter(this_x(end), this_y(end), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
end

colormap(turbo);

% Scale bar definition
scale_len = 1.0; 
x_limit_max = (num_instances-1)*spacing + 0.6;
x_limit_min = -0.6;
y_limit_min = -0.8;
y_limit_max = 0.8;

scale_x_end = x_limit_max - 0.1;
scale_x_start = scale_x_end - scale_len;
scale_y = y_limit_min + 0.2;

line([scale_x_start, scale_x_end], [scale_y, scale_y], 'Color', 'k', 'LineWidth', 4);
text((scale_x_start + scale_x_end)/2, scale_y + 0.18, '1 \mum', ...
     'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', ...
     'FontSize', 14, 'FontWeight', 'bold');

axis equal;
xlim([x_limit_min, x_limit_max]);
ylim([y_limit_min - 0.2, y_limit_max + 0.2]);

cb = colorbar;
cb.Label.String = 'Time (s)';
cb.FontName = 'Times New Roman';
cb.FontSize = 11;

hold off;