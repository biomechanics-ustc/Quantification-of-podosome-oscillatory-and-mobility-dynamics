%% === Initialization & Physical Simulation ===
clear; clc; close all;

R = 0.8; c = 100; k_spring = 0.85;
asymptote_um = 60.0/1000; 
tau_corr = 20; num_springs = 150; ron = 0.15; roff = 0.1;
F_brown_mag = 2.2; dt = 1.0; T = 600; t = 0:dt:T; N_steps = length(t);

theta = linspace(0, 2*pi, num_springs+1); theta = theta(1:end-1);
anchor_points_base = [R * cos(theta); R * sin(theta)];
xr0 = -0.25 + 0.2 * sin(2*pi*(1/300.0)*t);

pos_all = zeros(2, N_steps); 
clutch_state_all = false(num_springs, N_steps);
force_mag_all = zeros(num_springs, N_steps);

current_pos = [0; 0]; 
pos_all(:, 1) = current_pos;
spring_connected = rand(1, num_springs) < (ron/(ron+roff));
current_brown_force = [0; 0];
exp_factor = exp(-dt/tau_corr);
noise_prefactor = F_brown_mag * sqrt(1 - exp_factor^2);

fprintf('Running physical simulation...\n');

for i = 1:N_steps-1        
    % O-U Brownian force update
    white_noise = randn(2, 1);
    current_brown_force = current_brown_force * exp_factor + noise_prefactor * white_noise;
    
    % Markov transitions
    spring_connected = (rand(1, num_springs) < (1-exp(-ron*dt)) & ~spring_connected) | ...
                       (~(rand(1, num_springs) < (1-exp(-roff*dt))) & spring_connected);
    
    clutch_state_all(:, i) = spring_connected; 
    
    % Elastic force evaluation
    F_spring_net = [0; 0];
    step_forces = zeros(num_springs, 1);
    
    idx = find(spring_connected);
    if ~isempty(idx)
        r_vecs = current_pos - anchor_points_base(:,idx);
        r_mags = sqrt(sum(r_vecs.^2, 1));
        valid = r_mags > 0;
        if any(valid)
            valid_idx = idx(valid);
            f_scalar = -k_spring * (r_mags(valid) - (R + xr0(i)));
            step_forces(valid_idx) = abs(f_scalar);
            
            forces = f_scalar .* (r_vecs(:,valid) ./ r_mags(valid));
            F_spring_net = sum(forces, 2);
        end
    end
    force_mag_all(:, i) = step_forces; 
    
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
    pos_all(:, i+1) = current_pos;
end

clutch_state_all(:, N_steps) = clutch_state_all(:, N_steps-1);
force_mag_all(:, N_steps) = force_mag_all(:, N_steps-1);

fprintf('Physical simulation completed.\n');

%% === MSD & Diffusion Exponent Calculation ===
fprintf('Calculating MSD and alpha...\n');

% Time-averaged MSD
max_lag = floor(N_steps / 4);
msd = zeros(max_lag, 1);
lags = (1:max_lag)' * dt;

for lag = 1:max_lag
    dx = pos_all(1, 1+lag:end) - pos_all(1, 1:end-lag);
    dy = pos_all(2, 1+lag:end) - pos_all(2, 1:end-lag);
    sq_disp = dx.^2 + dy.^2;
    msd(lag) = mean(sq_disp);
end

% Power-law fitting
fit_idx = 1:floor(max_lag/2); 
log_lags = log10(lags(fit_idx));
log_msd = log10(msd(fit_idx));

p = polyfit(log_lags, log_msd, 1);
alpha = p(1);

fprintf('Calculation completed. Alpha = %.4f\n', alpha);
fprintf('Starting video rendering...\n');

%% === Video Rendering & Export ===
save_path = 'C:\Users\jkdhryhd\Desktop\reserch\Podosome gg\PPT'; 

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

video_filename = fullfile(save_path, 'Podosome_Clutch_Dynamics.mp4');

vidObj = VideoWriter(video_filename, 'MPEG-4');
vidObj.FrameRate = 30; 
open(vidObj);

fig = figure('Name', 'Podosome_Clutch_Video', 'Position', [100, 100, 850, 750], 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');

% Global bounds
max_force = max(force_mag_all(:));
if max_force < 1e-6; max_force = 1e-6; end 

lim_bound = R + max(max(abs(pos_all))) + 0.1; 
axis(ax, 'equal');
xlim(ax, [-lim_bound, lim_bound]);
ylim(ax, [-lim_bound, lim_bound]);
set(ax, 'LineWidth', 2.2, 'Box', 'on', 'XTick', [], 'YTick', [], 'XColor', 'k', 'YColor', 'k');

cmap = colormap(ax, turbo);
num_colors = size(cmap, 1);
caxis(ax, [0, max_force]);
cb = colorbar(ax);
cb.Label.String = 'Clutch Force Magnitude (pN)';
cb.FontName = 'Times New Roman';
cb.FontSize = 13;
cb.FontWeight = 'bold';

visual_scale = 0.6; 
visual_R = R * visual_scale;
visual_anchor_points = anchor_points_base * visual_scale;

theta_circle = linspace(0, 2*pi, 100);
plot(ax, visual_R*cos(theta_circle), visual_R*sin(theta_circle), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);

scatter(ax, visual_anchor_points(1,:), visual_anchor_points(2,:), 25, [0.5 0.5 0.5], 'filled');

h_springs = gobjects(num_springs, 1);
for s = 1:num_springs
    h_springs(s) = plot(ax, [NaN, NaN], [NaN, NaN], 'LineWidth', 2.5);
end

h_traj = plot(ax, NaN, NaN, '-', 'Color', [0, 0, 0, 0.6], 'LineWidth', 3);

h_core = scatter(ax, pos_all(1,1), pos_all(2,1), 500, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

h_title = title(ax, 'Podosome Dynamics: t = 0.0 s', 'FontSize', 18, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

text(ax, -lim_bound * 0.9, lim_bound * 0.9, sprintf('\\alpha = %.3f', alpha), ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
    'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

% Z-order enforcement
uistack(h_core, 'top');
uistack(h_traj, 'top');
uistack(h_core, 'top'); 

% Frame-by-frame rendering
for i = 1:N_steps
    h_title.String = sprintf('Podosome Dynamics: t = %.1f s', t(i));
    
    h_traj.XData = pos_all(1, 1:i);
    h_traj.YData = pos_all(2, 1:i);
    
    h_core.XData = pos_all(1,i);
    h_core.YData = pos_all(2,i);
    
    for s = 1:num_springs
        if clutch_state_all(s, i)
            f_val = force_mag_all(s, i);
            c_idx = round((f_val / max_force) * (num_colors - 1)) + 1;
            c_idx = max(1, min(num_colors, c_idx)); 
            
            h_springs(s).XData = [visual_anchor_points(1, s), pos_all(1,i)];
            h_springs(s).YData = [visual_anchor_points(2, s), pos_all(2,i)];
            h_springs(s).Color = cmap(c_idx, :);
            h_springs(s).Visible = 'on';
        else
            h_springs(s).Visible = 'off'; 
        end
    end
    
    drawnow limitrate;
    frame = getframe(fig);
    writeVideo(vidObj, frame);
    
    if mod(i, 60) == 0
        fprintf('Rendering frame: %d / %d...\n', i, N_steps);
    end
end

close(vidObj);
fprintf('Video saved successfully to: %s\n', video_filename);