clc; clear; close all;

% Parameters
N = 500;
Lx = 20; Ly = 20;
dt = 0.003;
TMAX = 10;
steps = TMAX / dt;

n = 250;            
add_every = 1500;   % add new cells 
max_cells = 1000;   % optional cap on total cells

ssteps = round(TMAX / dt);
region_counts = zeros(1, ssteps);

v0 = 1;
Dr = 0.1;           
alpha_m = 0;        
kii = 1;
lambda = 2;
dcutoff = 0.6;

gamma1 = 10;

% Source 
x_source = 10;
y_source = 10;
R_source = 1;

theta_circ = linspace(0, 2*pi, 100);
x_circ = x_source + R_source * cos(theta_circ);
y_circ = y_source + R_source * sin(theta_circ);


% Initialize
x = rand(1, N) * Lx;
y = rand(1, N) * Ly;

figure;

% vidObj = VideoWriter('cell_wall_simulation.avi');
% vidObj.FrameRate = 5;
% open(vidObj);

for t = 1:steps
    theta = 2 * pi * rand(1, N);

    force_x = zeros(1, N);
    force_y = zeros(1, N);
    for i = 1:N
        for j = 1:N
            if i == j, continue; end
            dx_ij = x(i) - x(j);
            dx_ij = dx_ij - Lx * round(dx_ij / Lx);  
            dy_ij = y(i) - y(j);
            dy_ij = dy_ij - Ly * round(dy_ij / Ly);  
            r = sqrt(dx_ij^2 + dy_ij^2);
            if r < dcutoff && r > 1e-6
                f = (r^lambda) + (dcutoff/r)^lambda - 0.5 - dcutoff^lambda;
                fx = kii * f * dx_ij / r;
                fy = kii * f * dy_ij / r;
                force_x(i) = force_x(i) + fx;
                force_y(i) = force_y(i) + fy;
            end
        end
    end

    x_next = x + gamma1 * (v0*cos(theta) + force_x) * dt;
    y_next = y + gamma1 * (v0*sin(theta) + force_y) * dt;

r_to_source = sqrt((x_next - x_source).^2 + (y_next - y_source).^2);
in_source = r_to_source < R_source;

valid = ~in_source;
x(valid) = x_next(valid);
y(valid) = y_next(valid);

    x = x_next;
    y = y_next;

    % boundary conditions
    x(x < 0)   = -x(x < 0);
    x(x > Lx)  = 2*Lx - x(x > Lx);
    y(y < 0)   = -y(y < 0);
    y(y > Ly)  = 2*Ly - y(y > Ly);

    in_region = (x > 7) & (x < 13) & (y > 7) & (y < 13);
    region_counts(t) = sum(in_region);

    if mod(t, add_every) == 0 && N < max_cells
        new_x = [];
        new_y = [];
        tries = 0; 
        target_new = min(max_cells - N, round(sum(in_region) * 1)); 
        while numel(new_x) < target_new && tries < 1000
            xt = 8 + (12 - 8) * rand;   
            yt = 6 + (14 - 6) * rand;   
            new_x(end+1) = xt; 
            new_y(end+1) = yt; 
            tries = tries + 1;
        end

        if ~isempty(new_x)
            x = [x new_x];
            y = [y new_y];
            theta = [theta 2*pi*rand(1, numel(new_x))]; 
            N = numel(x);
        end
    end

    
    % draw
    if mod(t, 50) == 0
        clf;
        plot(x, y, 'bo', 'MarkerFaceColor', 'b'); hold on;
        fill(x_circ, y_circ, [190 30 45]/255, 'EdgeColor', 'none');
        axis([0 Lx 0 Ly]); axis manual;
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
        title(sprintf('Time = %.2f', t * dt));
        drawnow;
        % frame = getframe(gcf);
        % writeVideo(vidObj, frame);
    end

    drawnow;
end

% close(vidObj);
