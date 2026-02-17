clc; clear; close all;

% Parameters
N = 500;                  
Lx = 20; Ly = 20;         
dt = 0.003;               
TMAX = 10;                
steps = TMAX / dt;

v0 = 1;                   
Dr = 0.1;                 
alpha_m = 0.5;              
chemo_m = 6;              
kii = 1;                 
lambda = 2;
dcutoff = 0.6;

gamma1 = 10;

x_source = 10;
y_source = 10;
R_source = 1;

theta_circ = linspace(0, 2*pi, 100);
x_circ = x_source + R_source * cos(theta_circ);
y_circ = y_source + R_source * sin(theta_circ);

x = [];
y = [];
while numel(x) < N
    xt = rand * Lx;
    yt = rand * Ly;
    if hypot(xt - x_source, yt - y_source) > R_source
        x(end+1) = xt; 
        y(end+1) = yt; 
    end
end

figure;

% vidObj = VideoWriter('cell_wall_simulation.avi');
% vidObj.FrameRate = 10;
% open(vidObj);

for t = 1:steps
    theta = 2 * pi * rand(1, N);

    dis = zeros(1, N);
    for i = 1:N
        dxs = x_source - x(i);
        dys = y_source - y(i);
        dis(i) = hypot(dxs, dys);
        if dis(i) < chemo_m
            angle_to_source = atan2(dys, dxs);
            dtheta = (alpha_m * (chemo_m - dis(i)) / chemo_m) ...
                     * sin(angle_to_source - theta(i)) ...
                     + sqrt(2 * Dr * dt) * randn;
            theta(i) = theta(i) + dtheta;
        end
    end

    force_x = zeros(1, N);
    force_y = zeros(1, N);
    for i = 1:N
        xi = x(i); yi = y(i);
        for j = 1:N
            if i == j, continue; end
            dx_ij = xi - x(j);
            dx_ij = dx_ij - Lx * round(dx_ij / Lx);
            dy_ij = yi - y(j);
            dy_ij = dy_ij - Ly * round(dy_ij / Ly);
            r = hypot(dx_ij, dy_ij);
            if r < dcutoff && r > 1e-6
                f = (r^lambda) + (dcutoff/r)^lambda - 0.5 - dcutoff^lambda;
                force_x(i) = force_x(i) + kii * f * dx_ij / r;
                force_y(i) = force_y(i) + kii * f * dy_ij / r;
            end
        end
    end

    x_old = x; y_old = y;

    x_next = x + gamma1 * (v0*cos(theta) + force_x) * dt;
    y_next = y + gamma1 * (v0*sin(theta) + force_y) * dt;

    r_to_source = hypot(x_next - x_source, y_next - y_source);
    in_source = r_to_source < R_source;

    valid = ~in_source;
    x(valid) = x_next(valid);
    y(valid) = y_next(valid);

    % boundary conditions
    x(x < 0)   = -x(x < 0);
    x(x > Lx)  = 2*Lx - x(x > Lx);
    y(y < 0)   = -y(y < 0);
    y(y > Ly)  = 2*Ly - y(y > Ly);

    disp_mag = hypot(x - x_old, y - y_old);

    displacement_threshold = 1.0e-4;
    num_active = sum(disp_mag > displacement_threshold);

    if t > 10
        max_cells = 600; 
        if num_active < 200 && N < max_cells
            num_to_add = min(200 - num_active, max_cells - N);
            for k = 1:num_to_add
                added = false;
                while ~added
                    xt = rand * Lx;
                    yt = rand * Ly;
                    if hypot(xt - x_source, yt - y_source) > R_source
                        x(end+1) = xt; 
                        y(end+1) = yt; 
                        added = true;
                    end
                end
            end
            N = numel(x);  
        end
    end

    % draw
    if mod(t, 50) == 0
        clf;
        plot(x, y, 'bo', 'MarkerFaceColor', 'b'); hold on;
        fill(x_circ, y_circ, [190 30 45]/255, 'EdgeColor', 'none');  % source
        axis([0 Lx 0 Ly]); axis manual;
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
        title(sprintf('Time = %.2f', t * dt));
        drawnow;
        % frame = getframe(gcf);
        % writeVideo(vidObj, frame);
    end
end

% close(vidObj);
