clc; clear; close all;

% Parameters
N = 300;                  % number of cells
Lx = 20; Ly = 20;         % domain size
dt = 0.003;                % time step
TMAX = 10;               % total simulation time
steps = TMAX / dt;

v0 = 1; %0.5;                 % self-propulsion speed
Dr = 0.1;                 % rotational diffusion
alpha_m = 1;            % chemotactic strength
chemo_m = 7; %7
kii = 1;                  % cell-cell repulsion
kiw = 5;                  % cell-wall repulsion
lambda = 2;
dcutoff = 0.6;

gamma1 = 10;
gamma2 = 1;

kbend = 1;     % bending stiffness
kpress = 0.1;     % internal pressure
kspring = 0.01;  % spring constant

% Source 
x_source = 10;
y_source = 10;
R_source = 1;

% Wall 
xw_center = 18;
yw_center = 10;
R_wall = 6;
Nw = 100;                              
theta_wall = linspace(0, 2*pi, Nw+1); 
theta_wall(end) = [];

xw = xw_center + R_wall * cos(theta_wall);
yw = yw_center + R_wall * sin(theta_wall);

xw2_center = 2;
yw2_center = 10;
R_wall2 = 6;
Nw2 = 100;
theta_wall2 = linspace(0, 2*pi, Nw2+1); theta_wall2(end) = [];
% theta_wall2 = fliplr(theta_wall2); % Ensure counter-clockwise like wall1
xw2 = xw2_center + R_wall2 * cos(theta_wall2);
yw2 = yw2_center + R_wall2 * sin(theta_wall2);

theta_circ = linspace(0, 2*pi, 100);
x_circ = x_source + R_source * cos(theta_circ);
y_circ = y_source + R_source * sin(theta_circ);



% Initialize
x = [];
y = [];
while length(x) < N
    xt = rand * Lx;
    yt = rand * Ly;
    if norm([xt - x_source, yt - y_source]) > R_source
        x(end+1) = xt;
        y(end+1) = yt;
    end
end

figure;

vidObj = VideoWriter('cell_wall_simulation.avi');
vidObj.FrameRate = 10;
open(vidObj);

% theta = 2 * pi * rand(1, N);
for t = 1:steps


theta = 2 * pi * rand(1, N);
    dis = zeros(1, N);
    for i = 1:N
        dx = x_source - x(i);
        dy = y_source - y(i);
        dis(i) = sqrt(dx^2 + dy^2);
        angle_to_source = atan2(dy, dx);
        if dis(i) < chemo_m
dtheta = (alpha_m * (chemo_m - dis(i)) / chemo_m) * sin(angle_to_source - theta(i)) + sqrt(2 * Dr * dt) * randn;
            theta(i) = theta(i) + dtheta;
        end
    end

    % repulsion
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

fw_x = zeros(1, Nw);
fw_y = zeros(1, Nw);

for i = 1:N
    for j = 1:Nw
        dxw = x(i) - xw(j);
        dyw = y(i) - yw(j);
        r = sqrt(dxw^2 + dyw^2);
        if r < dcutoff && r > 1e-6
f = kiw * (1 - r/dcutoff)^2;
            fx = kiw * f * dxw / r;
            fy = kiw * f * dyw / r;

            force_x(i) = force_x(i) + fx;
            force_y(i) = force_y(i) + fy;

            fw_x(j) = fw_x(j) - fx;
            fw_y(j) = fw_y(j) - fy;
        end
    end
end

fw2_x = zeros(1, Nw2);
fw2_y = zeros(1, Nw2);

for i = 1:N
    for j = 1:Nw2
        dxw = x(i) - xw2(j);
        dyw = y(i) - yw2(j);
        r = sqrt(dxw^2 + dyw^2);
        if r < dcutoff && r > 1e-6
            f = kiw * (1 - r/dcutoff)^2;
            fx = kiw * f * dxw / r;
            fy = kiw * f * dyw / r;

            force_x(i) = force_x(i) + fx;
            force_y(i) = force_y(i) + fy;

            fw2_x(j) = fw2_x(j) - fx;
            fw2_y(j) = fw2_y(j) - fy;
        end
    end
end

    x_next = x + gamma1 * (v0*cos(theta) + force_x) * dt;
    y_next = y + gamma1 * (v0*sin(theta) + force_y) * dt;
    

r_to_source = sqrt((x_next - x_source).^2 + (y_next - y_source).^2);
in_source = r_to_source < R_source;

in_wall = inpolygon(x_next, y_next, xw, yw);
in_wall2 = inpolygon(x_next, y_next, xw2, yw2);

valid = ~in_source & ~in_wall & ~in_wall2;
x(valid) = x_next(valid);
y(valid) = y_next(valid);

min_dist_to_wall = 0.3;
keep = true(1, N);
for i = 1:N
    dmin = min(sqrt((x(i) - xw).^2 + (y(i) - yw).^2));
    if dmin < min_dist_to_wall
        keep(i) = false;
    end
end
x = x(keep);
y = y(keep);
theta = theta(keep);
N = length(x);  

dx_disp = x_next(keep) - x;
dy_disp = y_next(keep) - y;
disp_mag = sqrt(dx_disp.^2 + dy_disp.^2);



displacement_threshold = 0.0001;  
active_idx = disp_mag > displacement_threshold;
num_active = sum(active_idx);


% BC:
x(x < 0)   = -x(x < 0);
x(x > Lx)  = 2*Lx - x(x > Lx);

y(y < 0)   = -y(y < 0);
y(y > Ly)  = 2*Ly - y(y > Ly);


    for i = 1:Nw
        ip = mod(i, Nw) + 1;
        im = mod(i - 2, Nw) + 1;

        % Bending
        vx1 = xw(i) - xw(im); vy1 = yw(i) - yw(im);
        vx2 = xw(ip) - xw(i); vy2 = yw(ip) - yw(i);
        len1 = sqrt(vx1^2 + vy1^2); len2 = sqrt(vx2^2 + vy2^2);
        if len1 > 1e-6 && len2 > 1e-6
            cos_theta = (vx1*vx2 + vy1*vy2) / (len1*len2);
            sin_theta = vx1*vy2 - vy1*vx2;
            torque = -kbend * sin_theta;
            fx = torque * (vy1/len1 + vy2/len2);
            fy = -torque * (vx1/len1 + vx2/len2);
            fw_x(i) = fw_x(i) + fx;
            fw_y(i) = fw_y(i) + fy;
        end

        % Pressure 
        dx1 = xw(ip) - xw(i);
        dy1 = yw(ip) - yw(i);
        edge_length = sqrt(dx1^2 + dy1^2);
        nx = dy1 / edge_length;
        ny = -dx1 / edge_length;
        fw_x(i) = fw_x(i) + kpress * nx;
        fw_y(i) = fw_y(i) + kpress * ny;

        % Spring force 
rest_length = 2 * R_wall * sin(pi / Nw);  
dx_s = xw(ip) - xw(i);
dy_s = yw(ip) - yw(i);
dist = sqrt(dx_s^2 + dy_s^2);
Fspring = -kspring * (dist - rest_length);
fx_s = Fspring * dx_s / dist;
fy_s = Fspring * dy_s / dist;

fw_x(i)  = fw_x(i)  + fx_s;
fw_y(i)  = fw_y(i)  + fy_s;
fw_x(ip) = fw_x(ip) - fx_s;
fw_y(ip) = fw_y(ip) - fy_s;
    end

    xw = xw + gamma2 * fw_x * dt;
    yw = yw + gamma2 * fw_y * dt;

    for i = 1:Nw2
    ip = mod(i, Nw2) + 1;
    im = mod(i - 2, Nw2) + 1;

    % Bending
    vx1 = xw2(i) - xw2(im); vy1 = yw2(i) - yw2(im);
    vx2 = xw2(ip) - xw2(i); vy2 = yw2(ip) - yw2(i);
    len1 = sqrt(vx1^2 + vy1^2); len2 = sqrt(vx2^2 + vy2^2);
    if len1 > 1e-6 && len2 > 1e-6
        cos_theta = (vx1*vx2 + vy1*vy2) / (len1*len2);
        sin_theta = vx1*vy2 - vy1*vx2;
        torque = -kbend * sin_theta;
        fx = torque * (vy1/len1 + vy2/len2);
        fy = -torque * (vx1/len1 + vx2/len2);
        fw2_x(i) = fw2_x(i) + fx;
        fw2_y(i) = fw2_y(i) + fy;
    end

    % Pressure
    dx1 = xw2(ip) - xw2(i);
    dy1 = yw2(ip) - yw2(i);
    edge_length = sqrt(dx1^2 + dy1^2);
    nx = dy1 / edge_length;
    ny = -dx1 / edge_length;
    fw2_x(i) = fw2_x(i) + kpress * nx;
    fw2_y(i) = fw2_y(i) + kpress * ny;

    % Spring
    rest_length = 2 * R_wall2 * sin(pi / Nw2);
    dx_s = xw2(ip) - xw2(i);
    dy_s = yw2(ip) - yw2(i);
    dist = sqrt(dx_s^2 + dy_s^2);
    Fspring = -kspring * (dist - rest_length);
    fx_s = Fspring * dx_s / dist;
    fy_s = Fspring * dy_s / dist;
    fw2_x(i)  = fw2_x(i)  + fx_s;
    fw2_y(i)  = fw2_y(i)  + fy_s;
    fw2_x(ip) = fw2_x(ip) - fx_s;
    fw2_y(ip) = fw2_y(ip) - fy_s;
end

xw2 = xw2 + gamma2 * fw2_x * dt;
yw2 = yw2 + gamma2 * fw2_y * dt;


if t > 10
max_cells = 600;
if num_active < 200 && N < max_cells
    num_to_add = min(200 - num_active, max_cells - N);
    for k = 1:num_to_add
        added = false;
        while ~added
            xt = rand * Lx;
            yt = rand * Ly;
            if norm([xt - x_source, yt - y_source]) > R_source && ~inpolygon(xt, yt, xw, yw)
                x(end+1) = xt;
                y(end+1) = yt;
                theta(end+1) = 2*pi*rand;
                added = true;
            end
        end
    end
    N = length(x);  % update cell count
end
end

    if mod(t, 50) == 0
        clf;
        plot(x, y, 'bo', 'MarkerFaceColor', 'b'); hold on;
        fill(xw, yw, [147 149 152]/255, 'EdgeColor', 'k', 'LineWidth', 2); % wall
        fill(xw2, yw2, [147 149 152]/255, 'EdgeColor', 'k', 'LineWidth', 2); % wall


fill(x_circ, y_circ, [190 30 45]/255, 'EdgeColor', 'none');
        % viscircles([x_source, y_source], R_source, 'Color', 'k', 'LineWidth', 1);
        axis([0 Lx 0 Ly]); axis manual;
        set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
        title(sprintf('Time = %.2f', t * dt));
        
        drawnow;
        frame = getframe(gcf);
        writeVideo(vidObj, frame);
    end
end

close(vidObj);