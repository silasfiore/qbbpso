clear all; close all;
addpath("utils/");

%initial attitude quaternion
q0 = randquat();
%magnitude of initial angular velocity
wb0_magnitude = rand();
%initial angular velocity (body frame)
wb0 = randuvec().*wb0_magnitude;

% for object with axial symmetry
It = rand()*500;
Ia = rand()*500;

% inertia ratios, symmetric case
c = [(It-Ia)/It, (Ia - It)/It, 0];
c = circshift(c,randi(3));

% general case (much slower!)
%I1 = rand()*500;
%I2 = rand()*500;
%I3 = rand()*500;
%I_sorted = sort([I1,I2,I3]);

% inertia ratios, general case
%c = [(I2 - I3)/I1,(I3 - I1)/I2,(I1-I2)/I3];



%inertia constants of body
c1 = c(1);
c2 = c(2);
c3 = c(3);

if c1 == 0
    %initial normalized angular momentum (body frame)
    lb0 = diag([(1-c2),1,1])*wb0;
    %normalized angular momentum
    l = quat2rotmat(q0)*lb0;
    %normalized angular momentum axis
    l0 = l./norm(l);
elseif c2 == 0
    %initial normalized angular momentum (body frame)
    lb0 = diag([1,(1-c3),1])*wb0;
    %normalized angular momentum
    l = quat2rotmat(q0)*lb0;
elseif c3 == 0
    %initial normalized angular momentum (body frame)
    lb0 = diag([1,1,(1-c1)])*wb0;
    %normalized angular momentum
    l = quat2rotmat(q0)*lb0;
else
    %initial normalized angular momentum (body frame)
    lb0 = diag([I1,I2,I3])*wb0./I_sorted(2);
    %normalized angular momentum
    l = quat2rotmat(q0)*lb0;
end

y0 = [wb0',q0'];

% simulate for one precession cycle
tspan = [0 2*pi/norm(l)];
options = odeset('RelTol',eps(1)*100,'AbsTol',eps(1)*100);
% numerical integration with renormalization
[t,y] = ode45(@(t,y)odefuneuler(t,y,c1,c2,c3),tspan,y0,options);

input_file = 'inputdata.txt';
output_file = 'particles.txt';

% number of measurements
n_measurements = 10;

% sample at random timestamps
isample = sort(randperm(length(t),n_measurements));
% sample at evenly spaced timestamps
%isample = round(linspace(1,length(t),n_measurements));


fileID = fopen(input_file,'w');
for k=1:n_measurements
    i = isample(k);
        
    trel = t(i) - t(1);

    % random vector in the body frame is observed
    b = randuvec();


    r = quat2rotmat(y(i,4:7)) * b;

    fprintf(fileID,'%f %f %f %f %f %f %f\n',[trel,r',b']);
end
fclose(fileID);


n_particles = 10000;
n_iterations = 100;

mag_l_prior = norm(l);
std_l_prior = 0.1*norm(l);
timestep_numint = 0.001; % only used in the general (numerically integrated) case

cmd = sprintf('../qbbpso %s %f %f %f %f %f -dt %f -n %d -i %d -o %s',...
    input_file,mag_l_prior,std_l_prior,c1,c2,c3,timestep_numint, n_particles, n_iterations,output_file);
tic;
% run qbbpso
system(cmd);
elapsedTime = toc;
fprintf('The system call took %.4f seconds.\n', elapsedTime);


% plotting
parts = readmatrix(output_file);
parts = sortrows(parts,1,'ascend');

mrps = zeros(3,length(parts));
for i=1:length(parts)
    mrps(:,i) = quat2mrp(parts(i,2:5));
end
mrp_gt = quat2mrp(q0);
figure; hold on; grid on;
view(3)
[X, Y, Z] = sphere(30);
srf = surf(X, Y, Z);
set(srf, 'FaceColor', 'none', ...
       'EdgeColor', [0.3 0.3 0.3], ...
       'LineWidth', 0.8);
scatter3(mrps(1,:), mrps(2,:), mrps(3,:), [],parts(:,1),'.');
plot3([-1,1], [mrp_gt(2),mrp_gt(2)], [mrp_gt(3), mrp_gt(3)], 'r');
plot3( [mrp_gt(1),mrp_gt(1)],[-1,1], [mrp_gt(3), mrp_gt(3)], 'r');
plot3( [mrp_gt(1),mrp_gt(1)],[mrp_gt(2),mrp_gt(2)],[-1,1], 'r');
axis equal
daspect([1 1 1]);
pbaspect([1 1 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Modified Rodriguez parameter of initial attitude quaternion');
clim([min(parts(:,1)), max(parts(:,1))]);
cb = colorbar;
cb.Label.String = 'attitude rms (°)';

figure; hold on; grid on;
view(240,70)
scatter3(parts(:,6).*parts(:,9), parts(:,7).*parts(:,9), parts(:,8).*parts(:,9), [], parts(:,1),'.');
plot3([0,l(1)],[0,l(2)],[0,l(3)],'r');
axis equal
daspect([1 1 1]);
pbaspect([1 1 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Inertial angular momentum normalized by intermediate MOI');
clim([min(parts(:,1)), max(parts(:,1))]);
colorbar;
cb = colorbar;
cb.Label.String = 'attitude rms (°)';


