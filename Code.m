clc, clear
data = readmatrix('cy1.xlsx');

% Extract relevant data
LTHI = data(:, 13:15);
LKNE = data(:, 16:18);
LTIB = data(:, 19:21);
LANK = data(:, 22:24);
LHEE_Z = data(:, 26);

% Calculate thigh and shank vectors and angles
thigh_vec = LKNE - LTHI;
shank_vec = LANK - LTIB;

angle_thigh = atan2(thigh_vec(:,2), thigh_vec(:,1));
angle_shank = atan2(shank_vec(:,2), shank_vec(:,1));

%% DRP (Discrete Relative Phase)
[~, peak_thigh] = findpeaks(angle_thigh);
[~, peak_shank] = findpeaks(angle_shank);

T1 = length(angle_thigh); 

if ~isempty(peak_thigh) && ~isempty(peak_shank)
    dt = abs(peak_thigh(1) - peak_shank(1));
    drp = dt / T1 * 360;
else
    drp = NaN;
end

%% CRP (Continuous Relative Phase)
dt = 1;  
vel_thigh = gradient(angle_thigh, dt);
vel_shank = gradient(angle_shank, dt);

% Normalize angles and velocities
thigh_n = 2 * (angle_thigh - min(angle_thigh)) / (max(angle_thigh) - min(angle_thigh)) - 1;
shank_n = 2 * (angle_shank - min(angle_shank)) / (max(angle_shank) - min(angle_shank)) - 1;
vel_thigh_n = vel_thigh / max(abs(vel_thigh));
vel_shank_n = vel_shank / max(abs(vel_shank));

% Calculate phase angles
phase_thigh = atan2(vel_thigh_n, thigh_n);
phase_shank = atan2(vel_shank_n, shank_n);
crp = rad2deg(wrapToPi(phase_thigh - phase_shank)); 

% Adjust CRP to be within [-180, 180]
crp(crp > 180) = crp(crp > 180) - 360;
crp(crp < -180) = crp(crp < -180) + 360;
crp = abs(crp);  

% Calculate cycle percentage and DRP over a specific range
cycle_percent = linspace(0, 100, length(angle_thigh));
time = linspace(0, 1, length(angle_thigh)); % Assuming 1-second duration for plotting
drp_range = cycle_percent >= 30 & cycle_percent <= 60;
drp = mean(crp(drp_range));

%% Plotting
figure;

% First column: Original time series (angles and velocities)
% Thigh angle
subplot(4, 4, 1)
plot(time, rad2deg(angle_thigh), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta_1 (deg)');
grid on;

% Thigh velocity
subplot(4, 4, 5)
plot(time, rad2deg(vel_thigh), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\omega_1 (deg/s)');
grid on;

% Shank angle
subplot(4, 4, 9)
plot(time, rad2deg(angle_shank), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\theta_2 (deg)');
grid on;

% Shank velocity
subplot(4, 4, 13)
plot(time, rad2deg(vel_shank), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\omega_2 (deg/s)');
grid on;

% Second column: Normalized phase planes
% Thigh phase plane
subplot(4, 4, 2)
plot(thigh_n, vel_thigh_n, 'b', 'LineWidth', 1.5);
xlabel('\theta_1_n');
ylabel('\omega_1_n');
grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);

% Shank phase plane
subplot(4, 4, 10)
plot(shank_n, vel_shank_n, 'r', 'LineWidth', 1.5);
xlabel('\theta_2_n');
ylabel('\omega_2_n');
grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);

% Third column: Phase angles over time
% Thigh phase angle
subplot(4, 4, 3)
plot(time, rad2deg(phase_thigh), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\phi_1 (deg)');
grid on;
ylim([-180 180]);

% Shank phase angle
subplot(4, 4, 11)
plot(time, rad2deg(phase_shank), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('\phi_2 (deg)');
grid on;
ylim([-180 180]);

% Fourth column: CRP
subplot(4, 4, 4)
plot(time, crp, 'g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('CRP (deg)');
grid on;
ylim([0 180]);

% Adjust figure size
set(gcf, 'Position', [100 100 1200 800]);
%% Vector Coding
gamma = zeros(length(angle_thigh)-1, 1);
for i = 1:length(gamma)
    dy = angle_shank(i+1) - angle_shank(i);
    dx = angle_thigh(i+1) - angle_thigh(i);
    gamma(i) = atan2d(dy, dx);
    if gamma(i) < 0
        gamma(i) = gamma(i) + 360;
    end
end

figure;
subplot(3,1,1)
plot(angle_thigh); hold on; plot(angle_shank); title('Segment Angles');
legend('Thigh','Shank')

subplot(3,1,2)
plot(crp); title('Continuous Relative Phase (CRP)');

subplot(3,1,3)
plot(gamma); title('Coupling Angle (Vector Coding)');

%% plots

bins = 0:45:360;
labels = {'0-45 (In-phase)', '45-90 (Thigh dominant)', '90-135 (Thigh leading)', ...
          '135-180 (Anti-phase)', '180-225 (Anti-phase)', '225-270 (Shank leading)', ...
          '270-315 (Shank dominant)', '315-360 (In-phase)'};

[counts, ~] = histcounts(gamma, bins);
frequencies = counts / sum(counts) * 100; 

figure;
polarhistogram(deg2rad(gamma), deg2rad(bins), 'Normalization', 'probability', 'FaceColor', 'b', 'EdgeColor', 'k');
title('Circular Plot of Coupling Angles (Vector Coding)', 'FontSize', 12);
set(gca, 'ThetaTick', 0:45:315, 'ThetaTickLabel', {'0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°'});


disp('Frq');
for i = 1:length(labels)
    fprintf('%s: %.2f%%\n', labels{i}, frequencies(i));
end


in_phase = sum(frequencies([1, end])); % 0-45 و 315-360
anti_phase = sum(frequencies([4, 5])); % 135-180 و 180-225
thigh_dominant = frequencies(2); % 45-90
shank_dominant = frequencies(7); % 270-315
thigh_leading = frequencies(3); % 90-135
shank_leading = frequencies(6); % 225-270

disp('coordintion:');
fprintf('In-phase : %.2f%%\n', in_phase);
fprintf('Anti-phase : %.2f%%\n', anti_phase);
fprintf('Thigh dominant : %.2f%%\n', thigh_dominant);
fprintf('Shank dominant : %.2f%%\n', shank_dominant);
fprintf('Thigh leading : %.2f%%\n', thigh_leading);
fprintf('Shank leading : %.2f%%\n', shank_leading);