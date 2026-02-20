% out = process_sonic_log_and_z0('C:\Users\ISMER\OneDrive\Desktop\Estrie 13 septembre\anemometer\anemometer_2025-09-13_08-55-15\anemometer_log.csv');

% Whole file (single block)
% out = process_sonic_log_blocks('C:\Users\ISMER\OneDrive\Desktop\Estrie\Data\anemometer_2025-09-13_08-45-12\anemometer_log.csv', Inf, 3.5);

% % 1-minute blocks
% out = process_sonic_log_blocks('C:\Users\ISMER\OneDrive\Desktop\Estrie 13 septembre\anemometer\anemometer_2025-09-13_08-45-12\anemometer_log.csv', 1, 3.5);
% 
% % 5-minute blocks
% clc
% close all
% clear all
out = process_sonic_log_blocks('C:\Users\ISMER\OneDrive\Desktop\Scripts\Estrie_Acquisition\anemometer_2025-09-21_11-01-10\anemometer_log.csv', 2, 3.5);


%% Compute and plot mean wind speed and wind direction from out = process_sonic_log_blocks(...)

% Assume  already did:
% out = process_sonic_log_blocks('file.csv', 5, 3.5);

tt      = out.tt;          % high-frequency timetable
blocks  = out.blocks;      % block summary table
nBlocks = height(blocks);
N       = height(tt);

% --- Recover the original block size in number of samples ---
Nwin = floor(N / nBlocks);

% --- Preallocate arrays for block-mean U and V ---
Ubar = nan(nBlocks,1);
Vbar = nan(nBlocks,1);

% --- Compute block-mean U and V from the raw data ---
for b = 1:nBlocks
    i0 = (b-1)*Nwin + 1;
    i1 = min(b*Nwin, N);

    Ubar(b) = mean(tt.U(i0:i1), 'omitnan');
    Vbar(b) = mean(tt.V(i0:i1), 'omitnan');
end

% --- Mean wind speed ---
WindSpeed = sqrt(Ubar.^2 + Vbar.^2);

% --- Mean wind direction (meteorological FROM direction) ---
theta_math = atan2(Vbar, Ubar);          % mathematical angle (rad)
WindDir = mod(270 - rad2deg(theta_math), 360);   % convert to met convention

% --- Add to blocks table (optional) ---
blocks.WindSpeed = WindSpeed;
blocks.WindDir   = WindDir;

% --- Plot wind speed ---
figure('Name','Mean Wind Speed','Color','w');
plot(blocks.time, WindSpeed, '-o', 'LineWidth', 1.5);
ylabel('Wind speed (m/s)');
xlabel('Time');
grid on;

% --- Plot wind direction ---
figure('Name','Mean Wind Direction','Color','w');
plot(blocks.time, WindDir, '-o', 'LineWidth', 1.5);
ylabel('Wind direction (° from)');
xlabel('Time');
ylim([0 360]);
yticks(0:45:360);
grid on;


%% For the whole timeseries now
clc

tt = out.tt;   % from function

% --- Mean horizontal wind components ---
Ubar = mean(tt.U, 'omitnan');
Vbar = mean(tt.V, 'omitnan');

% --- Mean horizontal wind speed ---
WindSpeed_mean = sqrt(Ubar^2 + Vbar^2);

% --- Mean wind direction (meteorological FROM direction) ---
theta_math = atan2(Vbar, Ubar);                 % radians
WindDir_mean = mod(270 - rad2deg(theta_math), 360);

% --- Display results ---
fprintf('Mean wind speed     = %.3f m/s\n', WindSpeed_mean);
fprintf('Mean wind direction = %.1f° (from)\n', WindDir_mean);

% --- Optional: plot instantaneous wind speed and direction ---
WindSpeed_inst = sqrt(tt.U.^2 + tt.V.^2);
theta_math_inst = atan2(tt.V, tt.U);
WindDir_inst = mod(270 - rad2deg(theta_math_inst), 360);

figure('Name','Instantaneous wind speed','Color','w');
plot(tt.t, WindSpeed_inst, 'LineWidth', 1.2);
ylabel('Wind speed (m/s)');
xlabel('Time'); grid on;

figure('Name','Instantaneous wind direction','Color','w');
plot(tt.t, WindDir_inst, 'LineWidth', 1.2);
ylabel('Wind direction (° from)');
xlabel('Time'); ylim([0 360]); grid on;



%% Filtering meaningless data (wind speed is dominated by boat motion)
clc
U  = tt.U;
V  = tt.V;
W = tt.W;

WindSpeed = sqrt(U.^2 + V.^2);
theta_math = atan2(V, U);
WD_inst = mod(270 - rad2deg(theta_math), 360);

% Remove low-speed garbage
WD_inst(WindSpeed < 1.0) = NaN;    % threshold 0.5–1 m/s works well


% Moving mean U and V over 2 minutes
window = round(2*60*out.Fs);  % 2 minutes
Umean = movmean(U, window, 'omitnan');
Vmean = movmean(V, window, 'omitnan');

% Compute direction of the smoothed mean flow
theta_mean = atan2(Vmean, Umean);
WD_mean = mod(270 - rad2deg(theta_mean), 360);


WD_smooth = movmedian(WD_mean, round(10*out.Fs), 'omitnan');

figure;
plot(tt.t, WD_smooth, 'LineWidth', 1.2);
ylabel('Wind direction (° from)');
xlabel('Time');
ylim([0 360]);
grid on;
title('Stabilized Wind Direction (no IMU)');


%% Unwrapping Direction
clc
U  = tt.U;
V  = tt.V;

WindSpeed = sqrt(U.^2 + V.^2);
theta_math = atan2(V, U);
WD_inst = mod(270 - rad2deg(theta_math), 360);

% Remove low-speed garbage
WD_inst(WindSpeed < 1.0) = NaN;

% Moving mean U and V over 2 minutes
window = round(2*60*out.Fs);  % 2 minutes
Umean = movmean(U, window, 'omitnan');
Vmean = movmean(V, window, 'omitnan');

% Direction of the smoothed mean flow (0..360)
theta_mean = atan2(Vmean, Umean);
WD_mean = mod(270 - rad2deg(theta_mean), 360);

% --- UNWRAP (handle NaNs) ---
wd_rad = deg2rad(WD_mean);
wd_unwrapped_rad = unwrap_nan(wd_rad);          % continuous (rad)
WD_unwrapped = rad2deg(wd_unwrapped_rad);       % continuous (deg), can exceed 360

% Smooth in the unwrapped domain (recommended)
WD_unwrapped_smooth = movmedian(WD_unwrapped, round(10*out.Fs), 'omitnan');

% Optional: wrap back to 0..360 for display
WD_wrapped_smooth = mod(WD_unwrapped_smooth, 360);

figure;
plot(tt.t, WD_unwrapped_smooth, 'LineWidth', 1.2);
ylabel('Wind direction (deg, unwrapped)');
xlabel('Time');
grid on;
title('Stabilized Wind Direction (unwrapped, no IMU)');

%% Direction vizualisation improvement
clc
th = deg2rad(WD_mean);
th(WindSpeed < 1.0) = NaN;   % si tu veux conserver le même masque vitesse

c = mean(cos(th), 'omitnan');
s = mean(sin(th), 'omitnan');
WD_ref = mod(270 - rad2deg(atan2(s, c)), 360);   % ref en "deg from" (0..360)

% Déviation signée autour de la ref: [-180, +180]
dWD = wrapTo180(WD_mean - WD_ref);

% (Optionnel) lissage sur la déviation
dWD_smooth = movmedian(dWD, round(10*out.Fs), 'omitnan');

% Reconstruire une "direction" continue autour de la ref (reste dans [ref-180, ref+180])
WD_rel = WD_ref + dWD_smooth;

% Plots
figure;
plot(tt.t, WD_rel, 'LineWidth', 1.2);
ylabel('Wind direction (deg, relative to ref)');
xlabel('Time');
grid on;
title('Wind Direction (relative / no wrap jumps)');

% Si tu préfères afficher seulement l'anomalie:
% figure; plot(tt.t, dWD_smooth); ylabel('\DeltaWD (deg)'); ylim([-180 180]); grid on;



%% Stabilized wind speed from the whole time series (no IMU)

tt = out.tt;   % full high-frequency timetable

U = tt.U;
V = tt.V;

%% Raw instantaneous wind speed
WindSpeed_inst = sqrt(U.^2 + V.^2);   % m/s

% 2. Remove unreliable values (optional)
% If speed is extremely low, direction is meaningless, but speed is still fine.
% Still, remove values that are physically absurd (optional).
WindSpeed_filt = WindSpeed_inst;
WindSpeed_filt(WindSpeed_inst < 0.1) = NaN;   % remove < 0.1 m/s if desired

% 3. Smooth the wind speed to make it stable
% Choose smoothing window depending on sampling rate (out.Fs)
% Example: 10-second moving average
smoothWindow = round(10 * out.Fs);

WindSpeed_smooth = movmean(WindSpeed_filt, smoothWindow, 'omitnan');

% Plot all versions
figure('Name','Wind Speed (raw and smoothed)','Color','w');
plot(tt.t, WindSpeed_inst, 'Color',[0.7 0.7 1]); hold on;
plot(tt.t, WindSpeed_smooth, 'b','LineWidth',1.5);
ylabel('Wind speed (m/s)');
xlabel('Time');
legend('Instantaneous','Smoothed','Location','best');
grid on;

%% Load Data From IMU
clc
%  READ IMU LOG + PLOT POSITION, SPEED, DIRECTION
clear D
D = read_gq7_data_flight("SensorConnectData_21_septembre_long.csv")
D.time = D.time - hours(4);


figure('Color','w','Name','GPS position vs time');

subplot(2,1,1)
plot(D.time, D.lat,'k');
grid on
ylabel('Latitude (deg)')
title('Latitude vs time')

subplot(2,1,2)
plot(D.time, D.lon,'k');
grid on
ylabel('Longitude (deg)')
xlabel('Time')
title('Longitude vs time')

R = 6371000;  % Earth radius (m)

lat0 = deg2rad(D.lat(1));
lon0 = deg2rad(D.lon(1));

xE = R * (deg2rad(D.lon) - lon0) .* cos(lat0);
yN = R * (deg2rad(D.lat) - lat0);

figure('Color','w','Name','Position vs time');

subplot(2,1,1)
plot(D.time, D.lat,'k');
grid on; ylabel('Latitude (deg)');

subplot(2,1,2)
plot(D.time, D.lon,'k');
grid on; ylabel('Longitude (deg)');
xlabel('Time');


figure('Color','w','Name','Platform speed');
plot(D.time, D.speed_imu,'k');
grid on
ylabel('Speed (m/s)')
xlabel('Time')
title('Platform horizontal speed (IMU)')

figure('Color','w','Name','Platform velocity components');

subplot(3,1,1)
plot(D.time, D.ve,'r');
grid on
ylabel('v_E (m/s)')

subplot(3,1,2)
plot(D.time, D.vn,'g');
grid on
ylabel('v_N (m/s)')

subplot(3,1,3)
plot(D.time, -D.vd,'b');  % Up = -Down
grid on
ylabel('v_U (m/s)')
xlabel('Time')


COG = mod(atan2d(D.ve, D.vn), 360);  % deg from North, clockwise

figure('Color','w','Name','Platform direction');
plot(D.time, COG);
ylim([0 360]); yticks(0:45:360)
grid on
ylabel('Direction (deg)')
xlabel('Time')
title('Course over ground (IMU)')


figure('Color','w','Name','IMU attitude');

subplot(3,1,1)
plot(D.time, rad2deg(D.roll),'k');
grid on
ylabel('Roll (deg)')

subplot(3,1,2)
plot(D.time, rad2deg(D.pitch),'k');
grid on
ylabel('Pitch (deg)')

subplot(3,1,3)
plot(D.time, rad2deg(D.yaw),'k');
grid on
ylabel('Yaw (deg)')
xlabel('Time')


%% Clean data from IMU
clc

S = hypot(D.vn, D.ve);
idx = find(S > 10);  % threshold: 10 m/s is already huge for a boat

fprintf('Number of speed spikes >10 m/s: %d\n', numel(idx));
if ~isempty(idx)
    k = idx(1:min(10,end));
    disp(table(D.time(k), D.vn(k), D.ve(k), S(k), 'VariableNames', {'time','vn','ve','speed'}));
end

fprintf('max speed_imu = %.2f\n', max(D.speed_imu,[],'omitnan'));
fprintf('max hypot(vn,ve) = %.2f\n', max(hypot(D.vn, D.ve),[],'omitnan'));



%% Interpolate IMU time series onto sonic time base  (MUST be before offset removal)
% This creates roll_i/pitch_i/yaw_i and VnVeVu aligned to sonic samples.
% It also defines:
%   t_sonic_dt : sonic time as datetime
%   tS, tI     : numeric seconds since a common t0 (for interp1)
%   dt         : sonic sample interval (s)
clc

% Sonic and IMU times
t_sonic = tt.t;
t_imu   = D.time;

% --- 0) Harmonize sonic time as datetime ---
if isdatetime(t_sonic)
    t_sonic_dt = t_sonic(:);
else
    % If sonic time is numeric hours-of-day
    t_sonic = t_sonic(:);
    t_sonic = fillmissing(t_sonic, 'linear', 'EndValues','nearest');
    day0 = dateshift(t_imu(1), 'start', 'day');
    t_sonic_dt = day0 + hours(t_sonic);

    % handle wrap across midnight
    wrapIdx = [false; diff(t_sonic) < 0];
    if any(wrapIdx)
        t_sonic_dt = t_sonic_dt + days(cumsum(wrapIdx));
    end
end

% --- 1) TimeZone mismatch fix ---
if isdatetime(t_sonic_dt) && isdatetime(t_imu)
    tzS = t_sonic_dt.TimeZone;
    tzI = t_imu.TimeZone;

    if isempty(tzS) && ~isempty(tzI)
        t_sonic_dt.TimeZone = tzI;
    elseif ~isempty(tzS) && isempty(tzI)
        t_imu.TimeZone = tzS;
    elseif ~strcmp(tzS, tzI)
        % safest: convert both to UTC
        t_sonic_dt.TimeZone = 'UTC';
        t_imu.TimeZone      = 'UTC';
    end
end

% --- 2) Build numeric time axes (seconds) for interp1 ---
t0 = min(t_sonic_dt(1), t_imu(1));        % datetime
tS = seconds(t_sonic_dt - t0);            % Nx1 double
tI = seconds(t_imu(:)      - t0);         % Mx1 double

% Ensure strictly increasing IMU time for interp1
[ tI, uniqIdx ] = unique(tI, 'stable');

% --- 3) IMU angles (radians) on IMU time base ---
roll_raw  = D.roll(:);
pitch_raw = D.pitch(:);
yaw_raw   = D.yaw(:);

roll_raw  = fillmissing(roll_raw,  'linear', 'EndValues','nearest');
pitch_raw = fillmissing(pitch_raw, 'linear', 'EndValues','nearest');
yaw_raw   = fillmissing(yaw_raw,   'linear', 'EndValues','nearest');

roll_raw  = roll_raw(uniqIdx);
pitch_raw = pitch_raw(uniqIdx);
yaw_raw   = yaw_raw(uniqIdx);

% --- 4) Interpolate IMU angles onto sonic time base (radians) ---
roll_i  = interp1(tI, roll_raw,  tS, 'pchip', 'extrap');
pitch_i = interp1(tI, pitch_raw, tS, 'pchip', 'extrap');
yaw_i   = interp1(tI, yaw_raw,   tS, 'pchip', 'extrap');

% --- 5) Interpolate IMU velocities (North/East/Up) onto sonic time base ---
Vn_imu = D.vn(:);
Ve_imu = D.ve(:);
Vu_imu = -D.vd(:);   % Up = -Down

Vn_imu = fillmissing(Vn_imu, 'linear', 'EndValues','nearest');
Ve_imu = fillmissing(Ve_imu, 'linear', 'EndValues','nearest');
Vu_imu = fillmissing(Vu_imu, 'linear', 'EndValues','nearest');

Vn_imu = Vn_imu(uniqIdx);
Ve_imu = Ve_imu(uniqIdx);
Vu_imu = Vu_imu(uniqIdx);

Vn_i = interp1(tI, Vn_imu, tS, 'pchip', 'extrap');
Ve_i = interp1(tI, Ve_imu, tS, 'pchip', 'extrap');
Vu_i = interp1(tI, Vu_imu, tS, 'pchip', 'extrap');

VnVeVu = [Vn_i, Ve_i, Vu_i];

% --- 6) Sonic dt (seconds) ---
dt = median(diff(tS), 'omitnan');

% --- Sanity prints ---
fprintf('Interpolation done: N_sonic=%d samples | N_imu=%d samples\n', numel(tS), numel(tI));
fprintf('roll_i  : finite=%.1f%%, min=%.3g, max=%.3g\n', 100*mean(isfinite(roll_i)),  min(roll_i,[],'omitnan'),  max(roll_i,[],'omitnan'));
fprintf('pitch_i : finite=%.1f%%, min=%.3g, max=%.3g\n', 100*mean(isfinite(pitch_i)), min(pitch_i,[],'omitnan'), max(pitch_i,[],'omitnan'));
fprintf('yaw_i   : finite=%.1f%%, min=%.3g, max=%.3g\n', 100*mean(isfinite(yaw_i)),   min(yaw_i,[],'omitnan'),   max(yaw_i,[],'omitnan'));
fprintf('dt (sonic) = %.4f s (Fs ~ %.2f Hz)\n', dt, 1/dt);

%% Lets remove the mean from each angle of IMU before motion correct
clc
% Goal: remove the fixed camera mounting tilt (e.g., ~55 deg) so that the
% lever-arm correction uses anemometer attitude variations rather than camera bias.



% --- robust "mean" = median (less sensitive to spikes)
roll0  = median(roll_i,  'omitnan');
pitch0 = median(pitch_i, 'omitnan');

% --- anemometer attitude variations (radians)
roll_anemometer  = roll_i  - roll0;
pitch_anemometer = pitch_i - pitch0;

% --- yaw: usually keep as-is (heading). Only remove a known mounting yaw offset.
yaw_anemometer = yaw_i;  % start with this

% OPTIONAL: if a fixed yaw mounting offset (deg), apply it:
% yaw_mount_deg = 0; % <-- set if known
% yaw_anemometer = yaw_anemometer - deg2rad(yaw_mount_deg);

% Build rotation matrix input for motion correction (radians)
rot_mat = [roll_anemometer, pitch_anemometer, yaw_anemometer];

% rot_mat(:,3) = -rot_mat(:,3);

% --- Then apply  time-lag shift AFTER this (so lag acts on the final angles used):
% rot_mat = shift_by_samples(rot_mat, lag);
% VnVeVu  = shift_by_samples(VnVeVu,  lag);

% --- sanity prints (optional)

fprintf('Removed fixed offsets: roll0=%.2f deg, pitch0=%.2f deg\n', rad2deg(roll0), rad2deg(pitch0));

fprintf('roll_anemometer : mean=%.3f deg | std=%.3f deg | min=%.3f deg | max=%.3f deg\n', ...
    mean(rad2deg(roll_anemometer),'omitnan'), std(rad2deg(roll_anemometer),'omitnan'), ...
    min(rad2deg(roll_anemometer),[],'omitnan'), max(rad2deg(roll_anemometer),[],'omitnan'));

fprintf('pitch_anemometer: mean=%.3f deg | std=%.3f deg | min=%.3f deg | max=%.3f deg\n', ...
    mean(rad2deg(pitch_anemometer),'omitnan'), std(rad2deg(pitch_anemometer),'omitnan'), ...
    min(rad2deg(pitch_anemometer),[],'omitnan'), max(rad2deg(pitch_anemometer),[],'omitnan'));

fprintf('yaw_anemometer  : mean=%.3f deg | std=%.3f deg | min=%.3f deg | max=%.3f deg\n', ...
    mean(rad2deg(yaw_anemometer),'omitnan'), std(rad2deg(yaw_anemometer),'omitnan'), ...
    min(rad2deg(yaw_anemometer),[],'omitnan'), max(rad2deg(yaw_anemometer),[],'omitnan'));

figure('Color','w','Name','Anemometer attitude (offset removed)');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile; plot(t_sonic_dt, rad2deg(roll_anemometer));  grid on; ylabel('roll\_an (deg)');
nexttile; plot(t_sonic_dt, rad2deg(pitch_anemometer)); grid on; ylabel('pitch\_an (deg)');
nexttile; plot(t_sonic_dt, rad2deg(yaw_anemometer));   grid on; ylabel('yaw\_an (deg)'); xlabel('Time');


%% ===================== MOTION CORRECTION + RESULTS (corr3b) =====================
%   tt.t   : Nx1 datetime (sonic)
%
%   D.time : Mx1 datetime (IMU)
%   D.roll, D.pitch, D.yaw : Mx1 angles in radians (IMU)
%   D.vn, D.ve, D.vd : Mx1 velocities (North, East, Down) in m/s at IMU
%
% And have these functions on  path:
%   - motcorr_deluxe_insvel.m
%   - rot3d.m
%   - diff1byf.m
%   - shift_by_samples.m

clc;

fs  = 20;      % sonic sampling rate
hf  = 2;       % angle low-pass / rate stability (Hz)
hp  = 0.02;    % high-pass for diagnostic correlations (Hz)
lag = 3;       % samples (20 Hz -> 0.15 s).  can change after diagnostics.

% -------------------- Inputs --------------------
% Sonic time
t_sonic = tt.t;

U = tt.U;
V = tt.V;
W = tt.W;


% IMU angles (radians) + velocities (N/E/Up)
roll  = roll_anemometer;    % rad
pitch = pitch_anemometer;   % rad
yaw   = yaw_anemometer;     % rad

Vn_imu = D.vn;     % m/s
Ve_imu = D.ve;     % m/s
Vu_imu = -D.vd;    % m/s (Up = -Down)

t_imu = D.time;
% -------------------- 0-2) Time harmonization + interpolation --------------------
% Already done earlier (before offset removal). Using existing variables:
%   t_sonic_dt, tS, tI, roll_i/pitch_i/yaw_i, VnVeVu, dt, rot_mat

% -------------------- 3)  yaw sign fix --------------------
rot_mat(:,3) = -rot_mat(:,3);

% -------------------- 4) Apply time-lag shift (IMU series -> aligned to sonic) --------------------
rot_mat = shift_by_samples(rot_mat, lag);
VnVeVu  = shift_by_samples(VnVeVu,  lag);

% -------------------- 5) Sonic wind matrix (BODY frame) --------------------
vel_mat = [double(U(:)), double(V(:)), double(W(:))];

% Force sonic W = Up (z-up convention used here)
vel_mat(:,3) = -vel_mat(:,3);

% -------------------- 6) Lever arm IMU -> sonic in BODY frame (meters) --------------------
% Anemometer is starboard (Ly>0) and ABOVE IMU (Lz>0) in this z-up convention
Lx = 0.0;
Ly = 0.5;
Lz = 2.5;
L  = [Lx, Ly, Lz];

% -------------------- 7) Motion correction --------------------
corr3b = motcorr_deluxe_insvel(vel_mat, rot_mat, VnVeVu, L, dt, hf);

% ===================== RESULTS / DIAGNOSTICS =====================
t = t_sonic_dt;

% Raw rotated (Earth) vs corrected (Earth, N/W/Up)
UrawN = corr3b.ue;     UrawW = corr3b.ve;     UrawU = corr3b.we;
UcorN = corr3b.U_N;    UcorW = corr3b.U_W;    UcorU = corr3b.U_up;

% Platform velocities in N/W/Up (W = -E)
Vn = VnVeVu(:,1);
Ve = VnVeVu(:,2);
Vu = VnVeVu(:,3);
Vw = -Ve;

% 1) Overlays
figure; plot(t,UrawN); hold on; plot(t,UcorN); grid on;
xlabel('Time'); ylabel('N (m/s)'); legend('Raw rotated','Corrected','Location','best');
title('Wind North: raw rotated vs motion-corrected');

figure; plot(t,UrawW); hold on; plot(t,UcorW); grid on;
xlabel('Time'); ylabel('W (m/s)'); legend('Raw rotated','Corrected','Location','best');
title('Wind West: raw rotated vs motion-corrected');

figure; plot(t,UrawU); hold on; plot(t,UcorU); grid on;
xlabel('Time'); ylabel('Up (m/s)'); legend('Raw rotated','Corrected','Location','best');
title('Wind Up: raw rotated vs motion-corrected');

% 2) Wind speed
Spd_raw = sqrt(UrawN.^2 + UrawW.^2 + UrawU.^2);
Spd_cor = sqrt(UcorN.^2 + UcorW.^2 + UcorU.^2);

figure; plot(t,Spd_raw); hold on; plot(t,Spd_cor); grid on;
xlabel('Time'); ylabel('|U| (m/s)');
legend('|U| raw rotated','|U| corrected','Location','best');
title('Wind speed');

% 3) Raw->Corrected correlation with platform velocity (finite-only) USING corrcoef (no corr() issues)
maskN = isfinite(UrawN) & isfinite(UcorN) & isfinite(Vn);
maskW = isfinite(UrawW) & isfinite(UcorW) & isfinite(Vw);
maskU = isfinite(UrawU) & isfinite(UcorU) & isfinite(Vu);

C = corrcoef(UrawN(maskN), Vn(maskN)); r_raw_N = C(1,2);
C = corrcoef(UcorN(maskN), Vn(maskN)); r_cor_N = C(1,2);

C = corrcoef(UrawW(maskW), Vw(maskW)); r_raw_W = C(1,2);
C = corrcoef(UcorW(maskW), Vw(maskW)); r_cor_W = C(1,2);

C = corrcoef(UrawU(maskU), Vu(maskU)); r_raw_U = C(1,2);
C = corrcoef(UcorU(maskU), Vu(maskU)); r_cor_U = C(1,2);

fprintf('Corr with platform vel (raw->corrected): N %.3f -> %.3f | W %.3f -> %.3f | Up %.3f -> %.3f\n', ...
    r_raw_N, r_cor_N, r_raw_W, r_cor_W, r_raw_U, r_cor_U);

% 4) High-pass correlation (diagnostic of motion imprint)
Vn2 = double(Vn(:)); Vw2 = double(Vw(:)); Vu2 = double(Vu(:));
UN2 = double(UcorN(:)); UW2 = double(UcorW(:)); UU2 = double(UcorU(:));

Vn2(~isfinite(Vn2))=NaN; Vn2=fillmissing(Vn2,'linear','EndValues','nearest'); Vn2(~isfinite(Vn2))=0;
Vw2(~isfinite(Vw2))=NaN; Vw2=fillmissing(Vw2,'linear','EndValues','nearest'); Vw2(~isfinite(Vw2))=0;
Vu2(~isfinite(Vu2))=NaN; Vu2=fillmissing(Vu2,'linear','EndValues','nearest'); Vu2(~isfinite(Vu2))=0;

UN2(~isfinite(UN2))=NaN; UN2=fillmissing(UN2,'linear','EndValues','nearest'); UN2(~isfinite(UN2))=0;
UW2(~isfinite(UW2))=NaN; UW2=fillmissing(UW2,'linear','EndValues','nearest'); UW2(~isfinite(UW2))=0;
UU2(~isfinite(UU2))=NaN; UU2=fillmissing(UU2,'linear','EndValues','nearest'); UU2(~isfinite(UU2))=0;

rN = corrcoef(highpass(Vn2,hp,fs), highpass(UN2,hp,fs)); rN = rN(1,2);
rW = corrcoef(highpass(Vw2,hp,fs), highpass(UW2,hp,fs)); rW = rW(1,2);
rU = corrcoef(highpass(Vu2,hp,fs), highpass(UU2,hp,fs)); rU = rU(1,2);

fprintf('High-pass corr(corrected, platform): N %.3f | W %.3f | Up %.3f\n', rN, rW, rU);

% 5) Removed terms magnitude (translation vs rotation)
Utrans = sqrt(corr3b.ube.^2    + corr3b.vbe.^2    + corr3b.wbe.^2);
Urot   = sqrt(corr3b.uberot.^2 + corr3b.vberot.^2 + corr3b.wberot.^2);

figure; plot(t,Utrans); hold on; plot(t,Urot); grid on;
xlabel('Time'); ylabel('Speed (m/s)');
legend('|V_{trans}| (INS)','|Ω×r|','Location','best');
title('Motion terms removed');


fprintf('D.roll: min=%.3f, max=%.3f\n', min(D.roll,[],'omitnan'), max(D.roll,[],'omitnan'));
fprintf('D.pitch: min=%.3f, max=%.3f\n', min(D.pitch,[],'omitnan'), max(D.pitch,[],'omitnan'));
fprintf('D.yaw: min=%.3f, max=%.3f\n', min(D.yaw,[],'omitnan'), max(D.yaw,[],'omitnan'));

%% Now compare results
clc
% ===================== COMPARE RAW vs MOTION-CORRECTED WIND =====================
% Requires: tt (with U,V,W), corr3b (from motcorr_deluxe_insvel), VnVeVu, t_sonic_dt (or tt.t)

% --- time vector ---
if exist('t_sonic_dt','var')
    t = t_sonic_dt;
else
    t = tt.t;
end

% --- RAW wind from sonic (choose a consistent convention) ---
Uraw_E = double(tt.U(:));          % check:  U,V are horizontal components
Uraw_N = double(tt.V(:));          % if  sonic is already N/E, swap accordingly!
Uraw_Up = -double(tt.W(:));        % Up = -Down (same convention used)

% --- MOTION-CORRECTED wind from corr3b ---
% corr3b fields (as used earlier):
Ucor_N  = double(corr3b.U_N(:));
Ucor_W  = double(corr3b.U_W(:));
Ucor_Up = double(corr3b.U_up(:));

% NOTE: raw uses (E,N,Up) above, corrected uses (N,W,Up) here.
% To compare apples-to-apples, we need matching axes. If raw U,V are actually (N,W),
% then set Uraw_N=tt.U, Uraw_W=tt.V, etc. Use the next block to pick one.

% --------- PICK CONSISTENT AXES (recommended) ----------
% If  sonic U,V are actually "North" and "West" already ( earlier labels suggest that),
% uncomment this and comment the E/N assignment above.

% Uraw_N  = double(tt.U(:));
% Uraw_W  = double(tt.V(:));
% Uraw_Up = -double(tt.W(:));

% If  keep the first assignment (E/N),  MUST map corrected W->E etc.
% For now, we assume raw U,V are (N,W) because  plots were "Wind North/West".

Uraw_N  = double(tt.U(:));      % <-- assumes tt.U corresponds to North (per  earlier labeling)
Uraw_W  = double(tt.V(:));      % <-- assumes tt.V corresponds to West
Uraw_Up = -double(tt.W(:));  

% --- speeds ---
Spd_raw = sqrt(Uraw_N.^2 + Uraw_W.^2 + Uraw_Up.^2);
Spd_cor = sqrt(Ucor_N.^2 + Ucor_W.^2 + Ucor_Up.^2);

% --------- PLOTS: components ----------
figure('Color','w','Name','Raw vs Motion-corrected wind (components)');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile
plot(t, Uraw_N, 'k'); hold on;
plot(t, Ucor_N, 'r');
grid on; ylabel('N (m/s)');
legend('Raw','Corrected','Location','best');

nexttile
plot(t, Uraw_W, 'k'); hold on;
plot(t, Ucor_W, 'r');
grid on; ylabel('W (m/s)');
legend('Raw','Corrected','Location','best');

nexttile
plot(t, Uraw_Up, 'k'); hold on;
plot(t, Ucor_Up, 'r');
grid on; ylabel('Up (m/s)'); xlabel('Time');
legend('Raw','Corrected','Location','best');

% --------- PLOT: speed ----------
figure('Color','w','Name','Raw vs Motion-corrected wind (speed)');
plot(t, Spd_raw, 'k'); hold on;
plot(t, Spd_cor, 'r');
grid on; ylabel('|U| (m/s)'); xlabel('Time');
legend('Raw','Corrected','Location','best');
title('Wind speed magnitude');

% --------- SUMMARY STATS ----------
fprintf('\n=== Summary stats ===\n');
fprintf('Mean N raw=%.3f | cor=%.3f   (m/s)\n', mean(Uraw_N,'omitnan'), mean(Ucor_N,'omitnan'));
fprintf('Mean W raw=%.3f | cor=%.3f   (m/s)\n', mean(Uraw_W,'omitnan'), mean(Ucor_W,'omitnan'));
fprintf('Mean Up raw=%.3f | cor=%.3f  (m/s)\n', mean(Uraw_Up,'omitnan'), mean(Ucor_Up,'omitnan'));
fprintf('Mean |U| raw=%.3f | cor=%.3f (m/s)\n', mean(Spd_raw,'omitnan'), mean(Spd_cor,'omitnan'));

fprintf('Std  N raw=%.3f | cor=%.3f   (m/s)\n', std(Uraw_N,'omitnan'), std(Ucor_N,'omitnan'));
fprintf('Std  W raw=%.3f | cor=%.3f   (m/s)\n', std(Uraw_W,'omitnan'), std(Ucor_W,'omitnan'));
fprintf('Std  Up raw=%.3f | cor=%.3f  (m/s)\n', std(Uraw_Up,'omitnan'), std(Ucor_Up,'omitnan'));

% --------- HOW MUCH PLATFORM VELOCITY GOT REMOVED? (correlation test) ----------
% IMU velocity interpolated to sonic base:
Vn = double(VnVeVu(:,1));
Ve = double(VnVeVu(:,2));
Vu = double(VnVeVu(:,3));

% If  corrected horizontal is (N,W), compare W with platform WEST velocity:
Vw = -Ve;

maskN = isfinite(Uraw_N) & isfinite(Ucor_N) & isfinite(Vn);
maskW = isfinite(Uraw_W) & isfinite(Ucor_W) & isfinite(Vw);
maskU = isfinite(Uraw_Up) & isfinite(Ucor_Up) & isfinite(Vu);

r_raw_N = corr(Uraw_N(maskN), Vn(maskN), 'rows','complete');
r_cor_N = corr(Ucor_N(maskN), Vn(maskN), 'rows','complete');

r_raw_W = corr(Uraw_W(maskW), Vw(maskW), 'rows','complete');
r_cor_W = corr(Ucor_W(maskW), Vw(maskW), 'rows','complete');

r_raw_U = corr(Uraw_Up(maskU), Vu(maskU), 'rows','complete');
r_cor_U = corr(Ucor_Up(maskU), Vu(maskU), 'rows','complete');

fprintf('\n=== Correlation with platform velocity (want these to drop after correction) ===\n');
fprintf('North: raw r=%.3f  -> corrected r=%.3f\n', r_raw_N, r_cor_N);
fprintf('West : raw r=%.3f  -> corrected r=%.3f\n', r_raw_W, r_cor_W);
fprintf('Up   : raw r=%.3f  -> corrected r=%.3f\n', r_raw_U, r_cor_U);

% --------- OPTIONAL: turbulence-band comparison (high-pass) ----------
% This checks if the correction mostly removed low-f motion without killing turbulence.
fs = 1/median(seconds(diff(t)),'omitnan');   % estimate from time axis if datetime
hp = 0.02;                                   % Hz (adjust)

UrawN_hp = highpass(fillmissing(Uraw_N,'linear','EndValues','nearest'), hp, fs);
UcorN_hp = highpass(fillmissing(Ucor_N,'linear','EndValues','nearest'), hp, fs);

UrawW_hp = highpass(fillmissing(Uraw_W,'linear','EndValues','nearest'), hp, fs);
UcorW_hp = highpass(fillmissing(Ucor_W,'linear','EndValues','nearest'), hp, fs);

fprintf('\n=== High-pass std (%.3f Hz) ===\n', hp);
fprintf('North: std raw=%.3f | cor=%.3f\n', std(UrawN_hp,'omitnan'), std(UcorN_hp,'omitnan'));
fprintf('West : std raw=%.3f | cor=%.3f\n', std(UrawW_hp,'omitnan'), std(UcorW_hp,'omitnan'));


%% ===== Wind SPEED + DIRECTION: RAW vs MOTION-CORRECTED (U=N, V=W) =====
% Requires: tt, corr3b, and either t_sonic_dt or tt.t
% Convention:
%   U = North (m/s)
%   V = West  (m/s)
% Direction computed as meteorological "FROM" (0..360, 0=N, 90=E, 180=S, 270=W)

% 5) Removed terms magnitude (translation vs rotation)
Utrans = sqrt(corr3b.ube.^2    + corr3b.vbe.^2    + corr3b.wbe.^2);
Urot   = sqrt(corr3b.uberot.^2 + corr3b.vberot.^2 + corr3b.wberot.^2);

figure; plot(t,Utrans); hold on; plot(t,Urot); grid on;
xlabel('Time'); ylabel('Speed (m/s)');
legend('|V_{trans}| (INS)','|Ω×r|','Location','best');
title('Motion terms removed');


% Time
if exist('t_sonic_dt','var')
    t = t_sonic_dt;
else
    t = tt.t;
end

% --- RAW horizontal components ---
Uraw = double(tt.U(:));       % North
Vraw = double(tt.V(:));       % West

% --- CORRECTED horizontal components ---
Ucor = double(corr3b.U_N(:)); % North

if isfield(corr3b,'U_W')
    Vcor = double(corr3b.U_W(:));       % West
elseif isfield(corr3b,'U_E')
    Vcor = -double(corr3b.U_E(:));      % East->West
else
    error('corr3b missing field U_W or U_E');
end

% --- Speeds (horizontal) ---
Spd_raw = hypot(Uraw, Vraw);
Spd_cor = hypot(Ucor, Vcor);

% --- Directions: meteorological FROM ---
% For NE components (N,E): WD_from = mod(270 - atan2d(E,N), 360)
% We have (N,W) so E = -W = -V
Eraw = -Vraw;
Ecor = -Vcor;

WD_raw = mod(270 - atan2d(Eraw, Uraw), 360);

WD_cor = mod(270 - atan2d(Ecor, Ucor), 360);

% --- Mask direction when speed is too small ---
minSpd = 1.0;  % m/s (change if  want)
WD_raw(Spd_raw < minSpd) = NaN;
WD_cor(Spd_cor < minSpd) = NaN;

% --- Optional: unwrap + smooth direction to avoid wrap jumps ---
% (smooth in unwrapped space, then wrap back to 0..360)
win = round(10 * out.Fs);  % if  have out.Fs; otherwise set e.g. win=200 for 20 Hz
if ~exist('out','var') || ~isfield(out,'Fs')
    win = 200; % fallback for ~20 Hz data
end

wd_raw_u = rad2deg(unwrap_nan(deg2rad(WD_raw)));
wd_cor_u = rad2deg(unwrap_nan(deg2rad(WD_cor)));

WD_raw_s = mod(movmedian(wd_raw_u, win, 'omitnan'), 360);
WD_cor_s = mod(movmedian(wd_cor_u, win, 'omitnan'), 360);

Spd_raw_s = movmean(Spd_raw, win, 'omitnan');
Spd_cor_s = movmean(Spd_cor, win, 'omitnan');

% ---- PLOTS ----
figure('Color','w','Name','Wind speed (horizontal): Raw vs Corrected');
plot(t, Spd_raw, 'g'); hold on;
plot(t, Spd_cor, 'r');
plot(t, Spd_raw_s, 'k', 'LineWidth', 1.5);
plot(t, Spd_cor_s, 'b', 'LineWidth', 1.5);
grid on; ylabel('Speed (m/s)'); xlabel('Time');
legend('Raw','Corrected','Raw smoothed','Corrected smoothed','Location','best');
title('Horizontal wind speed');

figure('Color','w','Name','Wind direction (FROM): Raw vs Corrected');
plot(t, WD_raw, 'g'); hold on;
plot(t, WD_cor, 'r');
plot(t, WD_raw_s, 'k', 'LineWidth', 1.5);
plot(t, WD_cor_s, 'b', 'LineWidth', 1.5);
grid on; ylabel('Direction (deg FROM)'); xlabel('Time');
ylim([0 360]); yticks(0:45:360);
legend('Raw','Corrected','Raw smoothed','Corrected smoothed','Location','best');
title('Wind direction (meteorological FROM)');

% ---- NUMBERS ----
fprintf('\n=== Horizontal wind speed ===\n');
fprintf('Mean speed: raw %.3f | cor %.3f (m/s)\n', mean(Spd_raw,'omitnan'), mean(Spd_cor,'omitnan'));
fprintf('Std  speed: raw %.3f | cor %.3f (m/s)\n', std(Spd_raw,'omitnan'),  std(Spd_cor,'omitnan'));

% Circular mean direction (robust)
th_raw = deg2rad(WD_raw); th_cor = deg2rad(WD_cor);
c_raw = mean(cos(th_raw),'omitnan'); s_raw = mean(sin(th_raw),'omitnan');
c_cor = mean(cos(th_cor),'omitnan'); s_cor = mean(sin(th_cor),'omitnan');
WD_raw_mean = mod(rad2deg(atan2(s_raw,c_raw)), 360);
WD_cor_mean = mod(rad2deg(atan2(s_cor,c_cor)), 360);

fprintf('\n=== Wind direction (FROM) ===\n');
fprintf('Circular mean dir: raw %.1f | cor %.1f (deg)\n', WD_raw_mean, WD_cor_mean);

%% Moving mean to filter wobbling
clc

% ===== define dt =====
if exist('t_sonic_dt','var')
    t = t_sonic_dt;
elseif ismember('t', tt.Properties.VariableNames)
    t = tt.t;
else
    error('No time vector found (expected t_sonic_dt or tt.t)');
end
dt = median(seconds(diff(t)), 'omitnan');

% ===== moving mean smoothing =====
Tw = 90;                    % seconds
Nw = max(5, round(Tw/dt));  % samples

UcN  = double(corr3b.U_N(:));
UcW  = double(corr3b.U_W(:));
UcUp = double(corr3b.U_up(:));

UcN_s  = movmean(UcN,  Nw, 'omitnan');
UcW_s  = movmean(UcW,  Nw, 'omitnan');
UcUp_s = movmean(UcUp, Nw, 'omitnan');

Spd_raw = sqrt(UcN.^2 + UcW.^2 + UcUp.^2);
Spd_s   = sqrt(UcN_s.^2 + UcW_s.^2 + UcUp_s.^2);

% Wind direction (meteorological FROM), using East = -West
E_s  = -UcW_s;
WD_s = mod(270 - atan2d(E_s, UcN_s), 360);

figure('Color','w');
subplot(2,1,1);
plot(t, Spd_raw, 'k'); hold on; plot(t, Spd_s, 'r', 'LineWidth', 1.2);
grid on; ylabel('|U| (m/s)'); legend('Corrected','Movmean','Location','best');

subplot(2,1,2);
plot(t, WD_s, 'r'); grid on; ylabel('Dir (deg FROM)'); xlabel('Time'); ylim([0 360]);

%% ===== Run process_sonic_log_blocks on MOVMEAN-filtered data (one block) =====

% Requires in workspace:
%   tt : timetable with tt.t (datetime) and variables U,V,W,T (already motion-corrected)

% ---- settings ----
Tw      = 60;      % moving-mean window [s]
blockMin= 10;      % block length [min] (Inf = whole series)
z_m     = 3;     % sonic height above surface [m]  <-- EDIT

% ---- dt and Fs from tt.t ----
dt = seconds(diff(tt.t));
dt = dt(~isnan(dt) & dt > 0);
dt = median(dt);
Fs = 1/dt;

% ---- moving mean window (samples) ----
Nw = max(5, round(Tw/dt));

% ---- build movmean-filtered timetable (keep same timestamps) ----
tt_mm = tt;

tt_mm.U = movmean(double(tt.U), Nw, 'omitnan');
tt_mm.V = movmean(double(tt.V), Nw, 'omitnan');
tt_mm.W = movmean(double(tt.W), Nw, 'omitnan');
tt_mm.T = movmean(double(tt.T), Nw, 'omitnan');

% ---- run MOST/u* blocks on smoothed data ----
out_mm = process_sonic_blocks_from_tt(tt_mm, blockMin, z_m);

% ---- quick sanity plots ----
Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

figure('Color','w');
subplot(2,1,1);
plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);
grid on; ylabel('|U| (m/s)');
legend('Raw (Motion Corrected)','Moving Mean Filter','Location','best');
title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));

subplot(2,1,2);
plot(out_mm.blocks.time, out_mm.blocks.u_star, 'o-');
grid on; ylabel('u* (m/s)'); xlabel('Time');
title(sprintf('Block stats (blockMin=%g min, z=%.2f m)', blockMin, z_m));

% ---- plot u*, L, z0 ----
figure('Color','w','Name','MOST outputs (movmean input)');
tB = out_mm.blocks.time;

subplot(3,1,1);
plot(tB, out_mm.blocks.u_star, 'o-'); grid on;
ylabel('u_* (m/s)'); title('Block MOST outputs');

subplot(3,1,2);
plot(tB, out_mm.blocks.L, 'o-'); grid on;
ylabel('L_{MO} (m)');

subplot(3,1,3);
plot(tB, out_mm.blocks.z0, 'o-'); grid on;
ylabel('z_0 (m)'); xlabel('Time');
set(gca,'YScale','log');  % z0 often spans decades

%% ZOOM
% --- Plot raw vs movmean, with a zoomed subplot (11:37 to 11:44) ---

Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

figure('Color','w');

% ===== Top: full series =====
subplot(2,1,1);
plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);
grid on; ylabel('|U| (m/s)');
legend('Raw (Motion Corrected)','Moving Mean Filter','Location','best');
title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));

% ===== Bottom: zoom 11:37–11:44 =====
subplot(2,1,2);

% Build zoom interval using the same date as your data
d0 = dateshift(tt.t(1), 'start', 'day');

t1 = d0 + hours(11) + minutes(37);
t2 = d0 + hours(11) + minutes(44);

% If tt.t has timezone, keep it consistent
if ~isempty(tt.t.TimeZone)
    t1.TimeZone = tt.t.TimeZone;
    t2.TimeZone = tt.t.TimeZone;
end

idx_raw = (tt.t    >= t1) & (tt.t    <= t2);
idx_mm  = (tt_mm.t >= t1) & (tt_mm.t <= t2);

plot(tt.t(idx_raw), Spd_raw(idx_raw), 'k'); hold on;
plot(tt_mm.t(idx_mm), Spd_mm(idx_mm), 'r', 'LineWidth', 1.2);
grid on; ylabel('|U| (m/s)'); xlabel('Time');
title('Zoom: 11:37 to 11:44');
xlim([t1 t2]);

% Optional: match y-limits to data in zoom for nicer scaling
ymin = min([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
ymax = max([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
if isfinite(ymin) && isfinite(ymax) && ymin < ymax
    ylim([ymin ymax] + 0.05*[-1 1]*(ymax - ymin));
end


%% Adding dots to plot to show measurements
clc
% ============================================================
% Raw vs movmean speed plot + zoom panel + green timestamp dots
% Requires in workspace:
%   tt    : timetable with tt.t (datetime) and variables U,V,W (motion-corrected)
%   tt_mm : timetable same timestamps, movmean-filtered U,V,W (as you built)
%   Tw, Fs : your moving-mean window [s] and sampling freq [Hz] (optional for title)
% ============================================================

% ---- speeds ----
Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

% ---- timestamps (from folder names) ----
markStr = [
"2025-09-21 10-55-19.320"
"2025-09-21 11-00-48.776"
"2025-09-21 11-09-12.843"
"2025-09-21 11-14-48.914"
"2025-09-21 11-22-14.107"
"2025-09-21 11-25-25.411"
"2025-09-21 11-31-20.321"
"2025-09-21 11-36-06.521"
"2025-09-21 11-39-14.973"
"2025-09-21 11-42-21.662"
"2025-09-21 11-45-26.115"
"2025-09-21 11-48-41.689"
"2025-09-21 11-51-52.353"
"2025-09-21 11-55-26.302"
"2025-09-21 11-58-37.140"
"2025-09-21 12-01-41.780"
];

tMark = datetime(markStr, "InputFormat","yyyy-MM-dd HH-mm-ss.SSS");
if ~isempty(tt.t.TimeZone)
    tMark.TimeZone = tt.t.TimeZone;
end

% keep only markers inside time span
tMark = tMark(tMark >= tt.t(1) & tMark <= tt.t(end));

% y-value at those times (nearest raw sample)
yMark = interp1(posixtime(tt.t), Spd_raw, posixtime(tMark), "nearest", NaN);

% ---- zoom interval (same day as data) ----
d0 = dateshift(tt.t(1), 'start', 'day');
t1 = d0 + hours(11) + minutes(37);
t2 = d0 + hours(11) + minutes(44);
if ~isempty(tt.t.TimeZone)
    t1.TimeZone = tt.t.TimeZone;
    t2.TimeZone = tt.t.TimeZone;
end

idx_raw = (tt.t    >= t1) & (tt.t    <= t2);
idx_mm  = (tt_mm.t >= t1) & (tt_mm.t <= t2);

% ---- figure ----
figure('Color','w');

% ===== Top: full series + green dots =====
subplot(2,1,1);
plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);
plot(tMark, yMark, '.', 'MarkerSize', 18, 'Color', [0 0.6 0], ...
     'DisplayName','Polarimetric Measurement (45s)');
grid on; ylabel('|U| (m/s)');
% Force x-limits to actual data span (no empty tail)
xlim([tt.t(1) tt.t(end)]);

legend('Raw (Motion Corrected)','Moving Mean Filter','Polarimetric Measurement (45s)','Location','best');
if exist('Tw','var') && exist('Fs','var')
    title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));
else
    title('Raw vs Moving Mean (full series)');
end

% ===== Bottom: zoom 11:37–11:44 =====
subplot(2,1,2);
plot(tt.t(idx_raw), Spd_raw(idx_raw), 'k'); hold on;
plot(tt_mm.t(idx_mm), Spd_mm(idx_mm), 'r', 'LineWidth', 1.2);
grid on; ylabel('|U| (m/s)'); xlabel('Time');
title('Zoom: 11:37 to 11:44');
xlim([t1 t2]);

% Optional: nicer y-lims in zoom
ymin = min([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
ymax = max([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
if isfinite(ymin) && isfinite(ymax) && ymin < ymax
    ylim([ymin ymax] + 0.05*[-1 1]*(ymax - ymin));
end


%% Adding green dots on both + zoom + redline
clc

% ============================================================
% Raw vs movmean speed plot + zoom panel + green dots on RED line
% Requires in workspace:
%   tt    : timetable with tt.t (datetime) and variables U,V,W (motion-corrected)
%   tt_mm : timetable same timestamps, movmean-filtered U,V,W (as you built)
%   Tw, Fs : optional for title
% ============================================================

% ---- speeds ----
Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

% ---- timestamps (from folder names) ----
markStr = [
"2025-09-21 10-55-19.320"
"2025-09-21 11-00-48.776"
"2025-09-21 11-09-12.843"
"2025-09-21 11-14-48.914"
"2025-09-21 11-22-14.107"
"2025-09-21 11-25-25.411"
"2025-09-21 11-31-20.321"
"2025-09-21 11-36-06.521"
"2025-09-21 11-39-14.973"
"2025-09-21 11-42-21.662"
"2025-09-21 11-45-26.115"
"2025-09-21 11-48-41.689"
"2025-09-21 11-51-52.353"
"2025-09-21 11-55-26.302"
"2025-09-21 11-58-37.140"
"2025-09-21 12-01-41.780"
];

tMark = datetime(markStr, "InputFormat","yyyy-MM-dd HH-mm-ss.SSS");
if ~isempty(tt.t.TimeZone)
    tMark.TimeZone = tt.t.TimeZone;
end

% keep only markers inside time span
tMark = tMark(tMark >= tt_mm.t(1) & tMark <= tt_mm.t(end));

% y-value at those times on the RED line (moving-mean series)
yMark_mm = interp1(posixtime(tt_mm.t), Spd_mm, posixtime(tMark), "nearest", NaN);

% ---- zoom interval (same day as data) ----
d0 = dateshift(tt_mm.t(1), 'start', 'day');
t1 = d0 + hours(11) + minutes(37);
t2 = d0 + hours(11) + minutes(44);
if ~isempty(tt_mm.t.TimeZone)
    t1.TimeZone = tt_mm.t.TimeZone;
    t2.TimeZone = tt_mm.t.TimeZone;
end

idx_raw = (tt.t    >= t1) & (tt.t    <= t2);
idx_mm  = (tt_mm.t >= t1) & (tt_mm.t <= t2);

% markers that fall in the zoom window
idx_mark_zoom = (tMark >= t1) & (tMark <= t2);

% ---- figure ----
figure('Color','w');

% ===== Top: full series + green dots on red line =====
subplot(2,1,1);
plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);

plot(tMark, yMark_mm, '.', 'MarkerSize', 18, 'Color', [0 0.6 0], ...
     'DisplayName','Polarimetric Measurement (45s)');

grid on; ylabel('|U| (m/s)');
legend('Raw (Motion Corrected)','Moving Mean Filter', ...
       'Polarimetric Measurement (45s)','Location','best');
% Force x-limits to actual data span (no empty tail)
xlim([tt.t(1) tt.t(end)]);

if exist('Tw','var') && exist('Fs','var')
    title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));
else
    title('Raw vs Moving Mean (full series)');
end

% ===== Bottom: zoom 11:37–11:44 + green dots on red line =====
subplot(2,1,2);
plot(tt.t(idx_raw), Spd_raw(idx_raw), 'k'); hold on;
plot(tt_mm.t(idx_mm), Spd_mm(idx_mm), 'r', 'LineWidth', 1.2);

plot(tMark(idx_mark_zoom), yMark_mm(idx_mark_zoom), '.', ...
     'MarkerSize', 18, 'Color', [0 0.6 0]);

grid on; ylabel('|U| (m/s)'); xlabel('Time');
title('Zoom: 11:37 to 11:44');
xlim([t1 t2]);

% Optional: nicer y-lims in zoom
ymin = min([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
ymax = max([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
if isfinite(ymin) && isfinite(ymax) && ymin < ymax
    ylim([ymin ymax] + 0.05*[-1 1]*(ymax - ymin));
end

%% Adding green tail for duration
clc
% ============================================================
% FULL BLOCK: Raw vs movmean speed + zoom + 45s recolor segments on RED line
% Green dots at folder timestamps ON the movmean line (red), in both panels
%
% Requires in workspace:
%   tt    : timetable with tt.t (datetime) and variables U,V,W (motion-corrected)
%   tt_mm : timetable with tt_mm.t (datetime) and movmean-filtered U,V,W
%   Tw, Fs : optional (for title)
% ============================================================

% ---- speeds ----
Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

% ---- folder timestamps (from folder names) ----
markStr = [
"2025-09-21 10-55-19.320"
"2025-09-21 11-00-48.776"
"2025-09-21 11-09-12.843"
"2025-09-21 11-14-48.914"
"2025-09-21 11-22-14.107"
"2025-09-21 11-25-25.411"
"2025-09-21 11-31-20.321"
"2025-09-21 11-36-06.521"
"2025-09-21 11-39-14.973"
"2025-09-21 11-42-21.662"
"2025-09-21 11-45-26.115"
"2025-09-21 11-48-41.689"
"2025-09-21 11-51-52.353"
"2025-09-21 11-55-26.302"
"2025-09-21 11-58-37.140"
"2025-09-21 12-01-41.780"
];

tMark = datetime(markStr, "InputFormat","yyyy-MM-dd HH-mm-ss.SSS");
if ~isempty(tt_mm.t.TimeZone)
    tMark.TimeZone = tt_mm.t.TimeZone;
end

% keep only markers inside time span
tMark = tMark(tMark >= tt_mm.t(1) & tMark <= tt_mm.t(end));

% y-value at those times on the movmean series (nearest sample)
yMark_mm = interp1(posixtime(tt_mm.t), Spd_mm, posixtime(tMark), "nearest", NaN);

% ---- 45s measurement duration ----
dur_s = 45;

% colors for each 45s segment (edit if you want)
segColors = repmat([0 0.6 0], numel(tMark), 1);

% ---- zoom interval (same day as data) ----
d0 = dateshift(tt_mm.t(1), 'start', 'day');
t1 = d0 + hours(11) + minutes(37);
t2 = d0 + hours(11) + minutes(44);
if ~isempty(tt_mm.t.TimeZone)
    t1.TimeZone = tt_mm.t.TimeZone;
    t2.TimeZone = tt_mm.t.TimeZone;
end

idx_raw = (tt.t    >= t1) & (tt.t    <= t2);
idx_mm  = (tt_mm.t >= t1) & (tt_mm.t <= t2);

% ============================================================
% PLOT
% ============================================================
figure('Color','w');

% ===================== TOP: full series =====================
subplot(2,1,1);

plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);

% green dots at timestamps ON the movmean series
plot(tMark, yMark_mm, '.', 'MarkerSize', 18, 'Color', [0 0.6 0], ...
     'DisplayName','Polarimetric Measurement (45s)');

% recolor red line for each 45s window (overlay segments)
for i = 1:numel(tMark)
    ts = tMark(i);
    te = tMark(i) + seconds(dur_s);
    idxSeg = (tt_mm.t >= ts) & (tt_mm.t <= te) & isfinite(Spd_mm);
    if any(idxSeg)
        plot(tt_mm.t(idxSeg), Spd_mm(idxSeg), '-', ...
            'LineWidth', 3.0, ...
            'Color', segColors(i,:), ...
            'HandleVisibility','off');
    end
end

grid on; ylabel('|U| (m/s)');
legend('Raw (Motion Corrected)','Moving Mean Filter','Polarimetric Measurement (45s)', ...
       'Location','best');

if exist('Tw','var') && exist('Fs','var')
    title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));
else
    title('Raw vs Moving Mean (full series)');
end

% IMPORTANT: avoid empty tail beyond data
xlim([tt_mm.t(1) tt_mm.t(end)]);

% ===================== BOTTOM: zoom 11:37–11:44 =====================
subplot(2,1,2);

plot(tt.t(idx_raw), Spd_raw(idx_raw), 'k'); hold on;
plot(tt_mm.t(idx_mm), Spd_mm(idx_mm), 'r', 'LineWidth', 1.2);

% green dots that fall in zoom
idxMarkZoom = (tMark >= t1) & (tMark <= t2);
plot(tMark(idxMarkZoom), yMark_mm(idxMarkZoom), '.', ...
     'MarkerSize', 18, 'Color', [0 0.6 0], 'HandleVisibility','off');

% recolor segments that overlap zoom
for i = 1:numel(tMark)
    ts = tMark(i);
    te = tMark(i) + seconds(dur_s);
    % segment indices, but only plot if it overlaps zoom window
    if te < t1 || ts > t2
        continue;
    end
    idxSeg = (tt_mm.t >= ts) & (tt_mm.t <= te) & isfinite(Spd_mm);
    if any(idxSeg)
        plot(tt_mm.t(idxSeg), Spd_mm(idxSeg), '-', ...
            'LineWidth', 3.0, ...
            'Color', segColors(i,:), ...
            'HandleVisibility','off');
    end
end

grid on; ylabel('|U| (m/s)'); xlabel('Time');
title('Zoom: 11:37 to 11:44');
xlim([t1 t2]);

% Optional: nicer y-lims in zoom
ymin = min([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
ymax = max([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
if isfinite(ymin) && isfinite(ymax) && ymin < ymax
    ylim([ymin ymax] + 0.05*[-1 1]*(ymax - ymin));
end

%% Green shaded area
clc

%% Adding green rectangle for each 45s measurement window
clc
% ============================================================
% FULL BLOCK: Raw vs movmean speed + zoom
% Green semi-transparent RECTANGLE over each 45s measurement window
%
% Requires in workspace:
%   tt    : timetable with tt.t (datetime) and variables U,V,W (motion-corrected)
%   tt_mm : timetable with tt_mm.t (datetime) and movmean-filtered U,V,W
%   Tw, Fs : optional (for title)
% ============================================================

% ---- speeds ----
Spd_raw = sqrt(double(tt.U).^2    + double(tt.V).^2    + double(tt.W).^2);
Spd_mm  = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

% ---- folder timestamps (from folder names) ----
markStr = [
"2025-09-21 10-55-19.320"
"2025-09-21 11-00-48.776"
"2025-09-21 11-09-12.843"
"2025-09-21 11-14-48.914"
"2025-09-21 11-22-14.107"
"2025-09-21 11-25-25.411"
"2025-09-21 11-31-20.321"
"2025-09-21 11-36-06.521"
"2025-09-21 11-39-14.973"
"2025-09-21 11-42-21.662"
"2025-09-21 11-45-26.115"
"2025-09-21 11-48-41.689"
"2025-09-21 11-51-52.353"
"2025-09-21 11-55-26.302"
"2025-09-21 11-58-37.140"
"2025-09-21 12-01-41.780"
];

tMark = datetime(markStr, "InputFormat","yyyy-MM-dd HH-mm-ss.SSS");
if ~isempty(tt_mm.t.TimeZone)
    tMark.TimeZone = tt_mm.t.TimeZone;
end

% keep only markers inside time span
tMark = tMark(tMark >= tt_mm.t(1) & tMark <= tt_mm.t(end));

% ---- 45s measurement duration ----
dur_s = 45;

% ---- zoom interval (same day as data) ----
d0 = dateshift(tt_mm.t(1), 'start', 'day');
t1 = d0 + hours(11) + minutes(37);
t2 = d0 + hours(11) + minutes(44);
if ~isempty(tt_mm.t.TimeZone)
    t1.TimeZone = tt_mm.t.TimeZone;
    t2.TimeZone = tt_mm.t.TimeZone;
end

idx_raw = (tt.t    >= t1) & (tt.t    <= t2);
idx_mm  = (tt_mm.t >= t1) & (tt_mm.t <= t2);

% ============================================================
% PLOT
% ============================================================
figure('Color','w');

% ===================== TOP: full series =====================
subplot(2,1,1);

plot(tt.t, Spd_raw, 'k'); hold on;
plot(tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2);

% --- green rectangles (full panel) ---
yl = ylim;  % after plotting, so limits exist
for i = 1:numel(tMark)
    ts = tMark(i);
    te = tMark(i) + seconds(dur_s);

    % clip to data range (optional but clean)
    ts2 = max(ts, tt_mm.t(1));
    te2 = min(te, tt_mm.t(end));
    if te2 <= ts2, continue; end

    patch([ts2 te2 te2 ts2], [yl(1) yl(1) yl(2) yl(2)], [0 0.6 0], ...
        'FaceAlpha', 0.12, 'EdgeColor', [0 0.6 0], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
end
% bring curves to front
uistack(findobj(gca,'Type','line'),'top');

grid on; ylabel('|U| (m/s)');
legend('Raw (Motion Corrected)','Moving Mean Filter', 'Location','best');

if exist('Tw','var') && exist('Fs','var')
    title(sprintf('Moving mean window Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));
else
    title('Raw vs Moving Mean (full series)');
end

% IMPORTANT: avoid empty tail beyond data
xlim([tt_mm.t(1) tt_mm.t(end)]);

% ===================== BOTTOM: zoom 11:37–11:44 =====================
subplot(2,1,2);

plot(tt.t(idx_raw), Spd_raw(idx_raw), 'k'); hold on;
plot(tt_mm.t(idx_mm), Spd_mm(idx_mm), 'r', 'LineWidth', 1.2);

% --- green rectangles (zoom panel) ---
yl = ylim;
for i = 1:numel(tMark)
    ts = tMark(i);
    te = tMark(i) + seconds(dur_s);

    % only draw if overlaps zoom
    if te < t1 || ts > t2
        continue;
    end

    ts2 = max(ts, t1);
    te2 = min(te, t2);
    if te2 <= ts2, continue; end

    patch([ts2 te2 te2 ts2], [yl(1) yl(1) yl(2) yl(2)], [0 0.6 0], ...
        'FaceAlpha', 0.12, 'EdgeColor', [0 0.6 0], 'LineWidth', 0.8, ...
        'HandleVisibility','off');
end
uistack(findobj(gca,'Type','line'),'top');

grid on; ylabel('|U| (m/s)'); xlabel('Time');
title('Zoom: 11:37 to 11:44');
xlim([t1 t2]);

% Optional: nicer y-lims in zoom (recompute after xlim)
ymin = min([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
ymax = max([Spd_raw(idx_raw); Spd_mm(idx_mm)], [], 'omitnan');
if isfinite(ymin) && isfinite(ymax) && ymin < ymax
    ylim([ymin ymax] + 0.05*[-1 1]*(ymax - ymin));
end


%% only smooth + 4 measurements
clc
% ============================================================
% ONE PANEL: movmean speed only (red) + 45s overlay segments
% Different color for EACH point + its 45s tail
%
% Points:
%  11-22-14.107, 11-36-06.521, 11-45-26.115, 11-58-37.140, 12-01-41.780
% ============================================================

Spd_mm = sqrt(double(tt_mm.U).^2 + double(tt_mm.V).^2 + double(tt_mm.W).^2);

d0 = dateshift(tt_mm.t(1), 'start', 'day');

tMark = d0 + [ ...
    %hours(11) + minutes(22) + seconds(14.107)
    hours(11) + minutes(36) + seconds( 6.521)
    hours(11) + minutes(45) + seconds(26.115)
    %hours(11) + minutes(58) + seconds(37.140)
    hours(12) + minutes( 1) + seconds(41.780)
];

if ~isempty(tt_mm.t.TimeZone)
    tMark.TimeZone = tt_mm.t.TimeZone;
end

% keep only markers inside data span
tMark = tMark(tMark >= tt_mm.t(1) & tMark <= tt_mm.t(end));

% y at marker times (nearest sample on movmean)
yMark_mm = interp1(posixtime(tt_mm.t), Spd_mm, posixtime(tMark), "nearest", NaN);

dur_s = 45;

% ===== choose a different color per window (Nx3 RGB in 0..1) =====
% Edit these if you want other colors:
segColors = [ ...
    0.00 0.60 0.00   % green
    0.00 0.45 0.70   % orange
    0.49 0.18 0.56   % purple
    0.00 0.45 0.70   % blue
    0.93 0.69 0.13   % yellow
];

% If some tMarks got clipped (outside data), clip colors too
segColors = segColors(1:numel(tMark),:);

figure('Color','w');
ax = axes; hold(ax,'on'); grid(ax,'on');

% ---- fonts ----
ax.FontName = 'Arial';   % ou 'Helvetica'
ax.FontSize = 18;        % ticks
ylabel(ax,'U (m/s)','FontSize',22);
xlabel(ax,'Time','FontSize',22);
title(ax, ax.Title.String, 'FontSize',22);


% base smoothed line
hmm = plot(ax, tt_mm.t, Spd_mm, 'r', 'LineWidth', 1.2, ...
    'DisplayName','Moving Mean Filter');

% points + 45s tails, each with its own color
hpts = gobjects(numel(tMark),1);

for i = 1:numel(tMark)
    ts = tMark(i);
    te = ts + seconds(dur_s);

    % point at start time
    hpts(i) = plot(ax, ts, yMark_mm(i), '.', ...
        'MarkerSize', 22, 'Color', segColors(i,:), ...
        'DisplayName', sprintf('45s window %d', i));

    % overlay tail segment on the red line
    idxSeg = (tt_mm.t >= ts) & (tt_mm.t <= te) & isfinite(Spd_mm);
    if any(idxSeg)
        plot(ax, tt_mm.t(idxSeg), Spd_mm(idxSeg), '-', ...
            'LineWidth', 3.2, ...
            'Color', segColors(i,:), ...
            'HandleVisibility','off');  % keep legend clean
    end
end
pbaspect(ax,[19 6 1]);   % large et bas

ylabel(ax,'$|U_{5}|$ (m/s)','Interpreter','latex');
xlabel(ax,'Time');

if exist('Tw','var') && exist('Fs','var')
    title(ax, sprintf('Smoothed Wind Magnitude Tw=%.1f s (Fs=%.2f Hz)', Tw, Fs));
else
    title(ax,'Moving mean speed with 45s colored windows');
end

xlim(ax, [tt_mm.t(1) tt_mm.t(end)]);

%% Mean wind + qualitative stability (with movmean) 
% Goal: clean mean wind time series + regime classification (stable/unstable)
% NOT for flux-grade u*, L, z0 magnitudes.

% REQUIRE:
%   tt : timetable with tt.t (datetime) and variables U,V,W,T (already motion-corrected)
%   process_sonic_blocks_from_tt.m on path

% --- user settings ---
Tw       = 60;      % moving mean window [s]
blockMin = 10;      % block size [min] (Inf = whole series)
z_m      = 5.0;     % sonic height [m]  <-- EDIT

% --- dt/Fs ---
dtv = seconds(diff(tt.t)); dtv = dtv(~isnan(dtv) & dtv>0);
dt  = median(dtv);
Fs  = 1/dt;

% --- movmean window in samples ---
Nw = max(5, round(Tw/dt));

% --- smooth components (keep timestamps) ---
ttA = tt;
ttA.U = movmean(double(tt.U), Nw, 'omitnan');
ttA.V = movmean(double(tt.V), Nw, 'omitnan');
ttA.W = movmean(double(tt.W), Nw, 'omitnan');
ttA.T = movmean(double(tt.T), Nw, 'omitnan');

% --- run MOST blocks (used only for qualitative stability) ---
outA = process_sonic_blocks_from_tt(ttA, blockMin, z_m);
B    = outA.blocks;

% --- mean wind diagnostics from smoothed components ---
Umean = ttA.U;  Vmean = ttA.V;  Wmean = ttA.W;

Spd_mean = sqrt(Umean.^2 + Vmean.^2 + Wmean.^2);

% Meteorological wind direction FROM (deg), with convention:
% N = U, W = V => East = -West
Emean = -Vmean;
WD_from = mod(270 - atan2d(Emean, Umean), 360);

% --- qualitative stability label from L (Monin–Obukhov length) ---
% Convention: L>0 stable, L<0 unstable, |L| very large ~ neutral
L = B.L;
stab = strings(height(B),1);
for i=1:height(B)
    if ~isfinite(L(i)) || abs(L(i)) > 500
        stab(i) = "neutral/weak";
    elseif L(i) > 0
        stab(i) = "stable";
    else
        stab(i) = "unstable";
    end
end

% --- display a compact summary table (what Option A reports) ---
summaryA = table(B.time, B.U_mean, B.u_star, B.L, stab, ...
    'VariableNames', {'time','U_mean','u_star_diagnostic','L_diagnostic','stability'});
disp(summaryA);

% --- plots: mean wind magnitude + direction + stability over blocks ---
figure('Color','w','Name','Option A: mean wind + qualitative stability');

subplot(3,1,1);
plot(ttA.t, Spd_mean, 'k'); grid on;
ylabel('|U|_{mean} (m/s)');
title(sprintf('Movmean Tw=%.1f s (Nw=%d, Fs=%.2f Hz)', Tw, Nw, Fs));

subplot(3,1,2);
plot(ttA.t, WD_from, 'k'); grid on;
ylabel('WD_{from} (deg)'); ylim([0 360]);

subplot(3,1,3);
plot(B.time, B.L, 'o-'); grid on;
ylabel('L_{MO} (m)'); xlabel('Time');
yline(0,'--');

% annotate stability labels
for i=1:height(B)
    text(B.time(i), B.L(i), "  "+stab(i), 'FontSize', 9);
end

%% ===== Center WD around a reference (no 0/360 wrapping spikes) =====
% Assumes ttA has ttA.t, ttA.U (North), ttA.V (West). East = -West.

U = double(ttA.U(:));        % North
W = double(ttA.V(:));        % West
E = -W;                      % East

% Meteorological direction FROM in degrees (0..360)
WD = mod(270 - atan2d(E, U), 360);

% Choose a center/reference direction (deg FROM)
theta = WD*pi/180;
C = mean(cos(theta),'omitnan');
S = mean(sin(theta),'omitnan');
WDref = mod(atan2d(S,C), 360);

% Compute signed difference in [-180, 180] without toolboxes
d = WD - WDref;
d = mod(d + 180, 360) - 180;     % now in [-180, 180]

% Centered direction (lives around WDref, can go below 0 or above 360)
WD_centered = WDref + d;

% Plot
figure('Color','w');
plot(ttA.t, WD_centered, 'k'); grid on;
ylabel(sprintf('WD centered around %g° (deg)', WDref));
xlabel('Time');
title('Wind direction centered (no wrap spikes)');

figure('Color','w');
plot(outA.blocks.U_mean, outA.blocks.u_star, 'o');
grid on;
xlabel('Mean wind speed (m/s)');
ylabel('u_* (m/s)');
title('Block u_* vs mean wind speed');

scatter(outA.blocks.U_mean, outA.blocks.u_star, 80, sign(outA.blocks.L), 'filled');
colormap([0 0 1; 1 0 0]);  % stable vs unstable
colorbar; caxis([-1 1]);
xlabel('Mean wind speed (m/s)'); ylabel('u_* (m/s)');


