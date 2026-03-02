function out = process_sonic_log_blocks(logFilePath, blockMin, z_m)
% CHristopher Bouillon 2025
% PROCESS_SONIC_LOG_BLOCKS
% Parse sonic CSV and compute u*, L, and z0 (MOST + Businger–Dyer) in blocks
% of a user-specified length.
%
% Usage:
%   % whole file as one block:
%   out = process_sonic_log_blocks('anemometer_log.csv', Inf, 3.5);
%
%   % 1-minute blocks:
%   out = process_sonic_log_blocks('anemometer_log.csv', 1, 3.5);
%
%   % 5-minute blocks:
%   out = process_sonic_log_blocks('anemometer_log.csv', 5, 3.5);
%
% Inputs
%   logFilePath : path to the CSV log
%   blockMin    : block length in minutes (use Inf for whole file)
%   z_m         : sonic height above mean water level [m]
%
% Output struct "out" contains:
%   out.tt       : timetable of parsed U,V,W,T
%   out.blocks   : table with time, u_star, U_mean, H, L, z0 per block
%   out.Fs       : estimated sampling rate [Hz]
%   out.figures  : figure handles
%
% Notes
% - Expects lines like:
%   2025-09-13 12:21:35.169, +000.62;-000.04;-000.63;+18.3;0E;41, NaN,-0.04,-0.63,+18.3,0E
%   First 4 numeric tokens after timestamp are interpreted as U,V,W,T.
% - Temperature T may be in °C; we auto-handle K conversion **only** for L.
% - Stability function: Businger–Dyer.

% --- constants
kappa = 0.40; g = 9.81; rho = 1.225; cp = 1004;

% --- read file as raw text
raw = fileread(logFilePath);
lines = regexp(raw, '\r\n|\n|\r', 'split');
lines = lines(~cellfun(@isempty,lines));

% --- parse lines
ts_ok = '^\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}(\.\d+)?$';
t = NaT(numel(lines),1,'Format','yyyy-MM-dd HH:mm:ss.SSS');
U = nan(size(t)); V = nan(size(t)); W = nan(size(t)); T = nan(size(t));

k = 0;
for i = 1:numel(lines)
    L = strtrim(lines{i});
    cidx = find(L==',',1,'first'); 
    if isempty(cidx), continue; end

    tstr = strtrim(L(1:cidx-1));
    if isempty(regexp(tstr, ts_ok, 'once'))
        % skip headers or malformed stamps
        continue
    end

    % parse timestamp (ms first, then seconds)
    try
        tparsed = datetime(tstr,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS','Locale','en_US_POSIX');
    catch
        tparsed = datetime(tstr,'InputFormat','yyyy-MM-dd HH:mm:ss','Locale','en_US_POSIX');
    end

    % strip control chars, then grab first 4 numeric tokens (U,V,W,T)
    payload = regexprep(L(cidx+1:end),'[\x00-\x1F\x7F]','');
    toks = regexp(payload,'([+\-]?\d+(?:\..\d+)?|NaN|nan)','tokens');
    if numel(toks) < 4
        % last resort: scan whole line
        toks = regexp(L,'([+\-]?\d+(?:\.\d+)?|NaN|nan)','tokens');
        if numel(toks) < 4, continue; end
    end

    vals = nan(1,4);
    for j = 1:4
        s = toks{j}{1};
        if any(strcmpi(s,{'NaN','nan'})), vals(j) = NaN; else, vals(j) = str2double(s); end
    end

    k = k + 1;
    t(k) = tparsed; U(k) = vals(1); V(k) = vals(2); W(k) = vals(3); T(k) = vals(4);
end

% trim to parsed rows & sort by time
t = t(1:k); U = U(1:k); V = V(1:k); W = W(1:k); T = T(1:k);
keep = ~isnat(t);
t = t(keep); U=U(keep); V=V(keep); W=W(keep); T=T(keep);
[ t, sortIdx ] = sort(t); U = U(sortIdx); V = V(sortIdx); W = W(sortIdx); T = T(sortIdx);

if isempty(t)
    error('No parsable data rows found in %s.', logFilePath);
end

tt = timetable(t, U, V, W, T);

% --- estimate sampling rate
dt = seconds(diff(tt.t));
dt = dt(~isnan(dt) & dt > 0);
if isempty(dt)
    Fs = 20;  % fallback
else
    Fs = round(1/median(dt));
end
fprintf('Parsed %d samples. Estimated Fs ≈ %d Hz\n', height(tt), Fs);

% --- mild cleaning (optional but helpful)
win = max(5, round(Fs*1)); % ~1 s movmedian window
tt.U = filloutliers(tt.U,'linear','movmedian',win);
tt.V = filloutliers(tt.V,'linear','movmedian',win);
tt.W = filloutliers(tt.W,'linear','movmedian',win);
tt.T = filloutliers(tt.T,'linear','movmedian',win);

% --- determine blocks
N = height(tt);
if isinf(blockMin)
    Nwin = N; 
    nBlocks = 1;
else
    Nwin = max(10, round(blockMin * 60 * Fs));
    nBlocks = floor(N / Nwin);
    if nBlocks < 1
        % keep at least one block using all data
        warning('File shorter than requested block; using whole file as 1 block.');
        Nwin = N; nBlocks = 1;
    end
end

% --- preallocate outputs
blkTime = NaT(nBlocks,1);
u_star  = nan(nBlocks,1);
H       = nan(nBlocks,1);
Lmo     = nan(nBlocks,1);
z0      = nan(nBlocks,1);
Uz_bar  = nan(nBlocks,1);

% --- per-block calculations
for b = 1:nBlocks
    i0 = (b-1)*Nwin + 1;
    i1 = min(b*Nwin, N);

    % means for rotations (use raw-ish series) 
    % **************************************
    % **************************************
    % **************************************
    % **************************************

    % NEED TO ADD THE INERTIAL STATION VALUES, BECAUSE BOAT MOTION 
    % LEAKS + SPEED LEAKS INTO THE WIND
    % PARAMETERS AND CONTAMINATES THE TRUE FLUXES
    
    % **************************************
    % **************************************
    % **************************************
    % **************************************

    ubar = mean(tt.U(i0:i1),'omitnan');
    vbar = mean(tt.V(i0:i1),'omitnan');
    Tbar = mean(tt.T(i0:i1),'omitnan');

    % rotate to mean-wind coordinates
%The mean wind direction becomes the +x direction.
%The cross-wind component becomes the new y direction.
%Vertical W is unchanged.
% Because turbulent fluxes rely on moments like 
% mean(𝑢′𝑤′). To avoid energy from u to go into v
% Aligns the coordinate axes so that U is parallel to the mean wind.

    alpha = atan2(vbar, ubar);              % yaw
    U1 =  tt.U(i0:i1)*cos(alpha) + tt.V(i0:i1)*sin(alpha);
    V1 = -tt.U(i0:i1)*sin(alpha) + tt.V(i0:i1)*cos(alpha);
    W1 =  tt.W(i0:i1);

% Because the sonic measures velocity components in its own tilted reference frame.
% If the instrument has even a 2–3° tilt: energy from x gows into w
    beta = atan2(mean(W1,'omitnan'), mean(U1,'omitnan'));  % pitch

    U2 =  U1*cos(beta) + W1*sin(beta);
    V2 =  V1;
    W2 = -U1*sin(beta) + W1*cos(beta);

    % fluctuations (remove block means)
    up = U2 - mean(U2,'omitnan');
    vp = V2 - mean(V2,'omitnan');
    wp = W2 - mean(W2,'omitnan');
    Tp = tt.T(i0:i1) - mean(tt.T(i0:i1),'omitnan');

    % covariances and fluxes
    uw = mean(up.*wp,'omitnan');     % u'w'
    vw = mean(vp.*wp,'omitnan');     % v'w'
    wT = mean(wp.*Tp,'omitnan');     % w'T'

    ustar = ((-uw)^2 + (-vw)^2)^(1/4);
    TmeanK = Tbar + (Tbar < 200)*273.15; % if T looked like °C, use K for L

    if isfinite(wT) && wT ~= 0
        L = -(ustar^3 * TmeanK) / (kappa * g * wT);
    else
        L = NaN;
    end

    Uz = sqrt(ubar^2 + vbar^2);

    if isfinite(L) && L ~= 0
        psi_m = psi_m_businger(z_m / L);
    else
        psi_m = 0;
    end

    z0_b = z_m * exp( -(kappa * Uz / max(ustar,1e-6)) - psi_m );

    % store
    blkTime(b) = tt.t(i0) + seconds((i1 - i0) / (2*Fs));
    u_star(b)  = ustar;
    H(b)       = rho * cp * wT;
    Lmo(b)     = L;
    z0(b)      = z0_b;
    Uz_bar(b)  = Uz;
end

% --- results table
blocks = table(blkTime, u_star, Uz_bar, H, Lmo, z0, ...
    'VariableNames', {'time','u_star','U_mean','H','L','z0'});

% --- plots
fig1 = figure('Name','U,V,W,T time series','Color','w');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
nexttile; plot(tt.t,tt.U); ylabel('U [m s^{-1}]'); grid on;
nexttile; plot(tt.t,tt.V); ylabel('V [m s^{-1}]'); grid on;
nexttile; plot(tt.t,tt.W); ylabel('W [m s^{-1}]'); grid on;
nexttile; plot(tt.t,tt.T); ylabel('T [°C or K]'); xlabel('Time'); grid on;

fig2 = figure('Name','u*, L, z0 (per block)','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile; plot(blocks.time,blocks.u_star,'-o'); ylabel('u_* [m s^{-1}]'); grid on;
nexttile; plot(blocks.time,blocks.L,'-o'); yline(0,'k-'); ylabel('L [m]'); grid on;
nexttile; semilogy(blocks.time,blocks.z0,'-o'); ylabel('z_0 [m]'); xlabel('Time'); grid on;

% --- save summary CSV next to the log
outCsv = fullfile(fileparts(logFilePath),'z0_blocks.csv');
writetable(blocks, outCsv);
fprintf('Saved block summary to %s\n', outCsv);

% --- pack outputs
out = struct();
out.tt       = tt;
out.blocks   = blocks;
out.Fs       = Fs;
out.figures  = [fig1 fig2];

end

% ===== helper =====
function ps = psi_m_businger(zL)
% Businger–Dyer stability correction for momentum - FROM ROMANIC 2020
% (NETHERLAND STORM)
ps = zeros(size(zL));
unst = zL < 0;
if any(unst)
    x = (1 - 16*zL(unst)).^(1/4);
    ps(unst) = 2*log((1+x)./2) + log((1+x.^2)./2) - 2*atan(x) + pi/2;
end
stab = zL > 0;
if any(stab)
    ps(stab) = -5*zL(stab);
end
end
