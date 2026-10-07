function res = allan_analysis(data, dt, varargin)
%ALLAN_ANALYSIS  Overlapping Allan deviation + Bias Instability (BI) check.
%   No toolbox required.
%
%   res = allan_analysis(data, dt)
%   res = allan_analysis(data, dt, 'Detrend','linear','Units','A', ...)
%
%   data : vector, uniformly sampled output (e.g. photocurrents)
%   dt   : sample interval (s)
%
%   Name-Value options
%     'Type'              'rate' (default) -> data is a rate, integrated internally
%                         'angle'          -> data is already integrated (phase/angle)
%     'Scale'             multiply data by this factor first (e.g. K_scale to deg/h)
%     'Detrend'           'none' | 'mean' (default) | 'linear' | 'poly2'
%                         Remove deterministic signal/thermal drift before ADEV.
%     'Units'             label for plots (unit of data after Scale)
%     'NumPoints'         number of log-spaced cluster sizes (default 150)
%     'MinClusters'       longest tau = record length / MinClusters (default 9)
%     'SlopeTol'          |log-log slope| below this counts as "flat" (default 0.15)
%     'MinPlateauDecades' flat region must span at least this many decades (default 0.3)
%     'Plot'              true (default) | false
%
%   Output struct res: tau, adev, ci_lo, ci_hi, slope, BI, tau_BI, plateau_decades,
%   hasMin, hasPlateau, hasPeak, tau_peak, tauc_est, N_arw, K_rrw, R_rr, verdict.

p = inputParser;
addParameter(p,'Type','rate',@(s)any(strcmpi(s,{'rate','angle'})));
addParameter(p,'Scale',1,@isscalar);
addParameter(p,'Detrend','mean',@(s)any(strcmpi(s,{'none','mean','linear','poly2'})));
addParameter(p,'Units','a.u.',@(s)ischar(s)||isstring(s));
addParameter(p,'NumPoints',150,@isscalar);
addParameter(p,'MinClusters',9,@isscalar);
addParameter(p,'SlopeTol',0.15,@isscalar);
addParameter(p,'MinPlateauDecades',0.3,@isscalar);
addParameter(p,'Plot',true,@islogical);
parse(p,varargin{:});
o = p.Results;
units = char(o.Units);

%% 1. Preprocess
y = double(data(:)) * o.Scale;
if any(~isfinite(y)), error('data contains NaN/Inf'); end
N = numel(y);
if N < 100, error('record too short'); end

switch lower(o.Detrend)
    case 'mean',   y = y - mean(y);
    case 'linear', x = linspace(-1,1,N).'; y = y - polyval(polyfit(x,y,1),x);
    case 'poly2',  x = linspace(-1,1,N).'; y = y - polyval(polyfit(x,y,2),x);
end

if strcmpi(o.Type,'rate')
    theta = [0; cumsum(y)*dt];
else
    theta = y;
end
Nth = numel(theta);

%% 2. Overlapping Allan deviation
mMax = floor(N / o.MinClusters);
if mMax < 4, error('record too short for MinClusters = %d', o.MinClusters); end
m = unique(round(logspace(0, log10(mMax), o.NumPoints))).';
m = m(2*m < Nth);
M = numel(m);
if M < 8, error('too few tau points; increase NumPoints or record length'); end

tau  = m * dt;
adev = zeros(M,1);
for i = 1:M
    mi = m(i);
    d = theta(1+2*mi:Nth) - 2*theta(1+mi:Nth-mi) + theta(1:Nth-2*mi);
    adev(i) = sqrt(sum(d.^2) / (2*(mi*dt)^2 * numel(d)));
end

% Approximate 1-sigma confidence band (relative error ~ 1/sqrt(2(N/m-1)))
relErr = 1 ./ sqrt(2*(N./m - 1));
ci_lo = adev .* (1 - relErr);
ci_hi = adev .* (1 + relErr);

%% 3. Log-log slope
ls   = log10(tau);
la_s = movmean(log10(adev), 5);
slope = movmean(gradient(la_s, ls), 5);
idx = (1:M).';

%% 4. Minimum / plateau (classical BI)
[~, k] = min(la_s);
hasMin = (k >= 3) && (k <= M-2) && (la_s(1) - la_s(k) > 0.04);
l = k; r = k;
while l > 1 && abs(slope(l-1)) <= o.SlopeTol, l = l-1; end
while r < M && abs(slope(r+1)) <= o.SlopeTol, r = r+1; end
plateauDec = ls(r) - ls(l);
hasPlateau = hasMin && plateauDec >= o.MinPlateauDecades;

BI = NaN; tauBI = NaN;
if hasMin
    BI = 10^la_s(k) / sqrt(2*log(2)/pi);   % sigma_min = 0.664 * B
    tauBI = tau(k);
end

%% 5. Peak (correlated-noise bump, e.g. first-order Gauss-Markov)
pk = find(slope(1:end-1) > 0 & slope(2:end) <= 0);
pk = pk(pk >= 3 & pk <= M-2);
keep = false(size(pk));
for q = 1:numel(pk)
    kk = pk(q);
    keep(q) = (la_s(kk) - min(la_s(1:kk)) > 0.06) && (la_s(kk) - min(la_s(kk:M)) > 0.06);
end
pk = pk(keep);
hasPeak = ~isempty(pk);
kp = NaN; tau_peak = NaN; tauc_est = NaN;
if hasPeak
    [~, j] = max(la_s(pk));
    kp = pk(j);
    tau_peak = tau(kp);
    tauc_est = tau_peak / 1.89;   % first-order Markov: peak near 1.89*Tc (IEEE 952)
end

%% 6. Other noise coefficients (indicative, from slope windows)
after = idx > k*hasMin;
firstRise = find(slope > 0.10, 1);
if isempty(firstRise), firstRise = M+1; end
selN = abs(slope + 0.5) < 0.10 & idx < firstRise;   % head region only (before any rise)
selK = abs(slope - 0.5) < 0.10 & after;
selR = abs(slope - 1.0) < 0.15 & after;
N_arw = NaN; K_rrw = NaN; R_rr = NaN;
if any(selN), N_arw = median(adev(selN) .* sqrt(tau(selN))); end
if any(selK), K_rrw = median(adev(selK) .* sqrt(3 ./ tau(selK))); end
if any(selR), R_rr  = median(adev(selR) .* sqrt(2) ./ tau(selR)); end

%% 7. Verdict
nq = max(1, round(M/4));
headSlope = mean(slope(1:nq));
tailSlope = mean(slope(M-nq+1:M));
if hasPlateau
    verdict = sprintf('BI PRESENT: flat minimum, %.2f decades wide', plateauDec);
elseif hasMin
    verdict = sprintf(['BI POSSIBLE: minimum exists but flat part is only %.2f decades ' ...
        '(masked by RRW/RR or record too short)'], plateauDec);
elseif hasPeak
    verdict = 'NO classical BI minimum; correlated-noise bump (peak) found';
elseif tailSlope < -0.3
    verdict = 'NO BI visible: white noise still dominates at longest tau (need longer record)';
elseif tailSlope > 0.3
    verdict = 'NO BI visible: drift/random-walk dominates, no minimum';
else
    verdict = 'NO clear BI: no interior minimum';
end

fprintf('\n===== Allan deviation summary =====\n');
fprintf('Samples: %d | dt = %g s | tau range: %g .. %g s\n', N, dt, tau(1), tau(end));
fprintf('Detrend: %s | Type: %s | Units: %s\n', o.Detrend, o.Type, units);
fprintf('Verdict: %s\n', verdict);
if hasMin
    fprintf('  min ADEV = %.4g %s at tau = %.4g s -> BI (min/0.664) = %.4g %s\n', ...
        10^la_s(k), units, tauBI, BI, units);
end
if hasPeak
    fprintf('  peak at tau = %.4g s -> estimated correlation time Tc ~ %.4g s\n', tau_peak, tauc_est);
end
if ~isnan(N_arw), fprintf('  ARW-like (slope -1/2): N = %.4g %s/sqrt(s)\n', N_arw, units); end
if ~isnan(K_rrw), fprintf('  RRW-like (slope +1/2): K = %.4g\n', K_rrw); end
if ~isnan(R_rr),  fprintf('  Rate-ramp-like (slope +1): R = %.4g\n', R_rr); end
fprintf('====================================\n');

%% 8. Plot
if o.Plot
    figure('Name','Allan deviation','Color','w');
    subplot(2,1,1);
    ax = gca; set(ax,'XScale','log','YScale','log'); hold(ax,'on');
    hl = []; ln = {};
    hb = fill([tau; flipud(tau)], [ci_lo; flipud(ci_hi)], [0.80 0.88 1], ...
        'EdgeColor','none','FaceAlpha',0.6);
    hl(end+1) = hb; ln{end+1} = '~1\sigma band';
    hl(end+1) = plot(tau, adev, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2); ln{end+1} = 'ADEV';
    if hasPlateau
        hl(end+1) = plot(tau(l:r), adev(l:r), 'g-', 'LineWidth', 3); ln{end+1} = 'flat region';
    end
    if hasMin
        hl(end+1) = plot(tau(k), 10^la_s(k), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        ln{end+1} = 'minimum';
        text(tau(k), 10^la_s(k)*0.7, sprintf('  min = %.3g %s', 10^la_s(k), units), ...
            'Color','r','VerticalAlignment','top');
    end
    if hasPeak
        hl(end+1) = plot(tau(kp), adev(kp), 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
        ln{end+1} = 'peak';
    end
    grid on; box on;
    xlabel('\tau (s)'); ylabel(['\sigma(\tau) (' units ')']);
    title(verdict, 'FontSize', 10);
    legend(hl, ln, 'Location','best');

    subplot(2,1,2);
    semilogx(tau, slope, 'k-', 'LineWidth', 1.2); hold on;
    plot(tau([1 end]), [ o.SlopeTol  o.SlopeTol], 'g--');
    plot(tau([1 end]), [-o.SlopeTol -o.SlopeTol], 'g--');
    plot(tau([1 end]), [-0.5 -0.5], ':', 'Color', [0.5 0.5 0.5]);
    plot(tau([1 end]), [ 0.5  0.5], ':', 'Color', [0.5 0.5 0.5]);
    plot(tau([1 end]), [ 1.0  1.0], ':', 'Color', [0.5 0.5 0.5]);
    grid on; ylim([-1.5 1.5]);
    xlabel('\tau (s)'); ylabel('d log\sigma / d log\tau');
    title('Local slope: -1/2 white, 0 flicker (BI), +1/2 RRW, +1 rate ramp');
end

%% Output
res = struct('tau',tau,'adev',adev,'ci_lo',ci_lo,'ci_hi',ci_hi,'slope',slope, ...
    'BI',BI,'tau_BI',tauBI,'plateau_decades',plateauDec, ...
    'hasMin',hasMin,'hasPlateau',hasPlateau,'hasPeak',hasPeak, ...
    'tau_peak',tau_peak,'tauc_est',tauc_est, ...
    'N_arw',N_arw,'K_rrw',K_rrw,'R_rr',R_rr,'verdict',verdict);
end
