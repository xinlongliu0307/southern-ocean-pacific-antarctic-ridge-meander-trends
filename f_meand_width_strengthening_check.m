function out = f_meand_width_strengthening_check(meand_width, meand_loc_speed, meand_lon)
% Meander width vs along-jet speed: regression, residual trends, visuals
% Inputs:
%   meand_width      [nLon x nT] km
%   meand_loc_speed  [nLon x nT] m s^-1
%   meand_lon        [nLon x 1] degrees (optional)
%
% Output: struct with diagnostics and region-mean series

%% ---- Basic checks
assert(nargin>=2, 'Provide meand_width and meand_loc_speed.');
[nLon,nT] = size(meand_width);
assert(isequal(size(meand_loc_speed),[nLon,nT]), 'Size mismatch between width and speed.');

if nargin>=3 && ~isempty(meand_lon) && numel(meand_lon)==nLon
    lon = meand_lon(:);
else
    lon = (1:nLon)'; % fallback index
end

%% ---- Settings
min_time_fraction = 0.80;   % >=80% valid points per longitude
use_robust_reg    = true;   % robust regression to limit outliers
export_base       = 'width_speed_regression_diagnostics';

%% ---- Time vector (monthly starting 1993-01)
t0   = datetime(1993,1,15);
t    = t0 + calmonths(0:nT-1);
t_dec= years(t - t(1))/10;  % decades since first month

%% ---- Remove monthly climatology (return anomalies)
W = remove_monthly_climatology(meand_width,     t);
U = remove_monthly_climatology(meand_loc_speed, t);

%% ---- Storage
slopeW   = nan(nLon,1);  pW   = nan(nLon,1);   % raw width trend (km/dec), p(MK)
slopeR   = nan(nLon,1);  pR   = nan(nLon,1);   % residual width trend (km/dec), p(MK)
slopeYh  = nan(nLon,1);                          % width trend explained by speed (km/dec)
frac_exp = nan(nLon,1);                          % fraction explained by speed
R2loc    = nan(nLon,1);                          % regression R^2
W_resid  = nan(nLon,nT);                         % residual width time series

%% ---- Regress width on speed per longitude; residuals and trends
for i = 1:nLon
    wi = W(i,:)'; ui = U(i,:)';
    good = isfinite(wi) & isfinite(ui);
    if nnz(good) < min_time_fraction*nT, continue; end

    % z-score predictor on valid times
    u_mu = mean(ui(good)); u_sd = std(ui(good));
    if ~isfinite(u_sd) || u_sd==0, continue; end
    uz = (ui(good)-u_mu)/u_sd;  yi = wi(good);

    % regression (width ~ speed)
    if use_robust_reg && exist('fitlm','file')
        mdl = fitlm(uz, yi, 'Intercept', true, 'RobustOpts','on');
    elseif exist('fitlm','file')
        mdl = fitlm(uz, yi, 'Intercept', true);
    else
        p = polyfit(uz, yi, 1); mdl = []; mdl.Rsquared.Ordinary = corr(yi, polyval(p,uz),'rows','complete')^2;
    end
    R2loc(i) = mdl.Rsquared.Ordinary;

    % residual width on valid times; map back
    ri = nan(size(wi));
    if ~isempty(mdl), ri(good) = yi - predict(mdl, uz); else, ri(good) = yi - polyval(p, uz); end
    W_resid(i,:) = ri';

    % predicted component (due to speed) across full series
    uz_full = (ui - u_mu)/u_sd;
    yhat = nan(size(wi));
    mask = isfinite(uz_full);
    if ~isempty(mdl), yhat(mask) = predict(mdl, uz_full(mask));
    else, yhat(mask) = polyval(p, uz_full(mask)); end

    % Theil–Sen slopes (km/decade) + MK p-values
    [sW, pMkW] = sen_mk(wi,  t_dec);
    [sR, pMkR] = sen_mk(ri,  t_dec);
    [sY, ~    ] = sen_mk(yhat,t_dec);

    slopeW(i)  = sW;   pW(i) = pMkW;
    slopeR(i)  = sR;   pR(i) = pMkR;
    slopeYh(i) = sY;

    % fraction of raw width trend explained by speed (sign-aware)
    if isfinite(sW) && sW~=0 && isfinite(sY)
        frac_exp(i) = sY / sW;
    end
end

%% ---- Region-mean series and trends
Yhat    = W - W_resid;                    % component explained by speed
W_mean  = mean(W,      1, 'omitnan');
WR_mean = mean(W_resid,1, 'omitnan');
Yh_mean = mean(Yhat,   1, 'omitnan');

[sWm, pWm]   = sen_mk(W_mean',  t_dec);
[sWRm,pWRm]  = sen_mk(WR_mean', t_dec);
[sYhm,~]     = sen_mk(Yh_mean', t_dec);

fprintf('--- Regression summary ---\n');
fprintf('Median raw width trend (km/dec):        %.3f\n', median(slopeW,'omitnan'));
fprintf('Median residual width trend (km/dec):   %.3f\n', median(slopeR,'omitnan'));
fprintf('Median fraction explained by speed:     %.2f\n', median(frac_exp,'omitnan'));
fprintf('Region-mean raw width trend (km/dec):   %.3f (p=%.3f)\n', sWm,  pWm);
fprintf('Region-mean resid width trend (km/dec): %.3f (p=%.3f)\n', sWRm, pWRm);

%% ===================== Visualization (JGR: Oceans style) ==================
col_raw = [0.00 0.45 0.74];   % blue
col_hat = [0.85 0.33 0.10];   % orange
col_res = [0.47 0.67 0.19];   % green
col_grey= [0.35 0.35 0.35];

% Trend lines centered at median time
t0c  = median(t_dec,'omitnan');
W_fit  = median(W_mean,'omitnan')  + sWm  *(t_dec - t0c);
Yh_fit = median(Yh_mean,'omitnan') + sYhm *(t_dec - t0c);
WR_fit = median(WR_mean,'omitnan') + sWRm *(t_dec - t0c);

fig = figure('Color','w','Units','pixels','Position',[100 100 1200 900]);
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');

% A) Raw width trend along longitude
nexttile(1);
plot(lon, slopeW, 'Color', col_raw, 'LineWidth', 1.5); hold on;
sigW = pW < 0.05 & isfinite(slopeW);
plot(lon(sigW), slopeW(sigW),'o','MarkerSize',3,'MarkerFaceColor',col_raw,'MarkerEdgeColor','k');
yline(0,'-','Color',col_grey);
xlabel('Longitude'); ylabel('Width trend (km decade^{-1})');
title('Raw width trend along stream'); grid on; box on;

% B) Speed-explained width trend along longitude
nexttile(2);
plot(lon, slopeYh, 'Color', col_hat, 'LineWidth', 1.5); hold on;
yline(0,'-','Color',col_grey);
xlabel('Longitude'); ylabel('Explained trend (km decade^{-1})');
title('Width trend explained by speed'); grid on; box on;

% C) Residual width trend along longitude
nexttile(3);
plot(lon, slopeR, 'Color', col_res, 'LineWidth', 1.5); hold on;
sigR = pR < 0.05 & isfinite(slopeR);
plot(lon(sigR), slopeR(sigR),'s','MarkerSize',3,'MarkerFaceColor',col_res,'MarkerEdgeColor','k');
yline(0,'-','Color',col_grey);
xlabel('Longitude'); ylabel('Residual trend (km decade^{-1})');
title('Residual width trend (speed removed)'); grid on; box on;

% D) Region-mean time series with Sen fits
nexttile(4);
plot(t, W_mean,  '-', 'Color', col_raw, 'LineWidth', 1.0); hold on;
plot(t, Yh_mean, '-', 'Color', col_hat, 'LineWidth', 1.0);
plot(t, WR_mean, '-', 'Color', col_res, 'LineWidth', 1.0);
plot(t, W_fit,  '--', 'Color', col_raw, 'LineWidth', 1.5);
plot(t, Yh_fit, '--', 'Color', col_hat, 'LineWidth', 1.5);
plot(t, WR_fit, '--', 'Color', col_res, 'LineWidth', 1.5);
xlabel('Time'); ylabel('Width anomaly (km)');
title(sprintf('Region-mean width and Sen trends (raw %.2f, expl %.2f, resid %.2f km/dec)', ...
              sWm, sYhm, sWRm));
legend({'Raw','Explained by speed','Residual','Raw trend','Explained trend','Residual trend'}, ...
        'Location','best','Box','off'); grid on; box on;

% E) Scatter of width vs speed anomalies (pooled)
nexttile(5);
Wi = W(:); Ui = U(:);
ok = isfinite(Wi) & isfinite(Ui);
Wi = Wi(ok); Ui = Ui(ok);
scatter(Ui, Wi, 8, 'MarkerFaceColor',[0.75 0.75 0.75], 'MarkerEdgeColor','none'); hold on;
if exist('fitlm','file')
    mdl_sc = fitlm(Ui, Wi, 'RobustOpts','on');
    xg = linspace(min(Ui), max(Ui), 200)'; yg = predict(mdl_sc, xg);
    R2txt = sprintf('R^2 = %.2f', mdl_sc.Rsquared.Ordinary);
else
    p = polyfit(Ui, Wi, 1); xg = linspace(min(Ui), max(Ui), 200)'; yg = polyval(p, xg);
    R2txt = '';
end
plot(xg, yg, '-', 'Color', 'k','LineWidth',1.5);
xlabel('Along-jet speed anomaly (m s^{-1})'); ylabel('Width anomaly (km)');
title(sprintf('Width vs speed anomalies %s', R2txt)); grid on; box on;

% F) Fraction of width trend explained by speed
nexttile(6);
plot(lon, frac_exp, 'k-', 'LineWidth', 1.2); hold on;
yline(0,'-','Color',col_grey); yline(1,'--','Color',col_grey);
xlabel('Longitude'); ylabel('Fraction explained');
title('Fraction of width trend explained by speed'); grid on; box on;

% Export figure
try
    exportgraphics(fig, [export_base '.pdf'], 'ContentType','vector');
    exportgraphics(fig, [export_base '.png'], 'Resolution', 300);
catch
    print(fig, [export_base '.pdf'], '-dpdf', '-painters');
    print(fig, [export_base '.png'], '-dpng', '-r300');
end

%% ---- Outputs
out = struct('slopeW',slopeW,'pW',pW,'slopeR',slopeR,'pR',pR, ...
             'slopeExpl',slopeYh,'fracExplained',frac_exp,'R2',R2loc, ...
             'W_anom',W,'U_anom',U,'W_resid',W_resid,'Yhat',Yhat, ...
             't',t,'t_dec',t_dec,'lon',lon, ...
             'regionMean',struct('W',W_mean,'WR',WR_mean,'Yh',Yh_mean, ...
                                  'slopeW',sWm,'pW',pWm,'slopeWR',sWRm,'pWR',pWRm,'slopeYh',sYhm));
end

%% ===================== Local helper functions =====================
function Xa = remove_monthly_climatology(X, t)
    % Remove monthly climatology (lon x time). Returns anomalies.
    if isrow(X), X = X'; end
    [n,m] = size(X);
    Xa = nan(n,m);
    mm = month(t(:));                % column
    for k = 1:12
        idx = (mm==k).';             % row logical to match 1 x m
        Xm  = X(:,idx);
        clim = mean(Xm,2,'omitnan');
        Xa(:,idx) = Xm - clim;
    end
end

function [slope_dec, pMK] = sen_mk(y, t_dec)
    % Theil–Sen slope (y per decade) and Mann–Kendall p-value.
    % Robust to orientation mismatches (forces column vectors).
    y  = y(:);
    td = t_dec(:);
    n  = min(numel(y), numel(td));   % ensure equal length
    y  = y(1:n);
    td = td(1:n);

    good = isfinite(y) & isfinite(td);   % 1-D logical mask
    y  = y(good);
    td = td(good);
    slope_dec = nan; pMK = nan;
    n = numel(y); if n < 24, return; end

    % Pairwise slopes
    % (n=~365 -> ~66k pairs; OK)
    idx_i = []; idx_j = [];
    for i = 1:n-1
        k = (i+1):n;
        idx_i = [idx_i; i*ones(numel(k),1)]; %#ok<AGROW>
        idx_j = [idx_j; k'];                 %#ok<AGROW>
    end
    slopes = (y(idx_j)-y(idx_i)) ./ (td(idx_j)-td(idx_i));
    slope_dec = median(slopes,'omitnan');

    % Mann–Kendall with tie correction
    S = 0;
    for i = 1:n-1
        S = S + sum(sign(y((i+1):end) - y(i)));
    end
    [~,~,cnt] = unique(y);
    ties = accumarray(cnt,1);
    varS = n*(n-1)*(2*n+1)/18 - sum(ties.*(ties-1).*(2*ties+1))/18;
    if varS<=0, pMK = NaN; return; end
    if S > 0
        Z = (S-1)/sqrt(varS);
    elseif S < 0
        Z = (S+1)/sqrt(varS);
    else
        Z = 0;
    end
    % Two-sided p via erfc (no toolboxes)
    pMK = 2*0.5*erfc(abs(Z)/sqrt(2));
end
