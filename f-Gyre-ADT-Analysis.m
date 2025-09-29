%% ======================================================================
%  CMEMS monthly SLA   (42°–38° S , 150° E–70° W)   1993-03 → 2023-06
%  Variable order in file:  LON × LAT × TIME ; units = metres
% ======================================================================

cmems_gyre_adt_19930101_20031231 = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
                  'Honours Research Project\third_meander', ...
                  'cmems_gyre_adt_19930101_20031231.nc');

cmems_gyre_adt_20040101_20131231 = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
                  'Honours Research Project\third_meander', ...
                  'cmems_gyre_adt_20040101_20131231.nc');

cmems_gyre_adt_20140101_20230630 = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
                  'Honours Research Project\third_meander', ...
                  'cmems_gyre_adt_20140101_20230630.nc');


%% 1.  Read time, lon, lat, ADT ------------------------------------------

lon     = ncread(cmems_gyre_adt_19930101_20031231,'longitude'); % 1120 x 1
lat     = ncread(cmems_gyre_adt_19930101_20031231,'latitude'); % 32 x 1

time_s_19930101_20031231  = double(ncread(cmems_gyre_adt_19930101_20031231,'time'));           % seconds since 1970-01-01
time_19930101_20031231    = datetime(1970,1,1) + seconds(time_s_19930101_20031231);    % 4017 daily steps
adt_raw_19930101_20031231 = ncread(cmems_gyre_adt_19930101_20031231,'adt');                    % 1120 × 32 × 4017  (m)

time_s_20040101_20131231  = double(ncread(cmems_gyre_adt_20040101_20131231,'time'));           % seconds since 1970-01-01
time_20040101_20131231    = datetime(1970,1,1) + seconds(time_s_20040101_20131231);    % 3653 daily steps
adt_raw_20040101_20131231 = ncread(cmems_gyre_adt_20040101_20131231,'adt');                    % 1120 × 32 × 3653  (m)

time_s_20140101_20230630  = double(ncread(cmems_gyre_adt_20140101_20230630,'time'));           % seconds since 1970-01-01
time_20140101_20230630    = datetime(1970,1,1) + seconds(time_s_20140101_20230630);    % 3468 daily steps
adt_raw_20140101_20230630 = ncread(cmems_gyre_adt_20140101_20230630,'adt');                    % 1120 × 32 × 3468  (m)


%% 2.  Area-mean monthly SLA (m)
% average first over longitude (dim 1), then latitude (dim 2)
adt_raw_19930101_20031231_mean = squeeze( mean( mean(adt_raw_19930101_20031231,1,'omitnan'), 2,'omitnan') ); % 4017×1
adt_raw_20040101_20131231_mean = squeeze( mean( mean(adt_raw_20040101_20131231,1,'omitnan'), 2,'omitnan') ); % 3653×1
adt_raw_20140101_20230630_mean = squeeze( mean( mean(adt_raw_20140101_20230630,1,'omitnan'), 2,'omitnan') ); % 3468×1

% Concate all the average-mean ADT data time series
% Ensure column shape and concatenate
adt_raw_19930101_20230630_mean = [
    adt_raw_19930101_20031231_mean(:);
    adt_raw_20040101_20131231_mean(:);
    adt_raw_20140101_20230630_mean(:)
];

time_19930101_20230630 = [
    time_19930101_20031231;...
    time_20040101_20131231;...
    time_20140101_20230630
    ];


%% 3.  Remove seasonal cycle (12-month climatology)
mo           = month(time_19930101_20230630);                         % 1…12
clim12       = accumarray(mo, adt_raw_19930101_20230630_mean, [12 1], @mean, NaN);
adt_raw_19930101_20230630_anom     = adt_raw_19930101_20230630_mean - clim12(mo);               % deseasoned anomaly (m)


%% 4.  Trend diagnostics (Sen, MK, modified MK, OLS, R²) ------------------
timeNum  = datenum(time_19930101_20230630);
N        = numel(timeNum);
monIdx   = (0:N-1)';                                % months since start
xC       = monIdx - mean(monIdx);                   % centred months

% Sen slope
sen_mo   = Sen_Slope(adt_raw_19930101_20230630_anom);                     % m per month
sen_dec  = sen_mo * 120;                            % m per decade

% Mann–Kendall tests
[~, pMK]  = Mann_Kendall(adt_raw_19930101_20230630_anom);
[~, pMKm] = Modified_Mann_Kendall(adt_raw_19930101_20230630_anom);

% OLS slope & R²
beta_mo  = (xC' * (adt_raw_19930101_20230630_anom - mean(adt_raw_19930101_20230630_anom))) / (xC' * xC);
beta_dec = beta_mo * 120;
fitLine  = beta_mo * xC + mean(adt_raw_19930101_20230630_anom);
R2       = 1 - sum((adt_raw_19930101_20230630_anom - fitLine).^2) / sum((adt_raw_19930101_20230630_anom - mean(adt_raw_19930101_20230630_anom)).^2);

%% 5.  Console summary
fprintf('\n===  Area-mean ADT anomaly (metres)  42°–38°S  ===\n');
fprintf('Sen slope   : %+8.4f m dec⁻¹\n',  sen_dec);
fprintf('OLS slope   : %+8.4f m dec⁻¹\n',  beta_dec);
fprintf('R² (OLS)    : %.3f\n',            R2);
fprintf('M-K  p      : %.3f\n',            pMK);
fprintf('Mod M-K p   : %.3f\n',            pMKm);

%% 6.  Plot using existing helper (units = metres)
plotWithTrend(timeNum, monIdx, {adt_raw_19930101_20230630_anom}, ...
              {[0 0 0]}, {'Area-mean ADT'}, ...
              'Absolute Dynamic Topography (m)', ' m dec^{-1}');
