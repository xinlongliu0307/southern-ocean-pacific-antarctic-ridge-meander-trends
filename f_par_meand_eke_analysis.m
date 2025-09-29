%% Quick check: dimensions of ugos / vgos in one CMEMS file
% -----------------------------------------------------------
% Adjust the path below if your folder structure is different.
exFile = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
                  'Honours Research Project\third_meander', ...
                  'cmems_ugos_vgos_19930101_19971231.nc');

% read the two variables
ug = ncread(exFile,'ugos');    % TIME × LAT × LON
vg = ncread(exFile,'vgos');

% report their sizes
fprintf('ugos size = %s  (LON × LAT × TIME)\n', mat2str(size(ug)));
fprintf('vgos size = %s  (LON × LAT × TIME)\n', mat2str(size(vg)));


%% =====================================================================
%  CONCATENATE CMEMS UGOS / VGOS  (1993-01-01 … 2023-06-01)
%  ---------------------------------------------------------------------
%  • Each file: ugos(LON,LAT,TIME), vgos(LON,LAT,TIME)
%  • Convert to TIME × LAT × LON, then cat(1,…) along TIME
% =====================================================================
rootDir = 'C:\Users\xliu38\OneDrive - University of Tasmania\Honours Research Project\third_meander';
wildPat = fullfile(rootDir,'cmems_ugos_vgos_*.nc');
files   = dir(wildPat);

% ---------- keep & sort files with *_YYYYMMDD_YYYYMMDD.nc pattern
startDates = arrayfun(@(f) parseStartDate(f.name), files, 'UniformOutput', true);

isValid = ~isnat(startDates);      % logical index of good names
files   = files(isValid);
startDates = startDates(isValid);  % keep matching dates only

[~, ord] = sort(startDates);       % chronological order
files   = files(ord);

% read coordinate vectors once
lon = ncread(fullfile(rootDir,files(1).name),'longitude');   % 560
lat = ncread(fullfile(rootDir,files(1).name),'latitude');    %  96

ugosCell = {};  vgosCell = {};  timeCell = {};

fprintf('Reading %d files …\n',numel(files));
for k = 1:numel(files)
    fn = fullfile(rootDir,files(k).name);
    fprintf('  %s\n',files(k).name);

    % ---- TIME (seconds since 1970-01-01) -----------------------------
    tsec = double(ncread(fn,'time'));
    timeCell{k} = datetime(1970,1,1) + seconds(tsec);  % 1×Ntime

    % ---- read, permute to TIME × LAT × LON --------------------------
    ug = permute( double(ncread(fn,'ugos')), [3 2 1]); % T × 96 × 560
    vg = permute( double(ncread(fn,'vgos')), [3 2 1]); % same order
    ugosCell{k} = ug;
    vgosCell{k} = vg;
end

ugos_all = cat(1, ugosCell{:});         % [Nt_total × 96 × 560]
vgos_all = cat(1, vgosCell{:});
time_all = vertcat(timeCell{:}).';      % column datetime

fprintf('ugos_all size = %s (TIME×LAT×LON)\n', mat2str(size(ugos_all)));
fprintf('vgos_all size = %s (TIME×LAT×LON)\n', mat2str(size(vgos_all)));

% save to MAT-file
outMat = fullfile(rootDir,'cmems_ugos_vgos_concat_1993_2023.mat');
save(outMat,'ugos_all','vgos_all','time_all','lat','lon','-v7.3');
fprintf('Merged arrays saved to  %s\n', outMat);



%% Below, we calculate the 4-month running mean velocities of the zonal and
% meridional current velocities, respectively
load('C:\Users\xliu38\OneDrive - University of Tasmania\Honours Research Project\third_meander\cmems_ugos_vgos_concat_1993_2023.mat');

x_months = 4;

time_all_num = datenum(time_all);

% Define the monthly timeseries starting from the 15th of the first month
% to 15th of the last months of data:
t_temp = datetime(time_all_num(1)+14,'ConvertFrom','datenum')...
    :calmonths(1)...
    :datetime(time_all_num(end),'ConvertFrom','datenum');

time_all_num_monthly = datenum(t_temp);

% Find the index of each monthly step in daily time stamp:
[sharedvals, idx_time_daily] = intersect(time_all_num,...
    time_all_num_monthly, 'stable');

% Convert months into days (nearest even number)
x_days = 2*floor(x_months*30.417/2);

% Below, we adjust the strucutre of the daily zonal and meridional concated
% geostrophic current velocity variables
ugos_all = permute(ugos_all,[3 2 1]);   % 560 × 96 × 11109
vgos_all = permute(vgos_all,[3 2 1]);   % same order


% For each month make a running sum of the flag values over x_days 
ugos_all_mon = NaN(length(lon),...
    length(lat),...
    length(time_all_num_monthly));

vgos_all_mon = NaN(length(lon),...
    length(lat),...
    length(time_all_num_monthly));

for i = floor(x_months/2+1) : length(time_all_num_monthly)-floor(x_months/2)

    ugos_all_mon(:,:,i) = sum(ugos_all(:,:,idx_time_daily(i)-(x_days/2-1):idx_time_daily(i)+(x_days/2)),3);
    vgos_all_mon(:,:,i) = sum(vgos_all(:,:,idx_time_daily(i)-(x_days/2-1):idx_time_daily(i)+(x_days/2)),3);

end


%% Finally, we save the calculated monthly zonal and meridional surface
% zonal and meridional current velocities variables into one dataset

outFile = fullfile(rootDir,'cmems_ugos_vgos_mon_1993_2023.mat');
save(outFile,'ugos_all_mon','vgos_all_mon','time_all_num_monthly','lat','lon','-v7.3');
fprintf('Merged arrays saved to  %s\n', outFile);


%% Below, we calculate and plot the spatial and temporal analysis of the
% eddy kinetic energy in the Pacific-Antarctic Ridge region

%  MONTHLY EKE AND 1993–2023 MEAN  (Pacific–Antarctic Ridge window)
%  --------------------------------------------------------------
%  Prerequisites already produced by earlier code:
%     ugos_all_mon   % LON × LAT × Nmon   (m s-1)
%     vgos_all_mon   % LON × LAT × Nmon   (m s-1)
%     lon            % 560 × 1            (deg E, centred −180…180)
%     lat            %  96 × 1            (deg N, −80…80 ~)
% ===============================================================

% 1.  Convert running–sum fields to running–mean (divide by #days)
% x_days was chosen above:
ugos_all_mon_avg = ugos_all_mon / x_days;     % LON × LAT × Nmon
vgos_all_mon_avg = vgos_all_mon / x_days;


% 2.  Eddy–component anomalies u' , v'  (remove temporal mean at each grid)
u_clim = mean(ugos_all_mon_avg, 3, 'omitnan');    % LON × LAT
v_clim = mean(vgos_all_mon_avg, 3, 'omitnan');

u_prime = ugos_all_mon_avg - u_clim;              % anomaly cube
v_prime = vgos_all_mon_avg - v_clim;


% 3.  Monthly Eddy Kinetic Energy
EKE_mon = 0.5 * ( u_prime.^2 + v_prime.^2 );  % m² s⁻²


% 4.  1993–2023 time-mean EKE
EKE_mean = mean(EKE_mon, 3, 'omitnan');       % LON × LAT

% Transverse the latitude variable for later plotting
lat = lat.';

% Build matching grids for PCOLOR
[LonG, LatG] = meshgrid(lon, lat);   % 96 × 560

% Transpose EKE_mean so its dimensions are LAT × LON → 96 × 560
Z = EKE_mean.';                      % 96 × 560

% Plot of the 1993-2023 temporal mean spatial distribution of the monthly
% eddy kinetic energy
figure('Color','w');
pcolor(LonG, LatG, Z);
shading flat
colormap(cmocean('dense'));                 % energetic palette
set(gca,'YDir','normal','FontSize',9);
xlabel('Longitude (^{\circ})');
ylabel('Latitude (^{\circ})');
%title('Mean Eddy Kinetic Energy 1993–2023  (m^{2} s^{-2})','FontWeight','bold');

hold on
% 1. vertical longitude lines (western longitudes are negative)
xline(-135.1 ,'k--','LineWidth',1.5);   % Upper / Lower boundary
xline(-109.4,'k--','LineWidth',1.5);    % Lower / Downstream boundary

% 2. Flat-Region rectangle (96.4–107.6 °W , 59.6–57.4 °S)
lonFlat = [-107.6  -96.4];
latFlat = [ -59.6   -57.4];
plot([lonFlat(1) lonFlat(2) lonFlat(2) lonFlat(1) lonFlat(1)], ...
     [latFlat(1) latFlat(1) latFlat(2) latFlat(2) latFlat(1)], ...
     'k--','LineWidth',1.5);

hold off

cb = colorbar('eastoutside');
cb.Label.String   = 'Eddy Kinetic Energy (m^{2} s^{-2})';
cb.Label.FontSize = 10;



%% Below, we continue to calculate and plot the spatial distributions of
% decadal-mean trends of eddy kinetic energy in the Pacific-Antarctic Ridge
% region

% 1.  Prepare time vector for regression
tDec   = year(t_temp) + (month(t_temp)-1)/12;   % decimal years
tD     = tDec - mean(tDec);                     % demean
Sxx    = sum(tD.^2);
Nmon   = numel(tDec);


% 2.  Trend & p-value at every grid cell
[LONn,LATn] = size(EKE_mon(:,:,1));     % 560 × 96

slope = NaN(LONn,LATn);                % m² s⁻² dec⁻¹
pVal  = NaN(LONn,LATn);

for i = 1:LONn
    for j = 1:LATn
        y = squeeze(EKE_mon(i,j,:));   % vector length = Nmon
        y = y(:);                      % make column

        good = isfinite(y);
        if sum(good) < 12,  continue,  end  % at least one year of data

        yGood = y(good);               yGood = yGood(:);         % column
        tGood = tD(good);              tGood = tGood(:);         % column
        tGood = tGood - mean(tGood);                          % demean

        % ----- OLS slope (per decade) -------------------------------
        beta  = (tGood' * yGood) / (tGood' * tGood);  % m² s⁻² yr⁻¹
        slope(i,j) = beta * 10;                       % → per decade

        % ----- two-sided t-test p-value -----------------------------
        yhat  = beta * tGood;
        RSS   = sum((yGood - yhat).^2);
        SEb   = sqrt( RSS / (numel(yGood)-2) / (tGood' * tGood) );
        tStat = beta / SEb;
        pVal(i,j) = 2*(1 - tcdf(abs(tStat), numel(yGood)-2));
    end
end


% 3.  Plot: decadal trend + significance (p < 0.05)
[LonG,LatG] = meshgrid(lon,lat);      % 96 × 560  (LAT × LON)

Z = slope.';                          % LAT × LON
P = pVal.';                           % same

% duplicate edges for pcolor
LonP = [LonG ,  LonG(:,end)]; LonP = [LonP ; LonP(end,:)];
LatP = [LatG ,  LatG(:,end)]; LatP = [LatP ; LatP(end,:)];
ZP   = [Z    ,  Z(:,end)  ];  ZP   = [ZP   ; ZP(end,:) ];
PP   = [P    ,  P(:,end)  ];  PP   = [PP   ; PP(end,:) ];

figure('Color','w');
pcolor(LonP,LatP,ZP);
shading flat
colormap(cmocean('balance'));         % diverging palette
caxis([-prctile(abs(Z(:)),95) prctile(abs(Z(:)),95)]);   % symmetric limits
set(gca,'YDir','normal','FontSize',9);
xlabel('Longitude (^{\circ})');  ylabel('Latitude (^{\circ})');
%title('Decadal trend of EKE 1993–2023  (m^{2} s^{-2} dec^{-1})','FontWeight','bold');

% overlay dots where p < 0.05
hold on
%sig = PP < 0.05;
%scatter(LonP(sig), LatP(sig), 8, 'k','filled');
% 1. vertical longitude lines (western longitudes are negative)
xline(-135.1 ,'k--','LineWidth',1.5);   % Upper / Lower boundary
xline(-109.4,'k--','LineWidth',1.5);    % Lower / Downstream boundary

% 2. Flat-Region rectangle (96.4–107.6 °W , 59.6–57.4 °S)
%lonFlat = [-107.6  -96.4];
%latFlat = [ -59.6   -57.4];
plot([lonFlat(1) lonFlat(2) lonFlat(2) lonFlat(1) lonFlat(1)], ...
     [latFlat(1) latFlat(1) latFlat(2) latFlat(2) latFlat(1)], ...
     'k--','LineWidth',1.5);
hold off

cb = colorbar('eastoutside');
cb.Label.String   = 'Eddy Kinetic Energy Decadal Trend (m^{2} s^{-2} dec^{-1})';
cb.Label.FontSize = 10;



%% Finally, we divide the eddy kinetic energy into the meander's different
% sections as well as plot their time series and estimate their linear
% trends and test their statistical significance

% 1.  Time-axis helpers
Nt        = numel(t_temp);           % 361 months  (1993-01 … 2023-06)
timeNum   = datenum(t_temp);         % datenum for plotting
monthIdx  = (0:Nt-1)';               % 0,1,2,… used by Sen_Slope


% 2.  Longitude & latitude masks
maskWholeLon = true(size(lon));
maskUpperLon =  lon >= -150   & lon <= -135.1;          % 135.1–150°W
maskLowerLon =  lon >  -135.1 & lon <= -109.4;          % 109.4–135.1°W
maskDownLon  =  lon >  -109.4 & lon <=  -80;            %  80.0–109.4°W
maskFlatLon  =  lon >= -107.6 & lon <=  -96.4;          %  96.4–107.6°W

maskAllLat   = true(size(lat));
maskFlatLat  = lat >= -59.6   & lat <=  -57.4;          % 59.6–57.4°S

lonMasks = {maskWholeLon, maskUpperLon, maskLowerLon, maskDownLon, maskFlatLon};
latMasks = {maskAllLat   , maskAllLat , maskAllLat  , maskAllLat , maskFlatLat };

sectLbl  = {'Whole Meander','Upper-Ridge Section','Lower-Ridge Section', ...
            'Downstream Section','Flat Region'};
nSec     = numel(lonMasks);


% 3.  Monthly EKE series for each section (mean over lon & lat)
series = NaN(Nt,nSec);

for s = 1:nSec
    Ek = EKE_mon(lonMasks{s}, latMasks{s}, :);         % L×M×T
    % average over lon (dim-1) and lat (dim-2)
    series(:,s) = squeeze( mean( mean(Ek, 2,'omitnan'), 1,'omitnan') );
end


% 4.  Trend table (Sen, OLS, MK, modified MK)
tbl = cell(nSec+1,6);
tbl(1,:) = {'Section','Sen slope (m² s⁻² dec⁻¹)','OLS slope','p_MK','p_MMK','R²'};

xDec = year(t_temp) + (month(t_temp)-1)/12;
xD   = xDec - mean(xDec);          % demeaned decimal-year time

for s = 1:nSec
    y = series(:,s);
    good = isfinite(y);

    % --- Sen slope -----------------------------------------------------
    sen_raw = Sen_Slope(y(good));       % per month
    sen_dec = sen_raw * 120;            % per decade

    % --- OLS slope -----------------------------------------------------
    xCol = xD(good);  xCol = xCol(:);   % force column
    yCol = y(good);   yCol = yCol(:);

    beta_yr  = (xCol' * (yCol - mean(yCol))) / (xCol' * xCol);  % m² s⁻² yr⁻¹
    beta_dec = beta_yr * 10;

    yhat = beta_yr * xCol + mean(yCol);
    RSS  = sum((yCol - yhat).^2);
    R2   = 1 - RSS / sum((yCol - mean(yCol)).^2);

    % --- MK & modified MK ---------------------------------------------
    [~,pMK]  = Mann_Kendall(yCol);
    [~,pMMK] = Modified_Mann_Kendall(yCol);

    % Fill row in table
    tbl{s+1,1} = sectLbl{s};
    tbl{s+1,2} = sprintf('%+.3e', sen_dec);
    tbl{s+1,3} = sprintf('%+.3e', beta_dec);
    tbl{s+1,4} = sprintf('%.3f',  pMK);
    tbl{s+1,5} = sprintf('%.3f',  pMMK);
    tbl{s+1,6} = sprintf('%.2f',  R2);
end

disp('EKE trend statistics (1993–2023)');
disp(tbl);


% 5.  Plot five series with Sen-line & annotations
colSet = {[0 0 0], [1 0 1], [1 0.5 0], [0 0 1], [0 1 1]};
plotWithTrend(timeNum, monthIdx, ...
              arrayfun(@(k) series(:,k), 1:nSec, 'UniformOutput',false), ...
              colSet, sectLbl, ...
              'Eddy Kinetic Energy (m^{2} s^{-2})', ' m^{2} s^{-2} dec^{-1}');





%% -------- local helper -------------------------------------------------
function dt = parseStartDate(fname)
    tok = regexp(fname,'_(\d{8})_\d{8}\.nc$','tokens','once');
    if isempty(tok)
        dt = NaT;
    else
        dt = datetime(tok{1},'InputFormat','yyyyMMdd');
    end
end

