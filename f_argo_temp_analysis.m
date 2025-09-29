%%  Pacific–Antarctic Ridge | ARGO mean temperature (climatology 2004–2018)
%  ---------------------------------------------------------------
%  Task :  (1) read RG_ArgoClim_Temperature_2004_2018.nc
%          (2) extract the 60° S–48° S, 150° W–80° W window
%          (3) visualise the 0-dbar (surface) climatological field
%  Note :  only the mean field, ARGO_TEMPERATURE_MEAN, is used.
%  ---------------------------------------------------------------

ncFile = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
    'Honours Research Project\third_meander\argo', ...
    'RG_ArgoClim_Temperature_2004_2018.nc');


%% -------------------------------------------------- 1.  Read coordinates
lon = ncread(ncFile,'LONGITUDE');          % 0 … 360
lat = ncread(ncFile,'LATITUDE');           % −90 … 90
prs = ncread(ncFile,'PRESSURE');           % 0 … 2000 dbar (58 levels)

% Convert longitude to conventional −180 … 180 for intuitive indexing
lon = mod(lon+180,360) - 180;


%% -------------------------------------------------- 2.  Spatial mask for 60 S–48 S, 150 W–80 W
ilat = find(lat >= -60 & lat <= -48);      % southern latitudes negative
ilon = find(lon >= -150 & lon <= -80);


%% -------------------------------------------------- 3.  Read mean temperature and put into P × Y × X order
Tmean = ncread(ncFile,'ARGO_TEMPERATURE_MEAN');   % returned order may be X × Y × P
Tmean(Tmean == -999) = NaN;                       % set fill to NaN

% The variable's dimensions are shown below:
% first dimension = LONGITUDE, second = LATITUDE, third = PRESSURE
dims = size(Tmean);

% Subset: pressure × selected-lat × selected-lon
Tsub = Tmean(ilon, ilat, :);

meand_temp_04_18 = Tsub;




%% Below, we combine the analyses over the entire 1993-2023 period to get the
% entire series of the temperature data in the Pacific-Antarctic Ridge
% meander region
ncFile_201901 = fullfile('C:\Users\xliu38\OneDrive - University of Tasmania', ...
    'Honours Research Project\third_meander\argo', ...
    'RG_ArgoClim_201901.nc');

Tanomaly_201901 = ncread(ncFile_201901, 'ARGO_TEMPERATURE_ANOMALY');   % returned order may be X × Y × P



%%  Pacific–Antarctic Ridge | RG Argo temperature anomalies
%  Concatenate ARGO_TEMPERATURE_ANOMALY from Jan-2019 … Jun-2023
%  File dimension order returned by ncread:
%          LONGITUDE × LATITUDE × PRESSURE  (× TIME=1)
%  We convert to  [PRESSURE × LATITUDE × LONGITUDE]
%  and then stack along a fourth (TIME) dimension.
%  ---------------------------------------------------------------

rootDir = 'C:\Users\xliu38\OneDrive - University of Tasmania\Honours Research Project\third_meander\argo';


%% 1.  Grid info from any reference file -------------------------------
refFile = fullfile(rootDir,'RG_ArgoClim_201901.nc');
lon = ncread(refFile,'LONGITUDE');     lon = mod(lon+180,360)-180;   % −180…180
lat = ncread(refFile,'LATITUDE');
prs = ncread(refFile,'PRESSURE');                                    % 58 levels

% PARidge window indices
ilat = find(lat >= -60 & lat <= -48);   nLat = numel(ilat);          % 12
ilon = find(lon >= -150 & lon <= -80);  nLon = numel(ilon);          % 70


%% 2.  Month list -------------------------------------------------------
allMon = datetime(2019,1,1) : calmonths(1) : datetime(2023,6,1);  % 54 months
nT     = numel(allMon);

% Pre-allocate   TIME × LON × LAT × PRESS
Tcat = NaN(nT, nLon, nLat, numel(prs), 'double');

skip = {};                      % record tags that could not be read
fprintf('Reading ARGO_TEMPERATURE_ANOMALY for %d months …\n', nT);

for k = 1:nT
    tag   = datestr(allMon(k),'yyyymm');
    fName = fullfile(rootDir, ['RG_ArgoClim_', tag, '.nc']);

    if ~isfile(fName)
        skip{end+1} = [tag '  (file missing)'];
        continue
    end

    try
        % Expected NetCDF order:  LONGITUDE × LATITUDE × PRESSURE × TIME(=1)
        Tanom = ncread(fName,'ARGO_TEMPERATURE_ANOMALY');
    catch ME
        skip{end+1} = [tag '  (' ME.message ')'];
        continue
    end

    Tanom(Tanom == -999) = NaN;           % replace fill value
    Tanom = squeeze(Tanom);               % → [360 145 58]

    % Check dimensions
    if size(Tanom,3) ~= numel(prs)
        skip{end+1} = [tag '  (unexpected dims)'];
        warning('Unexpected dimension order in %s (skipped)', tag);
        continue
    end

    % Subset: 70 × 12 × 58     (LON × LAT × PRESS)
    Tsub = Tanom(ilon, ilat, :);

    % Store into TIME-indexed slice
    Tcat(k, :, :, :) = double(Tsub);
end

fprintf('Finished.  anom_temp_19_23 size = %s   (skipped: %d)\n', ...
        mat2str(size(Tcat)), numel(skip));

if ~isempty(skip)
    fprintf('\nSkipped months:\n  • %s\n', skip{:});
end

% Rename final array as requested
anom_temp_19_23 = Tcat;



%% Below, we derive the 2019-2023 mean temperatures based on the 2006-2018 mean
% temperatures and then derive the 2006-2023 mean temperatures

% 1.  Ensure both variables are in double precision (optional but safe)
anom_temp_19_23  = double(anom_temp_19_23);
meand_temp_04_18 = double(meand_temp_04_18);

% 2.  Expand the 2004-18 mean so it broadcasts along the time dimension
%     (MATLAB automatically expands singleton dimensions when adding)
baseField = reshape(meand_temp_04_18, [1 70 12 58]);   % 1 × 70 × 12 × 58

% 3.  Add anomaly + mean  → absolute temperature
meand_temp_19_23 = anom_temp_19_23 + baseField;        % 54 × 70 × 12 × 58

% 4.  (Optional) verify basic statistics
%fprintf('2019-23 absolute T :  min = %.2f °C,  max = %.2f °C\n', ...
        %min(meand_temp_19_23(:),'omitnan'), ...
        %max(meand_temp_19_23(:),'omitnan'));



%% Next, we continue to load the temperature anomaly variable during the 2004-
% 2018 period and then add this temperature anomaly variable to the
% 2004-2018 mean temperature variable to get the time series of the
% 2004-2018 temperature together in the longitude, latitude and pressure

Tanom = ncread(ncFile,'ARGO_TEMPERATURE_ANOMALY');  % 360 × 145 × 58 × 180
Tanom(Tanom == -999) = NaN;

% Subset and permute → TIME × LON × LAT × PRESS   (180 × 70 × 12 × 58)
Tanom_sub = permute( Tanom(ilon, ilat, :, :), [4 1 2 3] );

baseField_04_18 = reshape(meand_temp_04_18, [1 nLon nLat numel(prs)]);  % 1 × 70 × 12 × 58
meand_temp_04_18_4D = baseField_04_18 + Tanom_sub;                        % 180 × 70 × 12 × 58



%% Finally, we continue to concate the 4-D 2004-2018 temperature variable with
% the 4-D 2019-2023 temperautre variable to get the 4-D 2004-2023
% temperature variable

% 1.  Ensure both arrays exist and share the same precision
assert(exist('meand_temp_04_18_4D','var')==1 , ...
       'meand_temp_04_18_4D not found in workspace');
assert(exist('meand_temp_19_23','var')==1     , ...
       'meand_temp_19_23 not found in workspace');

if ~strcmp(class(meand_temp_04_18_4D), class(meand_temp_19_23))
    % cast the second array to match the first
    meand_temp_19_23 = feval(class(meand_temp_04_18_4D), meand_temp_19_23);
end


% 2.  Concatenate along the TIME dimension (dimension 1)
meand_temp_04_23_4D_TIME_LON_LAT_PRES = cat(1, ...
        meand_temp_04_18_4D, ...
        meand_temp_19_23);

% 3.  Check the resulting size
outSize = size(meand_temp_04_23_4D_TIME_LON_LAT_PRES);
fprintf('Concatenated variable  meand_temp_04_23_4D_TIME_LON_LAT_PRES  size = %s\n', ...
        mat2str(outSize));



%% We then continue to save the derived temperatures-related variables into a
% matlab matrix data variable

%  1.  Create monthly datetime array  (Jan-2004 → Jun-2023)
% ===============================================================
meand_time_mon_04_23 = datetime(2004,1,1) : calmonths(1) : datetime(2023,6,1);

% sanity-check – should match TIME dimension of concatenated cube
assert(numel(meand_time_mon_04_23) == size(meand_temp_04_23_4D_TIME_LON_LAT_PRES,1), ...
      'Time axis length (%d) does not match data cube (%d).', ...
      numel(meand_time_mon_04_23), size(meand_temp_04_23_4D_TIME_LON_LAT_PRES,1));

% ===============================================================
%  2.  Save selected variables to MAT-file
% ===============================================================
outFile = 'meand_temp_mon_04_23_4D.mat';
save(outFile, 'ilat', 'ilon', 'meand_time_mon_04_23', ...
              'meand_temp_04_23_4D_TIME_LON_LAT_PRES', '-v7.3');

fprintf('Saved variables to  %s\n', fullfile(pwd,outFile));

