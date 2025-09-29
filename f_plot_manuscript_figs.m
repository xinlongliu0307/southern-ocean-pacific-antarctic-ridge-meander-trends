%% Figure 1: figure of the ADT gradients in the entire Southern Ocean and East Pacific Ridge Meander region

ncfile = 'C:\Users\xliu38\OneDrive - University of Tasmania\Honours Research Project\third_meander\cmems_30s_90s_180w_180e_20080320_adt_ugos_vgos.nc';

% To see the content of this .nc file:
alti_ADT_SO_20080320 = ncread(ncfile,'adt',[1 1 1],[Inf Inf Inf]); 
alti_lon_SO_20080320 = ncread(ncfile,'longitude');
alti_lat_SO_20080320 = ncread(ncfile,'latitude')';
alti_time_SO_20080320 = double(ncread(ncfile,'time')') + datenum('1993-01-01 00:00:00');


% We calculate the ADT gradients in the entire Southern Ocean on March
% 20th, 2008
[dx_SO_20080320, dy_SO_20080320] = gradient(alti_ADT_SO_20080320(:,:), 0.25); % 0.25 is the grid size in lat and lon
% Adjust the gradient to be in meters:
dx_ADT_SO_20080320 = dx_SO_20080320./(pi/180*Rearth*cos(deg2rad(alti_lat_SO_20080320)));
dy_ADT_SO_20080320 = dy_SO_20080320./(pi/180*Rearth);
% Overall gradient
d_ADT_SO_20080320(:,:) = sqrt(dx_ADT_SO_20080320.^2+dy_ADT_SO_20080320.^2); % in m/m


% Below we plot the spatial distributions of ADT in both the entire
% Southern Ocean and the East Pacific Ridge Meander region
i = 5558; % 5558 is the index number of March 20th, 2008
date_plot = datestr(alti_datetime(i));

figure('Name','1','Visible',fig);
hold on

% ADT in the entire Southern Ocean region
subplot(3,1,1)
pcolor(repmat(alti_lon_SO_20080320,1,length(alti_lat_SO_20080320)),repmat(alti_lat_SO_20080320,length(alti_lon_SO_20080320),1),d_ADT_SO_20080320(:,:)*100000)% from m/m to m/100km
% caxis([0 1.5]);
shading flat;
hcb2=colorbar;
xlim([-180 180]);
ylim([-90 -30]);
colormap(cmocean('dense'));
set(get(hcb2,'xlabel'),'string','ADT gradients (m/100km)','fontsize',10,'FontWeight','normal');
hold on
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[0.7 0.7],'r','Linewidth',1,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_min ADT_min],'k','Linewidth',0.5,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_max ADT_max],'k','Linewidth',0.5,'linestyle','-')
mapshow(coastlon, coastlat, 'DisplayType','polygon','FaceColor','black')
%xlabel('Longitude');
ylabel('Latitude');
% title('Gradient of ADT in m/100km')
% Add date
%annotation('textbox','String',['' date_plot ''],'FontSize',7,'FitBoxToText','off','LineStyle','none','FontWeight','normal','Position',[0.80 0.95 0.15 0.05],'backgroundcolor','none');
rectangle('Position',[-150 -70 90 40], ...   % [xmin ymin width height]
          'EdgeColor',[1 0.5 0], ...         % orange
          'LineWidth',2,'LineStyle','-');


% ADT gradients in the general East Pacific Ridge Meander region
subplot(3,1,2)
pcolor(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),d_ADT(:,:,i)*100000)% from m/m to m/100km
% caxis([0 1.5]);
shading flat;
hcb2=colorbar;
xlim([-150 -60]);
ylim([-70 -30]);
colormap(cmocean('dense'));
set(get(hcb2,'xlabel'),'string','ADT gradient (m/100km)','fontsize',10,'FontWeight','normal');
hold on
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[0.7 0.7],'r','Linewidth',1,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_min ADT_min],'k','Linewidth',0.5,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_max ADT_max],'k','Linewidth',0.5,'linestyle','-')
mapshow(coastlon, coastlat, 'DisplayType','polygon','FaceColor','black')
%xlabel('Longitude');
%ylabel('Latitude');
% title('Gradient of ADT in m/100km')
% Add date
%annotation('textbox','String',['' date_plot ''],'FontSize',7,'FitBoxToText','off','LineStyle','none','FontWeight','normal','Position',[0.80 0.95 0.15 0.05],'backgroundcolor','none');
rectangle('Position',[-150 -60 70 12], ...   % lon:‑150→‑80, lat:‑60→‑48
          'EdgeColor',[1 1 0], ...
          'LineWidth',2,'LineStyle','-');


% ADT gradients in the most specific East Pacific Ridge Meander area
subplot(3,1,3)
pcolor(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),d_ADT(:,:,i)*100000)% from m/m to m/100km
% caxis([0 1.5]);
shading flat;
hcb2=colorbar;
xlim([-150 -80]);
ylim([-60 -48]);
colormap(cmocean('dense'));
set(get(hcb2,'xlabel'),'string','ADT gradient (m/100km)','fontsize',10,'FontWeight','normal');
hold on
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[0.7 0.7],'r','Linewidth',1,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_min ADT_min],'k','Linewidth',0.5,'linestyle','-')
%contour(repmat(alti_lon,1,length(alti_lat)),repmat(alti_lat,length(alti_lon),1),alti_adt_day_1993_2023(:,:,i),[ADT_max ADT_max],'k','Linewidth',0.5,'linestyle','-')
mapshow(coastlon, coastlat, 'DisplayType','polygon','FaceColor','black')
xlabel('Longitude');
%ylabel('Latitude');
% title('Gradient of ADT in m/100km')
% Add date
annotation('textbox','String',['' date_plot ''],'FontSize',9,'FitBoxToText','off','LineStyle','none','FontWeight','normal','Position',[0.80 0.95 0.15 0.05],'backgroundcolor','none');





%% Figure 2: meander monthly and 1993-2023 mean positions; monthly position and width; transect of meander frequency

i=183; % i indicates the time: monthly from Mar 1993 to Mar 2023
ii=62; % ii indicates the specific longitude line: 150oW to 80oW

date_plot=datestr(meand_time(i),'mmm yyyy');


figure('Name','4','Visible',fig,'Position',[30 100 600 600]);

% subplot: meander monthly position and corresponding width plotting
subplot(3,1,2);
hold on
pcolor(repmat(meand_lon,1,length(meand_lat)),repmat(meand_lat,length(meand_lon),1),meand_sumfront(:,:,i))
plot([meand_lon(ii) meand_lon(ii)],[meand_S_lat(ii,i) meand_N_lat(ii,i)],'yellow','Linewidth',2.5,'linestyle','-')
%plot([meand_lon(ii) meand_lon(ii)],[-50 -30],'cyan','Linewidth',1.5,'linestyle','--')
plot(meand_loc_lon(ii,i),meand_loc_lat(ii,i),'pentagram','MarkerFaceColor',[1 0 0],'markersize',11,'MarkerEdgeColor',[0 0 0])
shading flat; hcb1=colorbar; caxis([0 x_days]);
mapshow(coastlon, coastlat, 'DisplayType','polygon','FaceColor','black')
set(get(hcb1,'xlabel'),'string','Meander frequency','fontsize',12,'FontWeight','normal');
xlim([meand_lon(ii)-8 meand_lon(ii)+8]); ylim([meand_loc_lat(ii,i)-4 meand_loc_lat(ii,i)+4])
xlabel('Longitude');
%ylabel('Latitude');
colormap(cmocean('dense'));
%title(['Meander position and width (' int2str(x_months) ' month mean)'])
annotation('textbox','String',['' date_plot ''],'FontSize',9,'FitBoxToText','off','LineStyle','none','FontWeight','normal','Position',[0.80 0.95 0.15 0.05],'backgroundcolor','none');
%annotation('textbox','String',['Relative detection threshold ' int2str(relat_thresh*100) '%'],'FontSize',8,'FitBoxToText','off','LineStyle','none','FontWeight','normal','Position',[0.65 0.50 0.35 0.05],'backgroundcolor','none');


% subplot: transect of meander frequency
subplot(3,1,3);
hold on
plot(meand_sumfront_masked(ii,:,i),meand_lat,'color','Black','linewidth',1)
%plot(pks,locs,'x','markersize',10,'MarkerEdgeColor',[1 0 0])
plot([meand_maxf(ii,i) meand_maxf(ii,i)],[meand_S_lat(ii,i) meand_N_lat(ii,i)],'r','Linewidth',2,'linestyle','-')
plot(meand_maxf(ii,i),meand_loc_lat(ii,i),'pentagram','MarkerFaceColor',[1 0 0],'markersize',11,'MarkerEdgeColor',[0 0 0])
xlabel('Meander frequency');
ylabel('Latitude');
xlim([0 x_days]);
%title(['Transect of meander frequency at ' num2str(meand_lon(ii)) '^{o}W ']);





%% Figure 3: three decadal-mean positions of the East-Pacific Ridge Meander during 1993-2023 period
%i=235;
%date_plot=datestr(meand_time(i),'mmm yyyy');

figure('Name','3','Visible',fig,'Position',[30 100 800 400]);
hold on

% Plot the colormap of the bathymetry of the East-Pacific Ridge meander
% study region
[X, Y] = meshgrid(mx, my);
pcolor(X, Y, mz2);
shading flat; % Optional: removes grid lines for cleaner appearance
%colormap(cmocean('deep'));
colormap('pink');
caxis([-5000 0]); % Correct usage: set color limits only
hcb1 = colorbar;
% Reverse the direction of the colorbar (labels from 0 to -6000)
set(hcb1, 'YDir', 'reverse');
set(hcb1, 'Ticks', -5000:500:0); % Optional: control tick intervals
ylabel(hcb1, 'Bathymetry (m)');   % Optional: label the colorbar


% Plot the three decadal-mean positions of the East-Pacific Ridge meander
% Monthly positions in every 18-month interval
plot(meand_loc_lon(:,1:18:end), meand_loc_lat(:,1:18:end),'Color',[0 1 1],'linewidth',0.001) % plot every 24 months
% 1st Decade: Jan 1993 – Dec 2002
plot(nanmean(meand_loc_lon(:,1:120),2), nanmean(meand_loc_lat(:,1:120),2),'Color', [1 1 1],'linewidth',1.5);
% 2nd Decade: Jan 2003 – Dec 2012
plot(nanmean(meand_loc_lon(:,121:240),2), nanmean(meand_loc_lat(:,121:240),2),'Color', [1 1 1],'linewidth',2.0);
% 3rd Period: Jan 2013 – May 2023
plot(nanmean(meand_loc_lon(:,241:end),2), nanmean(meand_loc_lat(:,241:end),2),'Color', [1 1 1],'linewidth',2.5);


% Coastline
mapshow(coastlon, coastlat, 'DisplayType','polygon','FaceColor','black')
xlabel('Longitude (^{\circ})', FontSize=12);
ylabel('Latitude (^{\circ})', FontSize=12);
xlim([-150 -80]);
ylim([-62 -48]);





%% Below, we first divide the meander's different sections as well as peaks and troughs
% into the corresponding meander's longitude dimensions separately.
meand_loc_lat_all_mean = nanmean(meand_loc_lat, 1);
meand_loc_lat_all_mean  = squeeze(meand_loc_lat_all_mean).';   % 1×365

meand_loc_lat_upp = meand_loc_lat(1:60, :);
meand_loc_lat_low = meand_loc_lat(61:163, :);
meand_loc_lat_down = meand_loc_lat(164:end, :);
meand_loc_lat_flat = meand_loc_lat(170:215, :);

meand_posi_pk1 = meand_loc_lat(42:48, :);
meand_posi_pk2 = meand_loc_lat(83:89, :);
meand_posi_pk3 = meand_loc_lat(97:101, :);

meand_posi_tr1 = meand_loc_lat(58:63, :);
meand_posi_tr2 = meand_loc_lat(110:114, :);
meand_posi_tr3 = meand_loc_lat(136:140, :);


% Below, we derive the one-dimensional time series of the four sections of
% the Pacific-Antarctic Ridge meander
meand_loc_lat_upp_mean = nanmean(meand_loc_lat_upp, 1);
meand_loc_lat_upp_mean  = squeeze(meand_loc_lat_upp_mean).';

meand_loc_lat_low_mean = nanmean(meand_loc_lat_low, 1);
meand_loc_lat_low_mean  = squeeze(meand_loc_lat_low_mean).';

meand_loc_lat_down_mean = nanmean(meand_loc_lat_down, 1);
meand_loc_lat_down_mean  = squeeze(meand_loc_lat_down_mean).';

meand_loc_lat_flat_mean = nanmean(meand_loc_lat_flat, 1);
meand_loc_lat_flat_mean  = squeeze(meand_loc_lat_flat_mean).';


% Initial plot of the time series of the latitude positions of the
% meander's different sections
figure('Color','w');
hold on;

plot(meand_time, meand_loc_lat_all_mean, '-k', 'LineWidth',0.5);
plot(meand_time, meand_loc_lat_upp_mean, '-m', 'LineWidth',0.5);
plot(meand_time, meand_loc_lat_low_mean, '-', 'Color', '[1 0.5 0]','LineWidth',0.5);
plot(meand_time, meand_loc_lat_down_mean, '-b', 'LineWidth',0.5);
plot(meand_time, meand_loc_lat_flat_mean, '-c', 'LineWidth',0.5);

xlabel('Year');
ylabel('Latitude (^{\circ}S)');
%title('Meander Location Trend (1993 – 2023)');
legend({'Whole Meander','Upper-Ridge Section','Lower-Ridge Section','Downstream Section', 'Flat Region'},...
    'Location','best');
grid on;
xlim([meand_time(1) meand_time(end)]);  % (optional) full x-axis extent

set(gca,'FontSize',9,'XTick',meand_time(1):365:meand_time(end));
datetick('x','yyyy','keepticks','keeplimits');



%% ----------------------------------------------------------------------
%  INPUT vectors (1×365) already prepared:
%    meand_time  (datetime or datenum)
%    meand_loc_lat_all_mean,  meand_loc_lat_upp_mean,
%    meand_loc_lat_low_mean,  meand_loc_lat_down_mean,
%    meand_loc_lat_flat_mean
% ----------------------------------------------------------------------

seriesList  = { meand_loc_lat_all_mean(:), ...
                meand_loc_lat_upp_mean(:), ...
                meand_loc_lat_low_mean(:), ...
                meand_loc_lat_down_mean(:), ...
                meand_loc_lat_flat_mean(:) };

plotColours = { [0 0 0] , [1 0 1] , [1 0.5 0] , [0 0 1] , [0 0.8 0.5] };
secLabel    = {'Whole Meander','Upper-Ridge Section','Lower-Ridge Section',...
    'Downstream Section', 'Flat Region'};

% 1. Prepare numeric time axis -----------------------------------------
if isdatetime(meand_time)
    timeNum = datenum(meand_time);          % serial date numbers
else
    timeNum = meand_time(:);                % already numeric
end
n         = numel(timeNum);
monthIdx  = (0:n-1)';                       % 0 … 364 (months from start)

% 2. Plot raw time series and capture handles --------------------------
figure('Color','w','Position',[120 120 900 500]); hold on;

hRaw = gobjects(numel(seriesList),1);       % preallocate for handles
for k = 1:numel(seriesList)
    hRaw(k) = plot(timeNum, seriesList{k}, '-', ...
                   'Color', plotColours{k}, 'LineWidth', 0.5);
end

xlabel('Year');  ylabel('Latitude (^{\circ}S)');
%title('Monthly Latitude Position of Meander Sections (1993–2023)');
lgd = legend(hRaw, secLabel,'Location','best');   % legend from raw handles
lgd.AutoUpdate = 'off';                           % freeze legend
grid on;

ax      = gca;
ax.XLim  = [timeNum(1) timeNum(end)];
ax.XTick = timeNum(1):365*5:timeNum(end);         % ~5-year spacing
datetick('x','yyyy','keepticks','keeplimits');
ax.FontSize = 9;

% 3. Loop over the five section series:  Sen slope, OLS trend, MK tests
% -------------------------------------------------------------------------
fprintf('\n=======  Trend statistics (1993–2023)  =======\n');
fprintf('Section               Sen-slope   MK-p  MKmod-p   OLS-slope   R²\n');

for k = 1:numel(seriesList)

    yFull = seriesList{k}(:);                % ensure column
    C     = plotColours{k};

    good  = isfinite(yFull);
    x     = monthIdx(good);                  % time in months (column)
    y     = yFull(good);

    % ---------- Sen slope ------------------------------------------------
    senSlope_mo = Sen_Slope(y);              % ° lat month⁻¹
    senSlope_dec = senSlope_mo * 120;        % ° lat decade⁻¹

    % ---------- Mann-Kendall --------------------------------------------
    [~,pMK]  = Mann_Kendall(y);
    [~,pMKm] = Modified_Mann_Kendall(y);

    % ---------- OLS slope + R² ------------------------------------------
    xC   = x - mean(x);                      % demean X for numerical stability
    beta_yr  = (xC' * (y - mean(y))) / (xC' * xC);     % ° lat month⁻¹
    beta_dec = beta_yr * 120;                            % per decade

    yFit = beta_yr * xC + mean(y);
    RSS  = sum((y - yFit).^2);
    TSS  = sum((y - mean(y)).^2);
    R2   = 1 - RSS/TSS;

    % ---------- Plot Sen trend line -------------------------------------
    b0      = median(yFull,'omitnan') - senSlope_mo*median(monthIdx,'omitnan');
    senLine = b0 + senSlope_mo*monthIdx;
    plot(timeNum, senLine, '--', 'Color',C, 'LineWidth',1.2);

    % ---------- Annotation ----------------------------------------------
    xText = timeNum(end) - 180;            % ~6 months from right edge
    yText = senLine(end-6);
    txt   = sprintf('%.2f^{\\circ}/dec  p=%.3f', senSlope_dec, pMK);
    text(xText, yText, txt, 'Color',C,'FontSize',7, ...
         'Interpreter','tex','Horiz','left','Vert','middle');

    % ---------- Console summary -----------------------------------------
    fprintf('%-20s %+7.2f  %6.3f  %7.3f   %+7.2f     %.3f\n', ...
            secLabel{k}, senSlope_dec, pMK, pMKm, beta_dec, R2);

end





%% Next, we continue to analyse the meridional movements of peaks and troughs
% of the meander

meand_posi_pk1_mean = nanmean(meand_posi_pk1, 1);
meand_posi_pk1_mean  = squeeze(meand_posi_pk1_mean).';

meand_posi_pk2_mean = nanmean(meand_posi_pk2, 1);
meand_posi_pk2_mean  = squeeze(meand_posi_pk2_mean).';

meand_posi_pk3_mean = nanmean(meand_posi_pk3, 1);
meand_posi_pk3_mean  = squeeze(meand_posi_pk3_mean).';


meand_posi_tr1_mean = nanmean(meand_posi_tr1, 1);
meand_posi_tr1_mean  = squeeze(meand_posi_tr1_mean).';

meand_posi_tr2_mean = nanmean(meand_posi_tr2, 1);
meand_posi_tr2_mean  = squeeze(meand_posi_tr2_mean).';

meand_posi_tr3_mean = nanmean(meand_posi_tr3, 1);
meand_posi_tr3_mean  = squeeze(meand_posi_tr3_mean).';


% Initial plot of the time series of the latitude positions of the
% meander's different sections
figure('Color','w');
hold on;

plot(meand_time, meand_posi_pk1_mean, '-k', 'LineWidth',0.5);
plot(meand_time, meand_posi_pk2_mean, '-m', 'LineWidth',0.5);
plot(meand_time, meand_posi_pk3_mean, '-', 'Color', '[1 0.5 0]','LineWidth',0.5);

plot(meand_time, meand_posi_tr1_mean, '-b', 'LineWidth',0.5);
plot(meand_time, meand_posi_tr2_mean, '-c', 'LineWidth',0.5);
plot(meand_time, meand_posi_tr3_mean, '-', 'Color', '[1 0.5 0.5]','LineWidth',0.5);

xlabel('Year');
ylabel('Latitude (^{\circ}S)');
%title('Meander Geostrophic-Speed Trend (1993 – 2023)');
legend({'Peak 1','Peak 2','Peak 3','Trough 1', 'Trough 2', 'Trough 3'},...
    'Location','best');
grid on;
xlim([meand_time(1) meand_time(end)]);  % (optional) full x-axis extent

set(gca,'FontSize',9,'XTick',meand_time(1):365:meand_time(end));
datetick('x','yyyy','keepticks','keeplimits');



%% ----------------------------------------------------------------------
%  INPUT (1×365)  — already averaged and squeezed above
%    meand_time               (datetime or numeric serial)
%    meand_posi_pk1_mean  …   meand_posi_tr3_mean
% ----------------------------------------------------------------------

% 0.  Pack series, colours, labels --------------------------------------
seriesList = { meand_posi_pk1_mean(:), meand_posi_pk2_mean(:), ...
               meand_posi_pk3_mean(:), meand_posi_tr1_mean(:), ...
               meand_posi_tr2_mean(:), meand_posi_tr3_mean(:) };

plotColours = { ...
    [0.0 0.0 0.0] ,   % black         – Peak 1
    [0.00 0.45 0.70] ,% blue          – Peak 2
    [0.80 0.47 0.00] ,% orange        – Peak 3
    [0.35 0.70 0.90] ,% sky-blue      – Trough 1
    [0.00 0.60 0.50] ,% bluish-green  – Trough 2
    [0.80 0.60 0.70]   % reddish-purple– Trough 3
};

secLabel = { 'Peak 1', 'Peak 2', 'Peak 3', ...
             'Trough 1', 'Trough 2', 'Trough 3' };

% 1.  Prepare numeric time axis ----------------------------------------
if isdatetime(meand_time)
    timeNum = datenum(meand_time);     % serial date numbers
else
    timeNum = meand_time(:);
end
nMonths = numel(timeNum);
monthIdx = (0:nMonths-1)';             % 0 … 364

% 2.  Plot raw curves and freeze legend --------------------------------
figure('Color','w','Position',[150 150 900 500]); hold on;
hRaw = gobjects(numel(seriesList),1);

for k = 1:numel(seriesList)
    hRaw(k) = plot(timeNum, seriesList{k}, '-', ...
                   'Color', plotColours{k}, 'LineWidth',0.5);
end

xlabel('Year');  ylabel('Latitude (^{\circ}S)');
%title('Monthly Latitude Positions of Peaks and Troughs (1993–2023)');
lgd = legend(hRaw, secLabel,'Location','best');
lgd.AutoUpdate = 'off';        % trends will not enter legend
grid on;

ax = gca;
ax.XLim  = [timeNum(1) timeNum(end)];
ax.XTick = timeNum(1):365*5:timeNum(end);   % ~5-year spacing
datetick('x','yyyy','keepticks','keeplimits');
ax.FontSize = 9;

% 3.  Sen slope, trend, annotation, console ----------------------------
fprintf('\n=======  Peaks & Troughs — Sen-slope / OLS statistics  =======\n');
fprintf('Series      Sen-slope  pMK  pMKmod   OLS-slope     R²\n');

for k = 1:numel(seriesList)

    % -------- raw series & colour ---------------------------------------
    y = seriesList{k}(:);       % column vector
    C = plotColours{k};

    good = isfinite(y);
    t    = monthIdx(good);      t = t(:);          % months since Jan-1993
    yG   = y(good);             yG = yG(:);

    % -------- Sen-slope & Mann–Kendall ----------------------------------
    senSlope_mo = Sen_Slope(yG);         % ° lat month⁻¹
    senSlope_dec = senSlope_mo * 120;    % ° lat decade⁻¹
    [~,pMK]  = Mann_Kendall(yG);
    [~,pMKm] = Modified_Mann_Kendall(yG);

    % -------- OLS slope & R² (for diagnostic only) ----------------------
    tC   = t - mean(t);                               % demean
    beta_mo = (tC' * (yG - mean(yG))) / (tC' * tC);   % ° lat month⁻¹
    beta_dec = beta_mo * 120;
    yFit = beta_mo * tC + mean(yG);
    R2   = 1 - sum((yG - yFit).^2) / sum((yG - mean(yG)).^2);

    % -------- Plot Sen trend line ---------------------------------------
    b0     = nanmedian(y) - senSlope_mo * nanmedian(monthIdx);
    senLine = b0 + senSlope_mo * monthIdx;
    plot(timeNum, senLine, '--', 'Color',C, 'LineWidth',1.2);

    % -------- Figure annotation -----------------------------------------
    txt = sprintf('%.2f^{\\circ}/dec  p=%.3f', senSlope_dec, pMK);
    xTxt = timeNum(end) - 180;       % ~6 months from right edge
    yTxt = senLine(end-6);
    text(xTxt, yTxt, txt, 'Color',C, 'FontSize',9, ...
         'Interpreter','tex', 'Horiz','left','Vert','middle');

    % -------- Console summary -------------------------------------------
    fprintf('%-11s %+7.2f  %.3f  %.3f   %+8.2f     %.3f\n', ...
            secLabel{k}, senSlope_dec, pMK, pMKm, beta_dec, R2);
end





%% Then, we continue to analyse the meander's width time series and their
% corresponding linear trends

% Below, we first divide the entire meander's width and speed into the five
% sections of the whole meander

% Firstly, width
meand_width_all = meand_width;
meand_width_upp = meand_width(1:60, :);
meand_width_low = meand_width(61:163, :);
meand_width_down = meand_width(164:end, :);
meand_width_flat = meand_width(170:215, :);


meand_width_all_mean = nanmean(meand_width_all, 1);
meand_width_all_mean  = squeeze(meand_width_all_mean).';

meand_width_upp_mean = nanmean(meand_width_upp, 1);
meand_width_upp_mean  = squeeze(meand_width_upp_mean).';

meand_width_low_mean = nanmean(meand_width_low, 1);
meand_width_low_mean  = squeeze(meand_width_low_mean).';

meand_width_down_mean = nanmean(meand_width_down, 1);
meand_width_down_mean  = squeeze(meand_width_down_mean).';

meand_width_flat_mean = nanmean(meand_width_flat, 1);
meand_width_flat_mean  = squeeze(meand_width_flat_mean).';



% Initial plot of the time series of the widths of the
% meander's different sections
figure('Color','w');
hold on;

plot(meand_time, meand_width_all_mean, '-k', 'LineWidth',0.5);
plot(meand_time, meand_width_upp_mean, '-m', 'LineWidth',0.5);
plot(meand_time, meand_width_low_mean, '-', 'Color', '[1 0.5 0]','LineWidth',0.5);
plot(meand_time, meand_width_down_mean, '-b', 'LineWidth',0.5);
plot(meand_time, meand_width_flat_mean, '-c', 'LineWidth',0.5);

xlabel('Year');
ylabel('Width (km)');
%title('Meander Location Trend (1993 – 2023)');
legend({'Whole Meander','Upper-Ridge Section','Lower-Ridge Section','Downstream Section', 'Flat Region'},...
    'Location','best');
grid on;
xlim([meand_time(1) meand_time(end)]);  % (optional) full x-axis extent

set(gca,'FontSize',9,'XTick',meand_time(1):365:meand_time(end));
datetick('x','yyyy','keepticks','keeplimits');



% The second step, we deal with the speed part of the meander's different
% sections
meand_loc_speed_all = meand_loc_speed;
meand_loc_speed_upp = meand_loc_speed(1:60, :);
meand_loc_speed_low = meand_loc_speed(61:163, :);
meand_loc_speed_down = meand_loc_speed(164:end, :);
meand_loc_speed_flat = meand_loc_speed(170:215, :);


meand_loc_speed_all_mean = nanmean(meand_loc_speed_all, 1);
meand_loc_speed_all_mean  = squeeze(meand_loc_speed_all_mean).';

meand_loc_speed_upp_mean = nanmean(meand_loc_speed_upp, 1);
meand_loc_speed_upp_mean  = squeeze(meand_loc_speed_upp_mean).';

meand_loc_speed_low_mean = nanmean(meand_loc_speed_low, 1);
meand_loc_speed_low_mean  = squeeze(meand_loc_speed_low_mean).';

meand_loc_speed_down_mean = nanmean(meand_loc_speed_down, 1);
meand_loc_speed_down_mean  = squeeze(meand_loc_speed_down_mean).';

meand_loc_speed_flat_mean = nanmean(meand_loc_speed_flat, 1);
meand_loc_speed_flat_mean  = squeeze(meand_loc_speed_flat_mean).';



% Initial plot of the time series of the widths of the
% meander's different sections
figure('Color','w');
hold on;

plot(meand_time, meand_loc_speed_all_mean, '-k', 'LineWidth',0.5);
plot(meand_time, meand_loc_speed_upp_mean, '-m', 'LineWidth',0.5);
plot(meand_time, meand_loc_speed_low_mean, '-', 'Color', '[1 0.5 0]','LineWidth',0.5);
plot(meand_time, meand_loc_speed_down_mean, '-b', 'LineWidth',0.5);
plot(meand_time, meand_loc_speed_flat_mean, '-c', 'LineWidth',0.5);

xlabel('Year');
ylabel('Speed (m s^{-1})');
%title('Meander Location Trend (1993 – 2023)');
legend({'Whole Meander','Upper-Ridge Section','Lower-Ridge Section','Downstream Section', 'Flat Region'},...
    'Location','best');
grid on;
xlim([meand_time(1) meand_time(end)]);  % (optional) full x-axis extent

set(gca,'FontSize',9,'XTick',meand_time(1):365:meand_time(end));
datetick('x','yyyy','keepticks','keeplimits');



%% ======================================================================
%  Common time axis ------------------------------------------------------
if isdatetime(meand_time)
    timeNum = datenum(meand_time);      % serial numbers for plotting/text()
else
    timeNum = meand_time(:);
end
n         = numel(timeNum);
monthIdx  = (0:n-1)';                   % 0 … 364 months from start


% ==================  WIDTH FIGURE  ====================================
fprintf('\n=== WIDTH (km) ===\n');
widthList  = { meand_width_all_mean(:),  meand_width_upp_mean(:), ...
               meand_width_low_mean(:), meand_width_down_mean(:), ...
               meand_width_flat_mean(:) };

plotColours = { [0 0 0] , [1 0 1] , [1 0.5 0] , [0 0 1] , [0 1 1] };
plotLabels  = { 'Whole','Upper','Lower','Downstream','Flat' };

plotWithTrend(timeNum, monthIdx, widthList, plotColours, plotLabels, ...
              'Width (km)',' km/dec');

% ==================  SPEED FIGURE  ====================================
fprintf('\n=== SPEED (m s^-1) ===\n');
speedList = { meand_loc_speed_all_mean(:),  meand_loc_speed_upp_mean(:), ...
              meand_loc_speed_low_mean(:), meand_loc_speed_down_mean(:), ...
              meand_loc_speed_flat_mean(:) };

plotWithTrend(timeNum, monthIdx, speedList, plotColours, plotLabels, ...
              'Speed (m s^{-1})',' m s^{-1}/dec');





%% Finally, we put together the three subplots of latitude position,
% width and geostrophic current speed in one figure together with the
% linear temporal trend analysis

% ========================  INPUT SERIES  =============================
%  meand_time (datetime or datenum 1×365)
%  Latitude:  meand_loc_lat_*_mean
%  Width   :  meand_width_*_mean
%  Speed   :  meand_loc_speed_*_mean
% ======================================================================

% ----------------------------------------------------------------------
%  1. Pack the three data “blocks” into cell arrays ---------------------
latSeries   = { meand_loc_lat_all_mean(:),  meand_loc_lat_upp_mean(:), ...
                meand_loc_lat_low_mean(:),  meand_loc_lat_down_mean(:), ...
                meand_loc_lat_flat_mean(:) };

widSeries   = { meand_width_all_mean(:)*100,    meand_width_upp_mean(:)*100, ...
                meand_width_low_mean(:)*100,    meand_width_down_mean(:)*100, ...
                meand_width_flat_mean(:)*100 };

spdSeries   = { meand_loc_speed_all_mean(:), meand_loc_speed_upp_mean(:), ...
                meand_loc_speed_low_mean(:), meand_loc_speed_down_mean(:), ...
                meand_loc_speed_flat_mean(:) };

blockNames  = {'Latitude (^{\circ}S)', 'Width (km)', 'Speed (m s^{-1})'};
blocks      = {latSeries, widSeries, spdSeries};
unitsTxt    = {' ^\circ/dec', ' km/dec', ' m s^{-1}/dec'};

plotCols = { [0 0 0], [1 0 1], [1 0.5 0], [0 0 1], [0 0.8 0.5] };
secLabels = {'Whole Meander','Upper-Ridge Section','Lower-Ridge Section',...
    'Downstream Section','Flat Region'};

% ----------------------------------------------------------------------
%  2. Build numeric time axis ------------------------------------------
if isdatetime(meand_time)
    tNum = datenum(meand_time);             % serial days
else
    tNum = meand_time(:);                   % already numeric
end
nMonths  = numel(tNum);
mIdx     = (0:nMonths-1)';                  % 0 … 364

% 3.  Prepare stats-output table ---------------------------------------
statsTbl = table('Size',[5 9],'VariableTypes',repmat("double",1,9), ...
                 'VariableNames',{'LatSlope','Lat_pMK','Lat_pMKm', ...
                                  'WidSlope','Wid_pMK','Wid_pMKm', ...
                                  'SpdSlope','Spd_pMK','Spd_pMKm'}, ...
                 'RowNames',secLabels);

% 4.  Plot figure with three sub-plots ---------------------------------
figure('Color','w','Position',[120 60 900 1100]);

statsTbl = cell(6,13);                      % 5 sections + header
statsTbl(1,:) = {'Section', ...
                 'Lat  Sen','Lat  pMK','Lat  pMKm','Lat  R2', ...
                 'Wid  Sen','Wid  pMK','Wid  pMKm','Wid  R2', ...
                 'Spd  Sen','Spd  pMK','Spd  pMKm','Spd  R2'};

for blk = 1:3
    subplot(3,1,blk); hold on
    rawHandles = gobjects(5,1);

    for s = 1:5
        y = blocks{blk}{s}(:);                 % ensure column
        col = plotCols{s};

        % ---------- raw plot -------------------------------------------
        rawHandles(s) = plot(tNum, y,'-','Color',col,'LineWidth',0.5);

        % ---------- Sen slope & MK -------------------------------------
        good = isfinite(y);
        senSlope = Sen_Slope(y(good));         % units per month
        senDec   = senSlope * 120;
        [~,pMK]  = Mann_Kendall(y(good));
        [~,pMKm] = Modified_Mann_Kendall(y(good));

        % ---------- OLS slope & R²  ------------------------------------
        xG = mIdx(good); xC = xG - mean(xG);
        yC = y(good)  - mean(y(good));
        beta = (xC' * yC) / (xC' * xC);        % per month
        yFit = beta * xC;
        R2   = 1 - sum((yC - yFit).^2) / sum(yC.^2);

        % ---------- plot Sen line --------------------------------------
        b0    = nanmedian(y) - senSlope*nanmedian(mIdx);
        trend = b0 + senSlope*mIdx;
        plot(tNum, trend, '--','Color',col,'LineWidth',1.0);

        % ---------- annotate Sen slope ---------------------------------
        xTxt = tNum(end) - 180;
        yTxt = trend(end-6);
        text(xTxt, yTxt, sprintf('%.3f%s', senDec, unitsTxt{blk}), ...
             'Color',col,'FontSize',9,'Interpreter','tex', ...
             'Horiz','left','Vert','middle');

        % ---------- dump stats into table ------------------------------
        switch blk          % choose column offset
            case 1, c0 = 2;   % latitude stats start col 2
            case 2, c0 = 6;   % width stats   start col 6
            case 3, c0 = 10;  % speed stats   start col 10
        end
        statsTbl{s+1,1}   = secLabels{s};          % section name
        statsTbl{s+1,c0}  = senDec;
        statsTbl{s+1,c0+1}= pMK;
        statsTbl{s+1,c0+2}= pMKm;
        statsTbl{s+1,c0+3}= R2;
    end

    % ---------- axes cosmetics -----------------------------------------
    ylabel(blockNames{blk},'Interpreter','tex');
    grid on
    ax = gca;
    ax.XLim  = [tNum(1) tNum(end)];
    ax.XTick = tNum(1):365*5:tNum(end);
    datetick('x','yyyy','keepticks','keeplimits');
    ax.FontSize = 9;
    if blk==1
        legend(rawHandles,secLabels,'Location','best','AutoUpdate','off');
    end
    if blk<3, ax.XTickLabel = ''; else, xlabel('Year'); end
end

%sgtitle('Meander Evolution: Latitude, Width & Speed (1993–2023)','FontWeight','normal');

% --- build and display a table -----------------------------------------
statsMat  = cell2mat(statsTbl(2:end,2:end));           % numeric 5 × 12
varNames  = matlab.lang.makeValidName(statsTbl(1,2:end));
rowNames  = statsTbl(2:end,1);

T = array2table(statsMat, ...
                'VariableNames', varNames, ...
                'RowNames',     rowNames);

disp(T);





%% Figure 6: monthly time series of the meander's amplitude composed of pk 3 and tr 2
% and its corresponding linear trends over the 1993--2023 period
meand_wave_ampli = (meand_posi_tr2_mean - meand_posi_pk3_mean) / 2;

figure('Color','w');
hold on;

plot(meand_time, meand_wave_ampli, '-k', 'LineWidth',0.5);

xlabel('Year');
ylabel('Wave Amplitude (^{\circ}latitude)');
%title('Meander Location Trend (1993 – 2023)');
%legend({'Wave'},'Location','best');
grid on;
xlim([meand_time(1) meand_time(end)]);  % (optional) full x-axis extent

set(gca,'FontSize',9,'XTick',meand_time(1):365:meand_time(end));
datetick('x','yyyy','keepticks','keeplimits');



% 1.  Remove missing values and build demeaned time axis
good  = isfinite(waveAmp);
yG    = waveAmp(good);                       % cleaned series
xG    = mIdx(good);                          % corresponding months
xC    = xG - mean(xG);                       % demeaned for OLS

% 2.  Robust Sen slope and non-parametric tests
senSlope_mo  = Sen_Slope(yG);                % ° lat mo⁻¹
senSlope_dec = senSlope_mo * 120;            % ° lat dec⁻¹
[~,pMK]   = Mann_Kendall(yG);                % classical M–K
[~,pMKm]  = Modified_Mann_Kendall(yG);       % autocorr-adjusted M–K

% 3.  Ordinary-least-squares (diagnostic) slope & R²
beta_mo  = (xC' * (yG - mean(yG))) / (xC' * xC);   % ° lat mo⁻¹
beta_dec = beta_mo * 120;                          % ° lat dec⁻¹
yFit     = beta_mo * xC + mean(yG);
R2       = 1 - sum((yG - yFit).^2) / sum((yG - mean(yG)).^2);

% 4.  Sen trend line for plotting (use full time axis)
b0      = median(yG) - senSlope_mo * median(xG);
senLine = b0 + senSlope_mo * mIdx;           % for *all* months

% 5.  Plot
figure('Color','w'); hold on
plot(tNum, waveAmp, '-k', 'LineWidth',0.5);        % raw series
plot(tNum, senLine, '--k','LineWidth',1.2);        % Sen trend

xlabel('Year');  ylabel('Wave amplitude (^{\circ} latitude)');
grid on
ax = gca;
ax.XLim  = [tNum(1) tNum(end)];
ax.XTick = tNum(1):365*5:tNum(end);
datetick('x','yyyy','keepticks','keeplimits');
ax.FontSize = 9;

% Annotation (≈ six months before record end)
xTxt = tNum(end) - 180;
yTxt = senLine(end-6);
txt  = sprintf('%.2f^{\\circ}/dec  p=%.3f  R^{2}=%.2f', ...
               senSlope_dec, pMK, R2);
text(xTxt, yTxt, txt, 'Color','k','FontSize',8, ...
     'Interpreter','tex','Horiz','left','Vert','middle');

% 6.  Console summary
fprintf('\n=== Wave amplitude trend  (Pk-3 – Tr-2)/2  ===\n');
fprintf('Sen slope   : %+8.2f °/dec\n',  senSlope_dec);
fprintf('OLS slope   : %+8.2f °/dec\n',  beta_dec);
fprintf('R² (OLS)    : %.3f\n',           R2);
fprintf('M-K  p-value: %.3f\n',           pMK);
fprintf('Mod M-K p   : %.3f\n',           pMKm);



%% Figure 7 and 8: Oceanic thermal contour movements and structure changes

% Figure 7: ocean temperatures contours latitude x depth figure

%  SECTION-PLOT OF MEAN TEMPERATURE (isotherms)
%  Three equal 6-yr periods: 2006-11, 2012-17, 2018-23
%  Domain: Pacific–Antarctic Ridge (60° S–48° S , 150° W–80° W)
%  ----------------------------------------------------------------------
%  Inputs already in the workspace
%     •  meand_temp_04_23_4D_TIME_LON_LAT_PRES   % TIME × LON × LAT × PRESS
%     •  meand_time_mon_04_23                    % TIME
%     •  latSub   (12 × 1)   ,  lonSub (70 × 1)  ,  prs (58 × 1)
% ======================================================================

Tfull = meand_temp_04_23_4D_TIME_LON_LAT_PRES;
time  = meand_time_mon_04_23;          % datetime vector

% 1.  Keep months from Jan-2006 onward
idx06 = find(time >= datetime(2006,1,1));
Tfull = Tfull(idx06,:,:,:);
time  = time(idx06);


% 2.  Define three 6-yr periods (exactly the same length)
edges  = datetime([2006 01 01; 2012 01 01; 2018 01 01; 2024 01 01]);
labels = {'2006–2011','2012–2017','2018–2023'};
nPer   = 3;
idxPer = cell(nPer,1);
for p = 1:nPer
    idxPer{p} = find(time >= edges(p) & time < edges(p+1));
end


% 3.  Compute mean field for each period (PRESS × LAT)
meanField = cell(1,nPer);                    % each 58 × 12

for p = 1:nPer
    Tsel = Tfull(idxPer{p},:,:,:);           % TIME × LON × LAT × PRESS
    % mean over TIME (dim-1) and LON (dim-2)
    tmp  = squeeze( mean( mean(Tsel,1,'omitnan'), 2,'omitnan') ); % LAT × PRESS
    meanField{p} = tmp.';                    % → PRESS × LAT
end


%  Make sure latitude/longitude subset vectors exist
% ----------------------------------------------------------
if ~exist('latSub','var') || ~exist('lonSub','var')
    % full coordinate vectors lon, lat and subset indices ilon, ilat
    % were created earlier when you loaded the NetCDF file
    latSub = lat(ilat);     % 12 × 1   (−60 … −48)
    lonSub = lon(ilon);     % 70 × 1   (−150 … −80)
end


% 4.  Plot (isotherms)
figure('Color','w','Position',[50 100 1000 600]);

% colour-blind-safe palette
clr = {[0 0 0.7], [0.65 0 0.65], [0 0.6 0.6]};  % navy, purple, teal
ls  = {'--','--','--'};  lw = 1.5;

hold on

for p = 1:nPer
    Z = meanField{p};                         % 58 × 12
    [C,h] = contour(latSub, prs, Z, ...
                    'LevelStep',1, ...
                    'LineColor',clr{p}, ...
                    'LineStyle',ls{p}, ...
                    'LineWidth',lw);
    if p==2                                    % annotate middle period
        clabel(C, h, ...
       'Color',       'r', ...
       'FontSize',    10, ...
       'FontWeight', 'bold', ...
       'LabelSpacing', 250);   % ← numeric value only
    end
    phndl(p) = h(1);                          % store handle for legend
end
hold off

set(gca,'YDir','reverse','FontSize',10,...
        'XTick',-60:2:-48,'YTick',0:200:max(prs));
xlabel('Latitude (^{\circ})');
ylabel('Pressure (dbar)');
title('Whole Meander: temperature (^{\circ}C)');

legend(phndl,labels,'Location','southeast','Box','off');


%% Below, we continue to plot the same spatial structure of the figure but
% in the four different sections

% 1.  Trim to Jan-2006 onward (if not already)
keep = time >= datetime(2006,1,1);
Tfull = Tfull(keep,:,:,:);
time  = time(keep);


% 2.  Period indices (exactly 6 yr each)
edges  = datetime([2006 01 01; 2012 01 01; 2018 01 01; 2024 01 01]);
labels = {'2006–11','2012–17','2018–23'};
nPer   = 3;
idxPer = cell(nPer,1);
for p = 1:nPer
    idxPer{p} = find(time >= edges(p) & time < edges(p+1));
end


% 3.  Longitude masks for four meander sections (all W longitudes)
maskUpper = lonSub >= -150   & lonSub <= -135.1;       % 135.1–150°W
maskLower = lonSub >  -135.1 & lonSub <= -109.4;       % 109.4–135.1°W
maskDown  = lonSub >  -109.4 & lonSub <=  -80;         %  80.0–109.4°W
maskFlat  = lonSub >= -107.6 & lonSub <= -96.4;        %  96.4–107.6°W

masks  = {maskUpper, maskLower, maskDown, maskFlat};
titles = {'(a) Upper-Ridge Section', '(b) Lower-Ridge Section',...
    '(c) Downstream Section', '(d) Flat Region'};


% 4.  Compute mean temperature (PRESS × LAT) for each section/period
meanField = cell(numel(masks), nPer);                  % 58 × 12 each

for s = 1:numel(masks)
    for p = 1:nPer
        Tsel = Tfull(idxPer{p}, masks{s}, :, :);       % TIME × LON × LAT × PRESS
        Z    = squeeze( mean( mean(Tsel, 2,'omitnan'), 1,'omitnan') ); % LAT × PRESS
        meanField{s,p} = Z.';                          % → PRESS × LAT
    end
end


% 5.  Plot settings
%clr = {[0 0 0.7], [0.55 0 0.55], [0 0.6 0.6]};           % colour-blind safe
%ls  = {'--','--','--'};
%lw = 1.4;
%latTicks = -60:2:-48;
%prsTicks = 0:200:max(prs);


% 6.  Figure with 4 stacked panels
figure('Color','w','Position',[50 100 680 780]);

latTicks = -60:2:-48;     prsTicks = 0:200:max(prs);
clr  = {[0 0 0.70],[0.55 0 0.55],[0 0.60 0.60]};   % navy, violet, teal
ls   = {'--','--','--'};  lw = 1.4;

legH   = gobjects(1,3);    % legend handles (one per 6-yr period)
legLbl = labels;           % {'2006–11','2012–17','2018–23'}

for s = 1:numel(masks)
    ax = subplot(4,1,s);  hold(ax,'on');

    for p = 1:nPer
        Z = meanField{s,p};
        if all(~isfinite(Z(:))), continue, end        % nothing to plot

        % draw contours (single level if constant; otherwise 1 °C steps)
        if range(Z(:),'omitnan') == 0
            lvl = Z(find(isfinite(Z),1));
            [C,h] = contour(ax, latSub, prs, Z, ...
                            'LevelList',lvl, ...
                            'LineColor',clr{p},'LineStyle',ls{p}, ...
                            'LineWidth',lw);
        else
            [C,h] = contour(ax, latSub, prs, Z, ...
                            'LevelStep',1, ...
                            'LineColor',clr{p},'LineStyle',ls{p}, ...
                            'LineWidth',lw);
            if p==2
                clabel(C,h,'Color','r','FontSize',10, ...
                           'FontWeight','bold','LabelSpacing',250);
            end
        end

        % store first valid handle for this period (any panel)
        if isempty(legH(p)) && ishghandle(h(1))
            legH(p) = h(1);
        end
    end

    hold(ax,'off');
    set(ax,'YDir','reverse','XTick',latTicks,'YTick',prsTicks,'FontSize',9);
    ylabel(ax,'Pressure (dbar)');
    if s==numel(masks)
        xlabel(ax,'Latitude (^{\circ})');
    else
        ax.XTickLabel = [];
    end
    title(ax,titles{s},'FontWeight','normal');
    if exist('cmocean','file'), colormap(ax,cmocean('balance')); end
end

% ---------- build legend from collected handles ------------------------
valid = isgraphics(legH);                 % which periods produced lines
if any(valid)
    % place the legend below the bottom panel (southoutside)
    legend(legH(valid), legLbl(valid), ...
           'Orientation','horizontal', ...
           'Location','southoutside', ...
           'Box','off','FontSize',9);
end

%sgtitle('Pacific–Antarctic Ridge: mean isotherms for three 6-year periods');



%% Figure 8: latitude x depth structure of the annual-mean temperautre trends
% in the whole meander an its four different sections

%  0.  PREP:  restrict to Jan-2006 – Jun-2023
% ======================================================================
keep  = time >= datetime(2006,1,1) & time <= datetime(2023,6,30);
Tcube = double(meand_temp_04_23_4D_TIME_LON_LAT_PRES(keep,:,:,:));
time  = time(keep);

tyr = year(time) + (day(time,'dayofyear')-1)/365;     % decimal years
tCol = tyr(:);    tD = tCol - mean(tCol);             % demean
Nt   = numel(tCol);
Sxx  = sum(tD.^2);


% 1.  LONGITUDE MASKS
maskWhole = true(size(lonSub));
maskUpper = lonSub >= -150   & lonSub <= -135.1;   % 135.1–150°W
maskLower = lonSub >  -135.1 & lonSub <= -109.4;   % 109.4–135.1°W
maskDown  = lonSub >  -109.4 & lonSub <=  -80;     %  80.0–109.4°W
maskFlat  = lonSub >= -107.6 & lonSub <= -96.4;    %  96.4–107.6°W

masks  = {maskWhole,maskUpper,maskLower,maskDown,maskFlat};
titles = {'(a) Whole Meander','(b) Upper-Ridge Section','(c) Lower-Ridge Section',...
          '(d) Downstream Section','(e) Flat Region'};

nSec = numel(masks);  nLat = numel(latSub);  nPrs = numel(prs);
trend   = cell(1,nSec);      % slope (°C yr⁻¹)
sigMask = cell(1,nSec);      % p < 0.05


% 2.  TREND + SIGNIFICANCE
for s = 1:nSec
    slope  = NaN(nPrs,nLat);
    sig    = false(nPrs,nLat);

    % average over longitude for this section -> TIME × LAT × PRESS
    Tlon = squeeze( mean(Tcube(:,masks{s},:,:), 2, 'omitnan') );

    for j = 1:nLat
        for i = 1:nPrs
            y = Tlon(:,j,i);
            good = isfinite(y);
            if sum(good) < 4,  continue,  end
            y  = y(good);
            td = tCol(good) - mean(tCol(good));

            b  = sum(td .* (y - mean(y))) / sum(td.^2);
            yhat = b*td + mean(y);
            RSS = sum((y - yhat).^2);
            SEb = sqrt(RSS/(numel(y)-2) / sum(td.^2));
            tstat = b / SEb;
            pval  = 2*(1 - tcdf(abs(tstat), numel(y)-2));

            slope(i,j) = b;
            sig(i,j)   = pval < 0.05;
        end
    end
    trend{s}   = slope;
    sigMask{s} = sig;
end


% 3.  PLOT  (3 rows × 2 columns)
figure('Color','w','Position',[90 60 780 820]);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

cmap  = cmocean('balance',256);   % diverging ± palette
colormap(cmap);                  % apply once at figure level

clims    = [-0.05 0.05];         % common colour limits
latTicks = -60:2:-48;
prsTicks = 0:200:max(prs);

axList = gobjects(1,nSec);       % store axes handles for linking

for s = 1:nSec
    ax = nexttile(tl);   axList(s) = ax;  hold(ax,'on');

    % trend field (NaNs transparent)
    im = imagesc(ax, latSub, prs, trend{s});
    set(im,'AlphaData',~isnan(trend{s}));

    % significance dots
    [LATM,PRSM] = meshgrid(latSub,prs);
    scatter(ax, LATM(sigMask{s}), PRSM(sigMask{s}), ...
            12,'k','filled','MarkerFaceAlpha',0.75);

    % axes styling
    set(ax,'YDir','reverse','CLim',clims,...
           'FontSize',9,'XTick',latTicks,'YTick',prsTicks);
    if mod(s,2)==1, ylabel(ax,'Depth (m)'), end
    if s>4,        xlabel(ax,'Latitude (^{\circ})'), end
    title(ax,titles{s},'FontWeight','bold','FontSize',11);
    hold(ax,'off');
end

% --- link all axes to share the same CLim ------------------------------
linkprop(axList,{'CLim'});

% --- single colour-bar, outside the tiled layout -----------------------
cb = colorbar(axList(1), 'eastoutside');   % attach to an axes, not to tl
cb = colorbar('eastoutside');   % attach to an axes, not to tl
cb.Ticks          = -0.04:0.02:0.04;
cb.Label.String   = 'Temperature trend (°C yr^{-1})';
cb.Label.FontSize = 10;
cb.FontSize       = 9;

%title(tl,'Pacific–Antarctic Ridge • Annual temperature trends 2006–23', ...
      %'FontWeight','bold','FontSize',12);



%%

% 3.  PLOT  –  annual temperature‐trend sections (5 panels + 1 blank)
%     single colour-bar shared by all panels, placed at the far right
% ----------------------------------------------------------------------
figure('Color','w','Position',[90 60 780 820]);

tl = tiledlayout(3,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');        % leaves room for colour-bar

cmap  = cmocean('balance',256);
colormap(cmap);                      % apply once to the whole figure

clims    = [-0.05 0.05];             % °C yr⁻¹
latTicks = -60:2:-48;
prsTicks = 0:200:max(prs);

axList = gobjects(1,nSec);           % save axes for linking limits

for s = 1:nSec
    ax = nexttile(tl);   axList(s) = ax;  hold(ax,'on');

    % --- temperature-trend field (imagesc) -----------------------------
    im = imagesc(ax, latSub, prs, trend{s}, clims);
    set(im,'AlphaData',~isnan(trend{s}));      % NaNs transparent

    % --- overlay significance dots ------------------------------------
    [LATM,PRSM] = meshgrid(latSub,prs);
    scatter(ax, LATM(sigMask{s}), PRSM(sigMask{s}), ...
            12,'k','filled','MarkerFaceAlpha',0.75);

    % --- axes cosmetics -----------------------------------------------
    set(ax,'YDir','reverse', ...
           'FontSize',9,'XTick',latTicks,'YTick',prsTicks);
    if mod(s,2)==1, ylabel(ax,'Depth (m)'), end
    if s>4,        xlabel(ax,'Latitude (^{\circ})'), end
    title(ax,titles{s},'FontWeight','bold','FontSize',11);
    hold(ax,'off');
end

% ----------------------------------------------------------------------
% Link all CLim so one colour-bar represents every panel
% ----------------------------------------------------------------------
linkprop(axList,'CLim');             % shares colour limits

% ----------------------------------------------------------------------
% Single global colour-bar, east-outside of the tiled layout
% ----------------------------------------------------------------------
cb = colorbar(axList(1));          % create, referencing the first axis
cb.Location      = 'eastoutside';  % move it outside the tiled layout
cb.Ticks         = -0.04:0.02:0.04;
cb.Label.String  = 'Temperature trend (°C yr^{-1})';
cb.FontSize      = 9;
cb.Label.FontSize = 10;
