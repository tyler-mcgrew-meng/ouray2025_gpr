%
clear all; close all;

%import tool
%double click on csv file, select numeric matrix option 
%select range of position row (row 2), and import selection, change variable name to line#_x
%slect range of time column (column B), and import selection, change
%variable name to line#_t


% load GPR and topography data 
load('csv/imo_long_1.mat');
load('csv/imo_long_2.mat');

% imo_long_x = 0:0.25:99; % 100 MHz


MergedLine2_x = MergedLine2_x + max(MergedLine_x);
% 
% figure(1);
% pcolor(MergedLine2_x,MergedLine2_t, MergedLine2); hold on; shading interp; axis ij; colormap bone;
% pcolor(MergedLine_x,MergedLine_t, MergedLine); shading interp;
% clim([-1e4 1e4]);
% 
% xlim([0 500]);

imo_long_x = horzcat(MergedLine_x,MergedLine2_x); % 25 MHz

% imo_long_x = 0:96; % 25 MHz
% imo_long_t = line1_t;
% amplitude = line1;

imo_long_topo = rmmissing(readtable('topo/imogene_long_topo_hires.csv')); topo_z = (imo_long_topo.z);
topo_x = imo_long_topo.x;

% debris_interp = readtable('interp/line17-18_debris.csv', 'NumHeaderLines',3);
% base_interp = readtable('interp/line17-18_base.csv', 'NumHeaderLines',3);

%assumed or measured velocity (m/ns)
v1 = 0.14;

%$ Topographic correction using GPR data and a separate elevation profile 
%depth correction from time to depth space
imo_long_depth1 = 0.5*v1.*MergedLine_t;
imo_long_depth2 = 0.5*v1.*MergedLine2_t;

%vertical sample interval in meters
dz1 = abs(imo_long_depth1(2) - imo_long_depth1(1));
dz2 = abs(imo_long_depth2(2) - imo_long_depth2(1));

% apply manual power law gain function fo plotting
for i = 1:size(MergedLine,1)
    imo_long_gain1(i, :) = MergedLine(i, :).*(0.2).*i.^(1.9);
end

for i = 1:size(MergedLine2,1)
    imo_long_gain2(i, :) = MergedLine2(i, :).*(0.3).*i.^(2);
end

%start at t = 0 
imo_long_gain1 = imo_long_gain1(imo_long_depth1 >= 0,:);
imo_long_depth1 = imo_long_depth1(imo_long_depth1 >= 0);

imo_long_gain2 = imo_long_gain2(imo_long_depth2 >= 0,:);
imo_long_depth2 = imo_long_depth2(imo_long_depth2 >= 0);

%number of samples per trace 
L1 = length(imo_long_depth1);
L2 = length(imo_long_depth2);

%datum to highest elevation in profile
datum = max(topo_z);

%interpolate topography data to match GPR sample spacing 
topo_x = topo_x(2:length(topo_x)-1); topo_z = topo_z(2:length(topo_z)-1);
F = griddedInterpolant(topo_x, topo_z)
topo_interp = F(imo_long_x);

datum = max(topo_interp);

% calculate elevation offset in meters and array indices
shift = datum-topo_interp;
shift_num1 = ceil(shift./dz1);
shift_num2 = ceil(shift./dz2);

%create grid in elevation space
imo_long_elev1 = min(topo_z)- max(imo_long_depth1)-max(shift):dz1:max(topo_interp);
imo_long_elev2 = min(topo_z)- max(imo_long_depth2)-max(shift):dz2:max(topo_interp);


%shift from depth to elevation space 
imo_shift1 = zeros(length(imo_long_elev1),length(MergedLine_x));

imo_shift2 = zeros(length(imo_long_elev2),length(MergedLine2_x));


for j = 1:length(MergedLine_x) 
   imo_shift1(shift_num1(j)+1:shift_num1(j)+L1,j) = imo_long_gain1(:,j);
end

for j = 1:length(MergedLine2_x) 
   imo_shift2(shift_num2(j+200)+1:shift_num2(j+200)+L2,j) = imo_long_gain2(:,j);
end

%plots
figure(1); 
plot(topo_x, topo_z, 'o'); hold on;
plot(imo_long_x, topo_interp,'.');

figure(2);
pcolor(MergedLine_x, imo_long_depth1, imo_long_gain1); hold on;
pcolor(MergedLine2_x, imo_long_depth2, imo_long_gain2);
shading interp; 
% plot(base_interp.Position_m_, base_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
% plot(debris_interp.Position_m_, debris_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
% ylim([0 20])
xlim([0 495]);
ylim([0 60]);
caxis([-5e6 5e6]);
axis ij;
colormap bone;

figure(3);
pcolor(MergedLine_x, imo_long_elev1, flipud(imo_shift1)); hold on;
pcolor(MergedLine2_x, imo_long_elev2, flipud(imo_shift2));
shading interp; 
plot(imo_long_x, topo_interp,'k', 'LineWidth', 2);
ylim([3450 max(topo_interp)]);
xlim([0 495]);
caxis([-5e6 5e6]);
colormap bone;
% cmocean('balance');
set(gca, 'FontSize', 20);
xlabel('Horizontal Position (m)','FontSize',30);
ylabel('Elevation (m)', 'FontSize', 30);
% title('Lines 17 & 18','FontSize',30);