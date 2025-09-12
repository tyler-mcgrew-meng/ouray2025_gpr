%
clear all; close all;

%import tool
%double click on csv file, select numeric matrix option 
%select range of position row (row 2), and import selection, change variable name to line#_x
%slect range of time column (column B), and import selection, change
%variable name to line#_t


%Meng, T.M., Petersen, E.I. and Holt, J.W. (2023) “Rock glacier composition and structure from radio wave speed analysis with dipping reflector correction,” Journal of Glaciology, 69(275), pp. 639–657. Available at: https://doi.org/10.1017/jog.2022.90.

%Meng, T.M. et al. (2025) “Effects of Rock Glacier Dynamics on Surface Morphology and Deformation,” Journal of Geophysical Research: Earth Surface, 130(3), p. e2024JF008106. Available at: https://doi.org/10.1029/2024JF008106.


% load GPR and topography data 
load('csv/line14.mat');
% imo_long_x = 0:0.25:99; % 100 MHz

imo_long_x = line14_x; % 25 MHz

% imo_long_x = 0:96; % 25 MHz
imo_long_t = line14_t;
amplitude = line14;

imo_long_topo = rmmissing(readtable('topo/line14_topo.csv')); topo_z = (imo_long_topo.z);
topo_x = imo_long_topo.x;

% debris_interp = readtable('interp/line147-18_debris.csv', 'NumHeaderLines',3);
% base_interp = readtable('interp/line147-18_base.csv', 'NumHeaderLines',3);

%assumed or measured velocity (m/ns)
v1 = 0.14;

%$ Topographic correction using GPR data and a separate elevation profile 
%depth correction from time to depth space
imo_long_depth = 0.5*v1.*imo_long_t;

%vertical sample interval in meters
dz = abs(imo_long_depth(2) - imo_long_depth(1));

% apply manual power law gain function fo plotting
for i = 1:size(amplitude,1)
    imo_long_gain(i, :) = amplitude(i, :).*(0.2).*i.^(1.7);
end

%start at t = 0 
imo_long_gain = imo_long_gain(imo_long_depth >= 0,:);
imo_long_depth = imo_long_depth(imo_long_depth >= 0);

%number of samples per trace 
L = length(imo_long_depth);

%datum to highest elevation in profile
datum = max(topo_z);

%interpolate topography data to match GPR sample spacing 
topo_x = topo_x(2:length(topo_x)-1); topo_z = topo_z(2:length(topo_z)-1);
F = griddedInterpolant(topo_x, topo_z)
topo_interp = F(imo_long_x);

datum = max(topo_interp);

% calculate elevation offset in meters and array indices
shift = datum-topo_interp;
shift_num = ceil(shift./dz);

%create grid in elevation space
imo_long_elev = min(topo_z)- max(imo_long_depth)-max(shift):dz:max(topo_interp);

%shift from depth to elevation space 
imo_long_shift = zeros(length(imo_long_elev),length(imo_long_x));

for j = 1:length(imo_long_x) 
   imo_long_shift(shift_num(j)+1:shift_num(j)+L,j) = imo_long_gain(:,j);
end

%plots
figure(1); 
plot(topo_x, topo_z, 'o'); hold on;
plot(imo_long_x, topo_interp,'.');

figure(2);
pcolor(imo_long_x, imo_long_depth, imo_long_gain); hold on;
shading interp; 
% plot(base_interp.Position_m_, base_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
% plot(debris_interp.Position_m_, debris_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
ylim([0 20]);
caxis([-5e6 5e6]);
axis ij;
colormap bone;

figure(3);
pcolor(imo_long_x, imo_long_elev, flipud(imo_long_shift)); 
shading interp; hold on;
plot(imo_long_x, topo_interp-3,'k', 'LineWidth', 2);
ylim([3730 max(imo_long_elev)]);
% xlim([0 100]);
caxis([-5e6 5e6]);
colormap bone;
% cmocean('balance');
set(gca, 'FontSize', 20);
xlabel('Horizontal Position (m)','FontSize',30);
ylabel('Elevation (m)', 'FontSize', 30);
% title('Lines 17 & 18','FontSize',30);