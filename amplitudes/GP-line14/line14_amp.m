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
load('line14.mat');
% imo_long_x = 0:0.25:99; % 100 MHz

imo_long_x = line14_x; % 25 MHz

% imo_long_x = 0:96; % 25 MHz
imo_long_t = line14_t;
amplitude = line14;

imo_long_topo = rmmissing(readtable('line14_topo_hires.csv')); topo_z = (imo_long_topo.z);
topo_x = imo_long_topo.x;
topo_E = imo_long_topo.E;
topo_N = imo_long_topo.N;

% debris_interp = readtable('interp/line147-18_debris.csv', 'NumHeaderLines',3);
base_interp = readtable('line14_base-Point.csv', 'NumHeaderLines',3);

%assumed or measured velocity (m/ns)
v1 = 0.15;

v_int = base_interp.Velocity_m_ns_(1); %ensure velocities are consistent between radargram and interpretation

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
topo_E = topo_E(2:length(topo_E)-1); topo_N = topo_N(2:length(topo_N)-1);

%%%for hi-resolution topography, I had to modify these specific x values
%because of the parsing of the QGIS-exported topography profile, leading to
%repeated values at these specific indices, throwing an error when calling
%"griddedInterpolant"
topo_x(286:313) = 10.03:0.03:10.84;

F = griddedInterpolant(topo_x, topo_z);
topo_interp = F(imo_long_x);
% topo_interp = interp1(topo_x,topo_z,imo_long_x);
datum = max(topo_interp);

base_interpolate = interp1(base_interp.Position_m_, base_interp.Depth_m_/v_int*v1, topo_x);
base_trace = interp1(base_interp.Position_m_, base_interp.Time_ns_, imo_long_x);
base_depth = interp1(base_interp.Position_m_, base_interp.Depth_m_/v_int*v1, imo_long_x);



power = zeros(size(imo_long_x));

dt = imo_long_t(2)-imo_long_t(1);

win_size = 20;

for i = 3:length(imo_long_x)-2
    window = imo_long_t(base_trace(i)-imo_long_t < win_size & base_trace(i)-imo_long_t > - win_size);
    
    amps1 = amplitude(:,i);
    amps = amps1(base_trace(i)-imo_long_t < win_size & base_trace(i)-imo_long_t > - win_size);
    % for k = length(window)
    %     for j = 1:length(imo_long_t)
    % 
    %     if imo_long_t(j) == window(k)
    %         amps(k) = amplitude(j,i)
    %     else 
    %     end
    %     end
    
    % end
        pow = amps.^2;
        power(i) = max(pow);
        t_window = window(pow == max(pow));
        t_max(i) = t_window(1);
        pow = []; amps = []; amps1 = []; window = []; t_window = [];
       
end

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
plot(base_interp.Position_m_, base_interp.Depth_m_/v_int*v1, 'r-','LineWidth',4);

plot(imo_long_x(2:96),0.5*t_max*v1, 'm*','LineWidth',4);
% plot(debris_interp.Position_m_, debris_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
ylim([0 50]);
caxis([-5e6 5e6]);
axis ij;
colormap bone;

figure(3);
pcolor(imo_long_x, imo_long_elev, flipud(imo_long_shift)); 
shading interp; hold on;
plot(imo_long_x, topo_interp,'k', 'LineWidth', 2);
p0 = plot(topo_x, topo_z-base_interpolate,'r-', 'LineWidth', 4);
ylim([3730 max(imo_long_elev)]);
% xlim([0 100]);
caxis([-5e6 5e6]);
colormap bone;
% cmocean('balance');
set(gca, 'FontSize', 20);
xlabel('Horizontal Position (m)','FontSize',30);
ylabel('Elevation (m)', 'FontSize', 30);
% title('Lines 17 & 18','FontSize',30);

trace_a = 10;
trace_b = 86;

figure(4);
% subplot 141
% plot(10*log10(line14(:,trace_a)/max(max(line14))),0.5*v1*line14_t,'k','LineWidth',2); hold on;
% axis ij; 
% ylim([0 50]);
% 
% subplot 142
% plot(10*log10(line14(:,trace_b)/max(max(line14))),0.5*v1*line14_t,'k','LineWidth',2); hold on;
% axis ij; 
% ylim([0 50]);

subplot 121
plot((10*log10(line14(:,trace_a).^2)),0.5*v1*line14_t,'k','LineWidth',2); hold on;
plot(10*log10(power(trace_a)),0.5*t_max(trace_a)*v1, 'r*');
axis ij; 
ylim([0 50]);

subplot 122
plot((10*log10(line14(:,trace_b).^2)),0.5*v1*line14_t,'k','LineWidth',2); hold on;

plot(10*log10(power(trace_b)),0.5*t_max(trace_b)*v1, 'r*');axis ij; 
ylim([0 50]);

figure(5);
pcolor(line14_x, 0.5*v1*line14_t, 10*log10(line14.^2)); hold on;
plot([trace_a trace_a],[0 50],'r','LineWidth',2); hold on;
plot([trace_b trace_b],[0 50],'r','LineWidth',2);
plot(imo_long_x(2:96), 0.5*v1*t_max,'r*', 'LineWidth', 4);

shading flat; 
% plot(base_interp.Position_m_, base_interp.Depth_m_/v_int*v1, 'r-','LineWidth',4);
% plot(debris_interp.Position_m_, debris_interp.Depth_m_/0.1*v1, 'r-','LineWidth',2);
ylim([0 50]);
caxis([60 90]);
axis ij;
colormap parula;
colorbar;

%nadir oriented --> R^2 might be closer approximation 
geometric_correction2 = power.*base_depth.^2;
geometric_correction3 = power.*base_depth.^3;
geometric_correction4 = power.*base_depth.^4;

%normalize
power = power./max(power);
geometric_correction2 = geometric_correction2./max(geometric_correction2);
geometric_correction3 = geometric_correction3./max(geometric_correction3);
geometric_correction4 = geometric_correction4./max(geometric_correction4);



figure(6);
scatter(base_depth(3:78),10*log10(power(3:78)),'k');hold on;
scatter(base_depth(3:78),10*log10(geometric_correction2(3:78)),'r');
scatter(base_depth(3:78),10*log10(geometric_correction3(3:78)),'b');
scatter(base_depth(3:78),10*log10(geometric_correction4(3:78)),'c');


fit = polyfit(base_depth(3:78),10*log10(power(3:78)),1);
fit2 = polyfit(base_depth(3:78),10*log10(geometric_correction2(3:78)),1);
fit3 = polyfit(base_depth(3:78),10*log10(geometric_correction3(3:78)),1);
fit4 = polyfit(base_depth(3:78),10*log10(geometric_correction4(3:78)),1);

X = 0:40;
plot(X, fit(1)*X+fit(2),'k');
plot(X, fit2(1)*X+fit2(2),'r');
plot(X, fit3(1)*X+fit3(2),'b');
plot(X, fit4(1)*X+fit4(2),'c');

xlabel('Depth (m)');
ylabel('10log_{10}(power)');
legend('No correction','R^{2}','R^{3}','R^{4}','Location','southwest');
xlim([10 40]);
ylim([-30 0]);


% %account for geometric spreading loss 
% debris_interpolate = interp1(debris_interp.Position_m_, debris_interp.Depth_m_/v_int*v1, topo_x);
% debris_out = horzcat(topo_x,topo_E,topo_N,topo_z,debris_interpolate);
% writematrix(debris_out,strcat('topo/line16_debris.csv'));