%Imogene and Gilpin Rock Glaciers
clear all; clf; clc;
C = 0.3; %these velocities are calculated in meters/nanosecond

imo = readtable("cmp3-Point.csv",NumHeaderLines=3);

% base = readtable("line14_base-Point.csv",NumHeaderLines=3);


[v1, t01,sig1, cov1,d_hat1] = cmpLS(imo.Position_m_, imo.Time_ns_);

VLS = [v1]; 
T0LS = [t01]; 

t1 = (t01.^2+4.*(imo.Position_m_./2).^2./v1^2).^(1/2);

t0a = min(imo.Time_ns_); %t0 pick (ns)

%NMO velocities
va = sqrt(imo.Position_m_.^2./(imo.Time_ns_.^2-t0a^2));

% t0s = min(L14_at); %air/surface wave

X = 0:40;
%test airwave speed
% [Sm, Sb] = polyfit(L14_ax, L14_at, 1); Surf = Sm(1)*X+Sm(2);
% Vs = 1/Sm(1);
% disp(Vs); 

T0 = [t0a];
V = [mean(va(length(va)-10:length(va)))];
%Real(dielectric permitivity)
ea = C^2./va.^2; % base 50 

sqrt1 = sqrt(sig1)*ones(length(t1),1);

%propagate travel time standard deviations to velocity uncertainty
%za, zb, zc are second order terms that are ignored
[v1p, t01p ,za, zb, zc] = cmpLS(imo.Position_m_, t1+sqrt1); %v and t0 due to positve time uncertainty
[v1n, t01n ,za, zb, zc] = cmpLS(imo.Position_m_, t1-sqrt1); %v and t0 due to negative time uncertainty 


%array of velocity error
%time uncertainty inversely related to velocity uncertainty
neg = [v1-v1p];  
pos = [v1n-v1];

d1 = 0.5*v1*t01;
d1n = 0.5*v1n*t01;
d1p = 0.5*v1p*t01;


figure(1); %Figure S6 B & C
subplot 121
axis ij; hold on;
% plot(L14_ax, L14_at,'o','Color','#EDB120','LineWidth',2);
%plot(X, Surf, X, Ref,'LineWidth',2); 
p1 = plot(imo.Position_m_,imo.Time_ns_,'o','Color','#4DBEEE','LineWidth',2);
p2 = plot(imo.Position_m_, t1, 'k-','LineWidth',2);

set(gca, 'FontSize', 20);
xlabel('x (m)', 'FontSize', 30); ylabel('t (ns)', 'FontSize', 30); ylim([0 800]);
legend([p1 p2], 'pick', 'best fit','Location','southwest');

subplot 122
p3 = plot(V, T0,'-o','LineWidth',2,'Color','#4DBEEE'); axis ij; hold on;
p4 = errorbar(VLS, T0LS,neg,pos,'horizontal','k-o','LineWidth',2);
set(gca, 'FontSize', 20);
xlabel('v (m/ns)', 'FontSize', 30); ylabel('t_0 (ns)', 'FontSize', 30);
ylim([0 800]); xlim([0.1 0.20]);
ticks = 0.1:0.01:0.20; xticks(ticks); 
xticklabels({'0.1','','','','','0.15','','','','','0.20'});
legend([p3 p4], 'high offset', 'best fit','Location','southwest');

% %compare dip angle equations from common offset and cmp
delta_t = 100; delta_x = 30;
vrange = [0.08:0.00001:0.24];
theta_cmp = acosd(vrange./v1);
theta_coh = asind(vrange*delta_t/(2*delta_x));

%arrays of uncertainties for wave speed/dip plot
NEG = ones(length(vrange),1)*neg(1);
POS = ones(length(vrange),1)*pos(1);

vice = [0 90;0.17 0.17]; %ice velocity for plotting (m/ns)

figure(2); %Figure S6 D
e = errorbar(vrange(real(theta_cmp)>0), real(theta_cmp(real(theta_cmp)>0)),...
    NEG(real(theta_cmp)>0),POS(real(theta_cmp)>0),'horizontal','y-s'); 
hold on;
plot(vrange, theta_coh, 'r-','LineWidth',5);
plot(vice(2,:), vice(1, :), 'b--', 'LineWidth',2);
set(gca, 'FontSize', 20);
xlabel('v (m/ns)', 'FontSize', 30); ylabel('Dip (degrees)', 'FontSize', 30);
legend('common midpoint','common offset');
xlim([0.12 0.2]); ylim([0 60]);
ticks = 0.12:0.01:0.2; xticks(ticks); 
xticklabels({'0.12','','0.14','','0.16','','0.18','','0.2'});
e.MarkerSize = 5; e.MarkerEdgeColor = 'k';
e.MarkerFaceColor = 'k';

format compact;
disp('Best Fit Velocity: '); disp(v1);
disp('Standard Deviation: '); disp((v1n-v1p)/2);
disp('Best Fit Dielectric Constant:'); disp((C/v1)^2);


%convert to dielectric constant
disp('Max. Dielectric Constant:'); disp((C/v1p)^2);
disp('Min. Dielectric Constant:'); disp((C/v1n)^2);
disp('Best Fit T0:'); disp(t01);
disp('Best Fit Depth:'); disp(d1);