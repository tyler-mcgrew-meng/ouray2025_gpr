
%%make a plot for mixing model results
%Supplementary Figure S5
clear all; clf; clc;

C = 3e-1; %speed of vacuum light [m/ns]

eps = [3:0.1:9]; %range of dielectric constants
rock_frac = zeros(length(eps),3); %array for rock fraction 

%unified mixing models from Sihvola, 2008 using mix_rule.m function
for i = 1:length(eps)
    for j = 1:3
        if j==1
        rock_frac(i,j) = mix_rule(3,9,eps(i),0);
        else
        rock_frac(i,j) = mix_rule(3,9,eps(i),j);
        end
        
    end
end

rg_vel = [0.14 0.15 0.168 0.165 0.16 0.147 0.17 0.109 0.128 0.112]; %wave speed results
% [gp19 sp19 gc20a gc19 gc16 sc20b sc20a egc23 sc23 iceland24]

rg_eps = C^2./rg_vel.^2;
rg_frac = 1-mix_rule(3,9,rg_eps,2); %fraction of ice according to each wavespeed measurement

crim = zeros(length(eps),1); 

%Two-phase CRIM model following Knight & Endres 1990
for i = 1:length(eps)
    crim(i) = (sqrt(eps(i))-sqrt(3))/(sqrt(9)-sqrt(3));
    crim(i) = 1-crim(i);
end

% 0.1378 0.143 0.1488

ice_frac = 1 - rock_frac;

figure(1); %Figure S7
plot(0.3./sqrt(eps),ice_frac(:,1),'LineWidth',4); hold on; grid on;
plot(0.3./sqrt(eps),ice_frac(:,2),'LineWidth',4);
plot(0.3./sqrt(eps),ice_frac(:,3),'LineWidth',4);
plot(0.3./sqrt(eps),crim,'LineWidth',4);
plot(rg_vel,rg_frac,'k+','MarkerSize',30,'LineWidth',5);
xlim([0.1 0.175]); ylim([0 1]);
legend('Maxwell Garnett','Brugemann','Coherent Potential','CRIM','CMP Result','Location','northwest','FontSize',30);
set(gca, 'fontsize', 25);
xlabel('v (m/ns)', 'fontsize', 50); ylabel('Ice Fraction','fontsize', 50); 
% title('Two-Phase Mixing Models');
