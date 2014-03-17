% This will calculate the gain distribution within an Yb:YAG-Crystal
% Monochromatic approach                Daniel Albach
%                   last edit:  2008/10/24

clear all;
tic;
%constants
const.N1per = 1.38e20;
const.c = 3e8;
const.h = 6.626e-34;

% crystal
crystal.doping = 2;
crystal.length = 0.7; % [cm]
crystal.tfluo = 9.5e-4;%1.0e-3 pour DA
crystal.nlexp = 1;
% nonlinearity of the doping - exponential - factor front/end of the
% crystal

% steps
steps.time = 1000;
steps.crys = 300;

% In later versions this should be replaced by a vector for the emission
% and absorption cross sections as well for the input spectrum

% pump
pump.s_abs = 0.76e-20; % absorption cross section in cm² (0.778e-20 pour DA)
pump.s_ems = 0.220e-20; % emission cross section in cm²(0.195e-20 pour DA)
pump.I = 10e3; % Pump intensity in W/cm²
pump.T = 1e-3; % pump duration in s
pump.wavelength = 940e-9;% pump wavelength in m

% laser
laser.s_abs = 1.10e-21; %cm2(1.16e-21 pour DA)
laser.s_ems = 2.40e-20; %cm2(2.48e-20)
laser.I = 1e6; % laser intensity
laser.T = 1e-8; % laser duration
laser.wavelength = 1030e-9; % lasing wavelength in m

% modus definitions
mode.BRM = 1; % 1 Backreflection, 0 Transmissionmode
mode.R = 1; % reflectivity of the coating

% constants for short use
c = const.c; %m/s
h = const.h; %Js

%  prepare the beta_crystal as a vector to be a parameter, global scope!
beta_crystal = zeros(steps.crys,1);
pulse = zeros(steps.time,1);


time_step = pump.T/(steps.time-1);
crystal_step = crystal.length/(steps.crys-1); 

gain_local = zeros(steps.crys,1);

% this is the function call to the pump
[beta_crystal,beta_store,pulse,Ntot_gradient] = beta_int(beta_crystal,pulse,const,crystal,steps,pump,mode);

% grids for plotting the lines
grid_t = 0:time_step:pump.T;
grid_z = 0:crystal_step:crystal.length;

% laser wavelength
sigma_abs_L=1.16e-21;
sigma_ems_L=2.48e-20;
beta_min = zeros(25,1);
beta_min(:) = sigma_abs_L/(sigma_abs_L+sigma_ems_L);
grid_z_beta_min = 0:crystal.length/24:crystal.length;

% integration yields the energy density, this has to be multiplied with the
% "average" surface of the pumped area
energy_density.before = trapz(grid_z,beta_crystal.*Ntot_gradient)*(h*c/pump.wavelength);
energy_density_ex.before = trapz(grid_z,(beta_crystal-beta_min(1)).*Ntot_gradient)*(h*c/pump.wavelength);

% overwrite the pulse with zeros
pulse = zeros(steps.time,1);

% this is the function call to the extraction
[beta_crystal_l,beta_store_l,pulse] = beta_int(beta_crystal,pulse,const,crystal,steps,laser,mode);

% integration yields the energy density, this has to be multiplied with the
% "average" surface of the pumped area
energy_density.after = trapz(grid_z,beta_crystal_l.*Ntot_gradient)*(h*c/pump.wavelength);

% get the optimum length at the end of the pump time
beta_fit_comp = polyfit(transpose(grid_z),beta_crystal,4);
beta_fit =polyval(beta_fit_comp,grid_z);

intersect=[beta_fit_comp(1) beta_fit_comp(2) beta_fit_comp(3) beta_fit_comp(4) beta_fit_comp(5)-beta_min(1)];
roots_beta = roots(intersect);

intersection = 0;

for iroots=1:length(roots_beta)
    if isreal(roots_beta(iroots)) && (roots_beta(iroots) > 0)
        intersection = roots_beta(iroots);
    end
end

time = toc
% now calculate the local gain at each point
gain_local(:) = -(sigma_abs_L - beta_crystal(:)*(sigma_abs_L+sigma_ems_L)).*Ntot_gradient(:);

% first image
figure(1);
fig_1=plot(grid_z,beta_crystal,'-',grid_z_beta_min,beta_min,'--');figure(gcf);
title('\beta distribution in absense of ASE','FontSize',16);
xlabel('crystal thickness [cm]','FontSize',16);
ylabel('beta','FontSize',16);
xlim([0 crystal.length]);
text(0.02,0.025,'\uparrow \beta_{min} for 1030nm','FontSize',14)
text(0.2,0.3,[num2str(crystal.doping),'%doping'],'FontSize',14);
set(fig_1,'LineWidth',2);
set(gca,'FontSize',13,'FontName','Calibri');
set(gca,...
    'TickDir', 'out',...
    'TickLength', [0.02 0.02],...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth', 1);

%now the beta after the run with the laser
figure(4);
plot(grid_z,beta_crystal_l,'-',grid_z_beta_min,beta_min,'--');figure(gcf);
title('\beta distribution after one pass');
xlabel('crystal length [cm]');
ylabel('beta');
xlim([0 crystal.length]);
text(0.02,0.025,'\uparrow \beta_{min} for 1030nm')
text(0.2,0.2,[num2str(crystal.doping),'%doping']);

figure(2);
[haxes,hline1,hline2] = plotyy(grid_z,beta_crystal,grid_z,gain_local);figure(gcf);
title('\beta distribution in absense of ASE','FontSize',16);
xlabel('crystal thickness [cm]','FontSize',16);
axes(haxes(1));
ylabel('\beta (z)','FontSize',14);
xlim([0 crystal.length]);
axes(haxes(2));
ylabel('gain (z)','FontSize',14);
xlim([0 crystal.length]);
set(hline1,'LineWidth',2);
set(hline2,'LineWidth',2);
set(haxes,'FontSize',13,'FontName','Calibri');
set(haxes,...
    'TickDir', 'out',...
    'TickLength', [0.02 0.02],...
    'XMinorTick', 'on',...
    'YMinorTick', 'on',...
    'LineWidth', 1);
% 
% 
% % grid for second plot
% grid_chart_x = 0:crystal_step:crystal.length;
% grid_chart_y = pump.T*1000:-time_step*1000:0;
% 
% % second image
% figure(3);
% beta_store = rot90(beta_store);
% imagesc(grid_chart_x,grid_chart_y,beta_store);figure(gcf);
% title('Evolution of \beta in absense of ASE as afunction of time and position within the crystal');
% xlabel('position in the crystal [cm]');
% ylabel('time [ms]');
% axis xy;