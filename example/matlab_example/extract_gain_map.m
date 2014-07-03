 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
 %
 % This file is part of HASENonGPU
 %
 % HASENonGPU is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, either version 3 of the License, or
 % (at your option) any later version.
 %
 % HASENonGPU is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 % GNU General Public License for more details.
 %
 % You should have received a copy of the GNU General Public License
 % along with HASENonGPU.
 % If not, see <http://www.gnu.org/licenses/>.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% extraction of gain map out of the calculation
% just make it in the center - position

% at first try with one distribution

clear all; 
tic
timeslices_tot = 150;


for i_s=1:timeslices_tot
    disp(['TimeSlice ' num2str(i_s)]);

    isOctave = exist('OCTAVE_VERSION') ~= 0;
    if (isOctave)
      fflush(stdout);
    else
      drawnow('update');
    end

    i_ss = i_s-1;
    
    filename=['save_' num2str(i_ss) '.mat'];
    
    load(filename);
    
    time_step = pump.T/(steps.time-1);
    crystal_step = crystal.length/(steps.crys-1);

    gainmap = zeros(size(p,1),1);

    for i_p=1:size(p,1)

        beta_crystal(1,:) = beta_cell(i_p,:);

        laser.I = 1;

        % grids for plotting the lines
        grid_t = 0:time_step:pump.T;
        grid_z = 0:crystal_step:crystal.length;

        time_step = pump.T/(steps.time-1);
        crystal_step = crystal.length/(steps.crys-1); 

        gain_local = zeros(steps.crys,1);

        pulse = zeros(steps.time,1);
        pulse(:) = laser.I;

        % extraction mode
        mode.BRM = 1; % 1 Backreflection, 0 Transmissionmode
        mode.R = 1; % reflectivity of the coating
        mode.extr =1;

        % laser wavelength
	[laser.max_ems, i] = max(laser.s_ems);
        laser.max_abs = laser.s_abs(i);

        sigma_abs_L=laser.max_abs;
        sigma_ems_L=laser.max_ems;
        beta_min = zeros(25,1);
        beta_min(:) = sigma_abs_L/(sigma_abs_L+sigma_ems_L);
        grid_z_beta_min = 0:crystal.length/24:crystal.length;

        % time constants for the pump
        time_step_ex = laser.T/(steps.time-1);
        grid_t_ex = 0:time_step_ex:laser.T;

	if (isOctave)
	  energy_pulse = trapz(grid_t_ex, pulse');
	else
          energy_pulse = trapz(grid_t_ex, pulse);
	end

        % energy dump for the plot
        energy_pulse_round(1,1)=energy_pulse;

        % integration yields the energy density, this has to be multiplied with the
        % "average" surface of the pumped area
        % energy_density.before = trapz(grid_z,beta_crystal.*Ntot_gradient)*(h*c/pump.wavelength);
        % energy_density_ex.before = trapz(grid_z,(beta_crystal-beta_min(1)).*Ntot_gradient)*(h*c/pump.wavelength);
    %     beta_crystal_l = beta_crystal;
    %     beta_store_l=beta_store;

    %     for iroundtrips=1:1
        % this is the function call to the extraction
        % we have to call it two times seperately for the aller-retour!
        [beta_crystal,beta_store,pulse] = beta_int(beta_crystal,pulse,const,crystal,steps,laser,mode);
%         gain_local(:) = -(sigma_abs_L - beta_crystal(:)*(sigma_abs_L+sigma_ems_L)).*Ntot_gradient(:);

        gainmap(i_p,1) = pulse(1)/laser.I;

        % now the pulse is at the backside

        % energy_pulse = trapz(grid_t_ex,pulse);
        % energy_pulse_round(iroundtrips+1,1)=energy_pulse;
    %     pulse_round(:,iroundtrips)=pulse(:);
        % after each pass he sees some losses e.g. 5%
        % now make the losses, but be carefull - the position changes the energy
        % output!!!
        % pulse(:)=pulse(:).*mode.R_ex;
    %     end
        % integration yields the energy density, this has to be multiplied with the
        % "average" surface of the pumped area
        % energy_density.after = trapz(grid_z,beta_crystal_l.*Ntot_gradient)*(h*c/pump.wavelength);

        % bar (energy_pulse_round); figure(gcf)
    end
    
%     figure(2)
    [x_grid,y_grid]=meshgrid(-3:0.5:3);
    gainmap_Interp(:,:,i_s) = griddata(p(:,1),p(:,2),gainmap(:,1),x_grid,y_grid);
    
%     imagesc(gainmap_Interp);
%     axis equal;
%     colorbar;
end
% now make an interpolation to make it plotable

%% extract data
for i_t=1:timeslices_tot
    gain_line(i_t,1)=gainmap_Interp(7,7,i_t);
end
x=fopen('gain_line.txt','w');
fprintf(x,'%.50f\n',gain_line);
fclose(x);
toc
