 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Copyright 2015 Erik Zenker, Carlchristian Eckert
 %
 % This file is part of HASEonGPU
 %
 % HASEonGPU is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, either version 3 of the License, or
 % (at your option) any later version.
 %
 % HASEonGPU is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 % GNU General Public License for more details.
 %
 % You should have received a copy of the GNU General Public License
 % along with HASEonGPU.
 % If not, see <http://www.gnu.org/licenses/>.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % This script should be started INSIDE a folder that contains the gain_line.txt file.
 % The output will be placed in files in the parent folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isOctave = exist('OCTAVE_VERSION') ~= 0;
if (isOctave)
   page_output_immediately(1);
end

tic

pwd
[pathstr,current_dir,ext] = fileparts(pwd);

data = load('gain_line.txt');
cd ..
baseline = load('gain_line_baseline.txt');

str = strcat(current_dir, '   max difference to baseline:');
fid = fopen('results.txt', 'at+');
fprintf(fid, '%s\n', str);
fprintf(fid, '%.10f\n',max(abs(baseline - data)));
fclose(fid);

%name = genvarname(strcat('gain_line_',current_dir));
%eval([ name ' = data ']);

x = (0:size(data)-1)';
curve = [ (20*x)' ; (data.*data * 1.0263)' ]';
curvename = strcat('gain_curve_', current_dir, '.txt');
fid = fopen(curvename , 'w');
fprintf(fid, '%s\n', '# time[us]        gain');
fclose(fid);
dlmwrite(curvename, curve, '-append', 'delimiter', ' '); 

cd(current_dir)
toc
return
