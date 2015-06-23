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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
isOctave = exist('OCTAVE_VERSION') ~= 0;
if (isOctave)
   page_output_immediately(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% running all 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for test_dir={'graybat','mpi','threaded'}
    cd(test_dir{1})
    pwd
    tic
    %laserPumpCladdingExample
    toc 
    %extract_gain_map
    cd ..
end

addpath(pwd)

for test_dir={'graybat','mpi','threaded'}
    cd(test_dir{1})
    gain_analysis
    cd ..
end

exit
