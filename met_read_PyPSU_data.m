function traces = met_read_PyPSU_data (file)

% function traces = met_read_PyPSU_data (file)
%
% Read data from PyPSUcurvetrace data file(s).
%
% INPUT:
% file: file path / name of PyPSUcurvetrace data file; or cell string containing multiple file names to be combined into one data set.
%
% OUTPUT:
% traces: array of structs with tube data (traces).
%
% EXAMPLE:
% t = met_read_PyPSU_data ("path/to/6H30P.dat");
%
% DISCLAIMER:
% This file is part of MATETRACER.
%
% MATETRACER is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
% MATETRACER is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with MATETRACE. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2020 Matthias S. Brennwald.
% Contact: info@audioroot.net
% Further information: http://www.audioroot.net/MATAA

% check input:
if ~iscellstr(file)
	file = {file};
end

% load data and combine into single data set:
x = [];
for i = 1:length(file)
	x = [ x ; load(file{i}) ];
end

% remove values with current limiter on:
k = find (x(:,5) == 0); x = x(k,:);

% find V2 values:
V2 = unique(x(:,6));

traces = [];
for i = 1:length(V2)

	% find i-th trace:
	k = find (x(:,6) == V2(i));
	[u,j] = sort(abs(x(k,3)));
	v.Ua  = x(k,3)(j);		% anode voltage
	v.Ia  = x(k,4)(j)*1000;	% anode current
	v.Ug  = median(x(k,6));	% grid voltage

	traces = [ traces , v ];
end
