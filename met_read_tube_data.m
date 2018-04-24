function traces = met_read_tube_data (file)

% function traces = met_read_tube_data (file)
%
% Read data from etracer CSV data files (vacuum tube anode voltage vs anode current as a function of grid voltage).
%
% INPUT:
% file: file path / name of etracer CSV data file
%
% OUTPUT:
% traces: array of structs with tube data (traces)
%
% EXAMPLE:
% t = met_read_tube_data ("path/to/300B.csv");
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
% Copyright (C) 2018 Matthias S. Brennwald.
% Contact: info@audioroot.net
% Further information: http://www.audioroot.net/MATAA

traces = [];
x = load (file);
k = unique (x(:,1)); % index to individual traces
for i = 1:length(k)

	% find i-th trace:
	j = find ( x(:,1) == k(i) );
	u.Ua = x(j(1),2:end); % anode voltage
	u.Ia = x(j(2),2:end); % anode current

	% find grid voltage:
	sweep_src = x(j(6),2:end); % sweep source (0: none / 1: NEGV1 / 2: HV2)
	Ug = repmat (NaN,size(u.Ua));
	for l = 1:length(Ug)
		if sweep_src(l) == 1
			Ug(l) =	x(j(5),l+1);
		elseif sweep_src(l) == 2
			Ug(l) =	x(j(3),l+1);
		end
	end
	u.Ug = mean(Ug(~isnan(Ug)));

	% remove NaN data:
	l = u.Ia .* u.Ug; l = find (~isnan(l));
	u.Ia = u.Ia(l); u.Ua = u.Ua(l);

	traces = [ traces , u ];
end
