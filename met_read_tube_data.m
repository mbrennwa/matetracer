function traces = met_read_tube_data (file,unit)

% function traces = met_read_tube_data (file,unit)
%
% Read data from data files with vacuum tube anode voltage vs anode current as a function of grid voltage (eTracer CSV or uTracer UTD format). The file format is determined from the file extension.
%
% INPUT:
% file: file path / name of data file (eTracer *.CSV or uTracer *.UTD format)
% unit (optional): index to tube unit (for multi-unit tubes, such as double triodes). Default: unit = 1.
%
% OUTPUT:
% traces: array of structs with tube data (traces). If length(unit)>1, traces is a cell array, where each cell corresponds to one unit.
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

if ~exist ('unit','var')
	unit = 1;
end

if length (unit) > 1
	for k = 1:length(unit)
		traces{k} = met_read_tube_data (file,unit(k));
	end
else

	traces = [];

	[DIR, NAME, EXT] = fileparts(file);

	switch toupper(EXT)

		case '.CSV' % eTracer CSV file
			x = load (file);
			k = unique (x(:,1)); % index to individual traces
			for i = 1:length(k)

				% find i-th trace:
				j = find ( x(:,1) == k(i) );
				if unit == 1 % tube unit 1
					ku = 1; ki = 2;
				elseif unit == 2 % tube unit 2
					ku = 3; ki = 4;
				else
	unit
					error ('ETD files do not contain data for more than 2 tube units.')
				end			
				u.Ua = x(j(ku),2:end); % anode voltage
				u.Ia = x(j(ki),2:end); % anode current

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

		case '.UTD' % uTracer UTD file
			x = dlmread(file);
			k = unique (x(:,2)); % index to individual traces

			if unit > 1
				error ('UTD data files with more than one tube unit not supported!')
			end

			for i = 1:length(k)

				% find i-th trace:
				j = find ( x(:,2) == k(i) );

				u.Ua = x(j,6); % anode voltage
				u.Ia = x(j,3); % anode current
				u.Ug = mean(x(j,5)); % grid voltage

				% remove NaN data:
				l = u.Ia .* u.Ug; l = find (~isnan(l));
				u.Ia = u.Ia(l); u.Ua = u.Ua(l);

				traces = [ traces , u ];
			end


		otherwise
			error (sprintf("File %s: unknown file format.",file));

	end
end
