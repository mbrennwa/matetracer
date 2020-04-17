function fig = met_plot_curves (traces,Pmax,annot,fig_save)

% function fig = met_plot_curves (traces,Pmax,annot,fig_save)
%
% Plot data from etracer CSV data files (vacuum tube anode voltage vs anode current as a function of grid voltage).
%
% INPUT:
% traces: tube tracer data (struct array, see met_read_tube_data.m). Only single tube units are supported (no cell array with multiple units).
% Pmax (optional): max anode dissipation of tube (Watt). Default: Pmax = Inf.
% annot (optional): annotation text (top right on figure)
% fig_save (optional): name of file for saving the graphics. The file format is determined from the file suffix (*.eps, *.pdf, or *png)
%
% OUTPUT:
% fig: figure handle
%
% EXAMPLE:
% met_plot_curves ("path/to/300B.csv",40,"300B","300b.pdf");
% met_plot_curves ("path/to/12AT7.utd",2.5);
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

% check input:
if ~exist('Pmax','var')
	Pmax = Inf;
end
if ~exist('fig_save','var')
	fig_save = '';
end
if ~exist('annot','var')
	annot = '';
end

lstyle = '-';
lcolor = [1 1 1];
lwidth = 1.3;
fsize  = 10;
fname  = "Helvetica";

xmax = ymax = 0;


% get figure handle:
fig = gcf();
set (fig, "defaultaxesfontname", fname)
set (fig, "defaultaxesfontsize", fsize)
set (fig, "defaulttextfontname", fname)
set (fig, "defaulttextfontsize", fsize)

Ug_labels = repmat(NaN,length(traces),2); % text handle / trace slope pair (needed for moving labels later)
for i = 1:length(traces)

	if length(traces(i).Ua) > 1

		% make sure resolution gives smooth curve (use spline interpolation)
		if length(traces(i).Ua) < 100
			xx = linspace ( min(traces(i).Ua),max(traces(i).Ua),100 );
			traces(i).Ia = interp1 ( traces(i).Ua,traces(i).Ia,xx,'spline' );
			traces(i).Ua = xx;
		end

		% throw away data that is far beyond Pmax:
		l = find (traces(i).Ia/1000.*traces(i).Ua > 1.3 * Pmax);
		traces(i).Ia(l) = NaN; traces(i).Ua(l) = NaN;
		
		% remove NaN data:
		traces(i).Ua = traces(i).Ua(~isnan(traces(i).Ua));
		traces(i).Ia = traces(i).Ia(~isnan(traces(i).Ia));

		% add Zero:
		traces(i).Ua = [0,traces(i).Ua];
		traces(i).Ia = [0,traces(i).Ia];

		% plot trace:
		plot ( traces(i).Ua,traces(i).Ia , 'linestyle',lstyle , 'color',0*lcolor , 'linewidth',lwidth );
		hold on;

		% add grid voltage label (move label position later, when plot scaling is done):
		Ug_labels(i,1) = text (traces(i).Ua(end),traces(i).Ia(end),sprintf('%i V',traces(i).Ug) , 'verticalalignment','middle', 'horizontalalignment','center' , 'fontsize',0.8*fsize);

		s = diff(traces(i).Ia)./diff(traces(i).Ua);
		Ug_labels(i,2) = s(end);

		% check max. x/y extent of trace:
		xmax = max([max(traces(i).Ua),xmax]);
		ymax = max([max(traces(i).Ia),ymax]);

	end
end

% move Ug labels:
for i = 1:length(traces)
	if length(traces(i).Ua) > 1
		p = get (Ug_labels(i,1),'position');
		ext = get (Ug_labels(i,1),'extent');
		%rectangle('position',ext);
		s0 = ext(4)/ext(3);
		if Ug_labels(i,2) > s0 % move label above the end of the line, an a bit to the right
			sy = ext(4) / 2;
			sx = sy / Ug_labels(i,2);
		else % move label to the right of the end of the line, and a bit up
			sx = ext(3) / 2;
			sy = sx * Ug_labels(i,2);
		end
		p = [ p(1)+1.5*sx , p(2)+1.5*sy , p(3) ];
		set (Ug_labels(i,1),'position',p)

		% check max. x/y extent of labels:
		xmax = max([p(1)+ext(3)/2,xmax]);
		ymax = max([p(2)+ext(4)/2,ymax]);

	end
end

% determine grid spacing:
g  = [0.1 0.25 0.5 1 2.5 5 10 25 50 100 250 500 1000 2500 5000 10000];

[u,kx] = min (abs(xmax/8./g-1)); dx = g(kx);
[u,ky] = min (abs(ymax/8./g-1)); dy = g(ky);

% determine axes limits:
xmax = dx * ceil (xmax/dx);
ymax = dy * ceil (ymax/dy);

% plot Pmax curve:
x = linspace (0,xmax,1001); x = x(1:end);
y = 1000*Pmax ./ x;
% pa=patch ( [x x(end) x(1)],[y y(1) y(1)],'k','facealpha',0.1)
plot (x,y , 'linestyle','--' , 'color',0*lcolor , 'linewidth',lwidth )

% annotation text:
if ~isempty(annot)
	h = text ( 0.95*xmax,0.95*ymax , annot , 'verticalalignment','top', 'horizontalalignment','right' , 'fontsize',fsize );
end

% done with all the elements in the plot:
hold off

% add axes labels:
xlabel ( 'ANODE VOLTAGE (V)' , 'fontsize',fsize )
ylabel ( 'ANODE CURRENT (mA)' , 'fontsize',fsize )

% format plot area (axes limits, grid):
axis ([0 xmax 0 ymax]);
dn = 5;
set ( gca , 'xtick',[0:dx/dn:xmax] , 'ytick',[0:dy/dn:ymax] ); % grid line positions
% set ( gca , 'tickdir','out' ); % ticks out
set ( gca , 'ticklength', [0,0] ); % no tickmarks
set ( gca , 'linewidth',lwidth/4 , 'fontsize',fsize); % line width, font size
set ( gca , 'gridcolor',0.85*lcolor ); % grid line color

xt = get ( gca , 'xticklabel' );
for k = 1:length(xt)
	if rem(k-1,dn)
		xt{k} = "";
	end
end
set ( gca , 'xticklabel',xt );

yt = get ( gca , 'yticklabel' );
for k = 1:length(yt)
	if rem(k-1,dn)
		yt{k} = "";
	end
end
set ( gca , 'yticklabel',yt );

set ( gca , 'gridcolor',0.85*lcolor ); % grid line color
grid on;

% make a fat frame:
h = rectangle ( 'Position',[0 0 xmax ymax] , 'EdgeColor','k' ,'linewidth',lwidth );

% print to file:
if ~isempty(fig_save)
	H = 4; W = 6;
	set(fig,'PaperUnits','inches')
	set(fig,'PaperOrientation','portrait');
	set(fig,'PaperSize',[W,H])
	set(fig,'PaperPosition',[0,0,W,H])
	
	drawnow; pause(0.1) 

	[DIR, NAME, EXT] = fileparts (fig_save);
	switch toupper(EXT)
		case {'.PDF' '.PNG'}
			print(fig,fig_save);
		case '.EPS'
			print(fig,fig_save,'-depsc2');
		otherwise
			warning ('Unknown graphics type. Figure not saved to file...')
	end
end
