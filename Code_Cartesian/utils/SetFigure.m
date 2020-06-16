
function SetFigure(Hcf, Hca, FigPosition, FigFontSize, varargin)

% 13/07/2016 set some parameters of given figure;

% Hcf : the handle to the figure.
% Hca : the handle to the current asis.
% FigPosition: figure position to be located.
% FigFontSize: figure font size.


if ~isempty(varargin)
	SaveFig = varargin{1} ;
	FigName = datestr(now,'yyyymmddTHHMMSS') ;
	if length(varargin) >= 2
		FigName = varargin{2} ;
	end
end

% set(Hcf,'Position',[0 100 320 200]);
set(Hcf,'Position',FigPosition);
pause(0.2)

% set axis
set(Hca,... %			'Ylim',[0,0.14],...			'ytick',0:0.02:0.16,...
	'FontName','Times New Roman',...
	'FontSize',FigFontSize);
% set xlabel
h_XLabel=get(Hca,'XLabel');
set(h_XLabel,...
	'FontName','Times New Roman',...
	'FontSize',FigFontSize,...
	'Vertical','top');
% set ylabel
h_YLabel=get(Hca,'YLabel');
set(h_YLabel,...
	'FontName','Times New Roman',...
	'FontSize',FigFontSize,...
	'Vertical','baseline',... % While setting the 'Vertical' property, use 'baseline' | 'top' | 'cap' |'middle' | 'bottom' | 'baseline'..
	'Horizontal','center'); % While setting the 'Horizontal' property, use  'left' | 'center' | 'right'.
pause(0.2)

if SaveFig
% 	addpath('../export_fig')
	export_fig(Hcf, [FigName, '.fig'])  % export_fig(figure_handle, filename);
% 		export_fig(gcf, [FileName, '.pdf'])  % export_fig(figure_handle, filename);
	export_fig(Hcf, [FigName, '.eps'])  % export_fig(figure_handle, filename);
	pause(0.2)
end