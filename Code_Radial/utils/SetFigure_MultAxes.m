
function SetFigure_MultAxes(Hcf, Hca, varargin)
% 10/01/2017 add figure format opertion
% 12/08/2016 set some parameters of given figure with multiple axes;

% Hcf : the handle to the figure.
% Hca : the handle to the current asis.
% varargin{1}: FigPosition: figure position and size, e.g., [0 100 320 200]);
% varargin{2}: FigFontName: e.g. 'Times New Roman';
% varargin{3}: FigFontSize: figure font size, e.g. 10.
% varargin{4}: SaveFig
% varargin{5}: FigName
% varargin{6}: FigFormatCell, e.g., .eps, .fig, .pdf

% example 
% SetFigure_MultAxes(Hcf, Hca, ...
% 	[0 100 320 200],...
% 	'Times New Roman', ...
% 	10, ...
% 	0, ...
% 	'xxxx', ...
% 	{ '.eps', '.fig'} ) ;


FigPosition = 'default';
FigFontName = 'default';
FigFontSize = 'default';
SaveFig = 0;
FigName = datestr(now,'yyyymmddTHHMMSS') ;
FigFormatCell = { '.eps' };

if ~isempty(varargin)
	for i = 1 : numel(varargin)
		switch i 
			case 1
				FigPosition = varargin{1} ;
			case 2
				FigFontName = varargin{2} ;
			case 3
				FigFontSize = varargin{3} ;
			case 4
				SaveFig = varargin{4} ;
			case 5
				FigName = varargin{5} ;
			case 6
				FigFormatCell = varargin{6} ;
			otherwise,
				disp('undefined options !')
		end
	end
end

set(Hcf,'Position',FigPosition);
pause(0.2)


% Hcf_Ch(1): the last legend; Hcf_Ch(2): the last axis; ... 
% Hcf_Ch(end-1): the first legend; Hcf_Ch(end): the first axis;
Hcf_Ch = get(Hcf, 'Children') ; 

for i = 1 : length(Hcf_Ch)
	set(Hcf_Ch(i),... 
		'FontName',FigFontName,...
		'FontSize',FigFontSize); % this command will change the font for both axis and legend.
end

if SaveFig
% 	addpath('../export_fig')
	for i = 1 : numel(FigFormatCell)
		export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
	end
	pause(0.2)
end





