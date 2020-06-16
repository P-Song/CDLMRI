% generate radial sampling data and figure.

% add path to fessler's Michigan Image Reconstruction Toolbox (MIRT)
addpath('./utils');
addpath('../fessler/irt')
setup;

%% Generate radial sampling masks
NUMROWS = 256;
NUMCOLS = 256;
% Define undersampling ratio
% sampling_ratioA= 0.15;
sampling_ratioA= 0.5;
fold = round( 1/ sampling_ratioA);
% Generate sampling locations
numRadialsA=round(sampling_ratioA*NUMROWS); %number of required radial lines - calculated by the sampling ratio
numSamplesA = NUMROWS*2;
incr = @(ii) mod((ii-1)*180/1.618, 180);
k       =   gen_radial(incr, numSamplesA, numRadialsA,1,0,1);
k = reshape(k(:,1:numRadialsA,:),[],1,2);                               
norm_factor=max(abs(k(:)));
k=k*pi/norm_factor; %k is a matrix that holds the sampling locations in the k-space domain


%% show figures
x = squeeze(k(:,1,1));
y = squeeze(k(:,1,2));
figure; plot(x, y,'w.' , 'markersize', 2)
set(gcf, 'Color', 'k'); % set figure background to be black
set(gcf,'Position', [100,100,256, 256]); axis off
set(gca, 'position', [-0.5, -0.5, 2, 2]) ; % set(gca, 'Position', [0,0,1,1]);

%%
addpath('../export_fig');
FigFormat = ['.eps']; % ['.png'] ; % ['.jpg'] ; % ['.eps'];
FileName = ['RadialMask7Fold'];
% export_fig(gcf,  [FileName, FigFormat], '-nocrop', '-transparent')
export_fig(gcf,  [FileName, FigFormat], '-nocrop')



