% 29/11/2017 save the figure and data. Add p in the name.

clear;
saveMask = 1; % 1: save mask; 0: do not save mask;
% sampling_factor = 1/6;    % m = oversample_factor*n;
                                         % m: number of measurements; n: image dimensions

% sampling_factor = 1/5; 
% sampling_factor = 1/10; 
% sampling_factor = 1/20; 
sampling_factor = 1/40; 
		
% radius_fully_sampled = 0.6;
radius_fully_sampled = 0;

% N1 = 128;
% N2 = 128;

N1 = 256;
N2 = 256;
% =========================================================================
% Create FFT and sampling operators: use genPDF from Lustig's toolbox and
% also sparco

% %  e.g. [pdf,val] = genPDF(imSize,p,pctg,distType,radius,disp)
% for fold <= 20, set p = 7
% pdf_quartic_poly  = genPDF([N1, N2], 7, sampling_factor, 2, ...
% radius_fully_sampled*sampling_factor, 0);

% for fold = 30, set p = 10
p = 10; % 20;
pdf_quartic_poly  = genPDF([N1, N2], p, sampling_factor, 2, ...
radius_fully_sampled*sampling_factor, 0);

random_matrix_aux = rand(N1, N2);

sampling_mask = (random_matrix_aux < pdf_quartic_poly);

figure; imagesc(sampling_mask); colormap gray;
set(gcf, 'position', [100,100, 256, 256]);
set(gca, 'position', [0, 0, 1, 1]) ;
set(gcf, 'Color', 'w');

% Number of measurements
m = sum(sum(sampling_mask));
fold = N1*N2/m

if saveMask 
	FileName = ['Q_', num2str(round(1/sampling_factor)), 'fold', '_p', num2str(round(p))] ;
	save(FileName, 'sampling_mask') ;
	
	set(gcf, 'Color', 'w');
	export_fig(gcf,  [FileName, '.eps'], '-nocrop', '-transparent')
end


% Compose restriction operator with FFT2d
Phi = opFoG(opRestriction(n, find(sampling_mask)), opFFT2d(N1, N2));
% =========================================================================

