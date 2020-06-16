

function [sampling_mask, varargout] = GenMask(params)

saveMask = 0; % 1: save mask; 0: do not save mask;

sampling_factor = 1/params.fold;
% sampling_factor = 1/40;    % m = oversample_factor*n;
%                                          % m: number of measurements; n: image dimensions
										 
% radius_fully_sampled = 0.6;
radius_fully_sampled = 0;

N1 = params.N1;
N2 = params.N2;


% =========================================================================
% Create FFT and sampling operators: use genPDF from Lustig's toolbox and
% also sparco

% %  e.g. [pdf,val] = genPDF(imSize,p,pctg,distType,radius,disp)
% for fold <= 20, set p = 7
% pdf_quartic_poly  = genPDF([N1, N2], 7, sampling_factor, 2, ...
% radius_fully_sampled*sampling_factor, 0);

% for fold = 30, set p = 10
p = 20;
pdf_quartic_poly  = genPDF([N1, N2], p, sampling_factor, 2, ...
radius_fully_sampled*sampling_factor, 0);

random_matrix_aux = rand(N1, N2);

sampling_mask = (random_matrix_aux < pdf_quartic_poly);

% Number of measurements
m = sum(sum(sampling_mask));

if saveMask 
	FileName = ['Q', num2str(round(1/sampling_factor)), 'fold'] ;
	save(FileName, 'sampling_mask') ;
end

fold = N1*N2/m;

% % Compose restriction operator with FFT2d
% Phi = opFoG(opRestriction(n, find(sampling_mask)), opFFT2d(N1, N2));
% =========================================================================

