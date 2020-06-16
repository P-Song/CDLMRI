% 20/09/2017 Generate Cartesian Masks.

clear;
addpath('./VariousMasks/undersampling');
				
%% generate a new Cartesian random sampling masks along x direction

% load and adjust true images
NUMROWS = 256;
NUMCOLS = 256;
center_block = [0 0]; % size of full sampling area in the center
cart_random = [6, 36]; % [0, 32]; % [0 8]; % # of fixed line in the center and # of random lines	
s_image = [NUMROWS, NUMCOLS] ;
coverage = s_image; % size of the k-space in which we want to sample.
n_fixed_lines = cart_random(1,1);  % number of fixed, centered lines in the phase encode direction
n_random_lines = cart_random(1,2); % number of random lines in the phase encode direction
cartesianX_random = MRI2_sample_cartesian_x_random_quadratic(coverage, [], n_fixed_lines, n_random_lines);
cartesianX_random = [cartesianX_random MRI2_sample_cartesian_x_standard(center_block)];
QA = MRI_sample_backward(ones(size(cartesianX_random,2),1), cartesianX_random, s_image);

indexQA=find(QA==1); %Index the sampled locations in sampling mask
fold = round( (NUMROWS * NUMCOLS) / length(indexQA) );
figure; imagesc(QA); colormap gray
set(gcf, 'position', [100, 100, NUMROWS, NUMCOLS]);
set(gca, 'position', [0,0,1,1]);
imwrite(255.*QA,colormap(gray(256)),sprintf('Q1.png', []));


%% load masks
QA = imread('./Cart42/Q1.png');
QA = (QA>0.1);
% QA = ifftshift( QA );
indexQA=find(QA==1); %Index the sampled locations in sampling mask
fold = round( (NUMROWS * NUMCOLS) / length(indexQA) );

disp('done!')