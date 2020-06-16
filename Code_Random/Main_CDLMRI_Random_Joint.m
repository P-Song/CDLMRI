% =========================================================================
% Coupled Dictionary Learning for Multi-contrast MRI Reconstruction
% =========================================================================
%
%This software is to perform Coupled Dictionary Learning based Multi-contrast MRI Reconstruction that jointly reconstructs two contrasts, 
% e.g. T1-weighted and T2-weighted contrasts from their under-sampled k-space measurements.
%Note that all input parameters need to be set prior to simulation. We provide some example settings for the input parameters below. However, the user is
%advised to choose optimal values for the parameters depending on the specific data or task at hand.
% 
% 
% Inputs -
%       Img.I1A and I1B : Input fully-sampled target and guidance MR Image (real valued, non-negative). 
%       QA, QB : Sampling Mask for 2D DFT data (Center of k-space corresponds to corners of mask Q1) with zeros at non-sampled locations and ones at sampled locations.
%       sigmaiA, sigmaiB : Noise level (standard deviation of complex noise) in the DFT space of the peak-normalized input image. To be set to 0 if input image data is assumed noiseless.
%       sigmaA, sigmaB : Simulated noise level (standard deviation of simulated complex noise - added in k-space during simulation).
%       DLMRIparams: Structure that contains the parameters of the CDLMRI algorithm. The various fields are as follows - 
%                   - num: Number of iterations of the CDLMRI algorithm (Example: about 15)
%                   - n: Patch size, i.e., Total # of pixels in square patch (e.g., n=36, n=64)
%                   - K: Number of dictionary atoms (e.g., K=n)
%                   - N: Number of signals used for training (e.g., N=200*K)
%                   - T0: Common sparsity settings of patch (e.g., T0=round((0.15)*n) )
%                   - Lambda: Sets the weight \nu in the algorithm (e.g., Lambda=140)
%                   - th: Thresholds used in sparse coding.
%                   - MAX_ITER: Number of iterations within the K-SVD algorithm. (e.g., numiterateKSVD = 20, 10)
%                   - weights: weights for exploiting guidance information;
%                   - variance_Thresh: Discard those patch pairs with too small variance during training.
%                   - r; distance between two adjacent image patches (e.g. r=1).
%					- fold: subsampling fold computed according to the sampling mask (e.g. fold = 10).
%					- lambda: lagrange multipliers for Lasso problem.
%					- s: sparsity constraints.
%					- weights: weights for guidance/reference info (e.g. 0.5).
%					- variance_Thresh: Discard those patch pairs with too small variance during training (e.g. 0.1).
%
% Outputs -
%       Iout1.A - reconstructed target MR image.
%       Iout1.B - reconstructed guidance MR image.
%       param1 - Structure containing various performance metrics, etc. from simulation for CDLMRI.
%                 - InputPSNR : PSNR of fully sampled noisy reconstruction after adding simulated noise.
%                 - PSNR0 : PSNR of normalized zero-filled reconstruction
%                 - PSNR : PSNR of the reconstruction at each iteration
%                 - HFEN : HFEN of the reconstruction at each iteration
%                 - itererror : norm of the difference between the reconstructions at successive iterations
%                 - Dictionary : Final dictionary
%
% References:
% ------------
% [1] P. Song, L. Weizman, J. F. C. Mota, Y. C. Eldar and M. R. D. Rodrigues, "Coupled Dictionary Learning for Multi-Contrast MRI Reconstruction," in IEEE Transactions on Medical Imaging, vol. 39, no. 3, pp. 621-633, March 2020, doi: 10.1109/TMI.2019.2932961.
% [2] P. Song, L. Weizman, J. F. C. Mota, Y. C. Eldar and M. R. D. Rodrigues, "Coupled Dictionary Learning for Multi-Contrast MRI Reconstruction," IEEE International Conference on Image Processing (ICIP), 2018, pp. 2880-2884, doi: 10.1109/ICIP.2018.8451341.
%
% Codes written & compiled by:
% ----------------------------
% Pingfan Song
% Electronic and Electrical Engineering,
% Imperial College London
% p.song@imperial.ac.uk
%
% =========================================================================


addpath('./Dicts');
addpath(genpath('./utils'))

clear;

% load and adjust true images
NUMROWS = 256;
NUMCOLS = 256;

% Construct testing dataset.
XpathCell = glob('./TestImages_MRI', '*_T1.png' );
Xcell = load_images( XpathCell );
for i =1 : length(Xcell)
	Xcell{i} = imresize(Xcell{i}, [NUMROWS, NUMCOLS], 'bicubic');
end

YpathCell = glob('./TestImages_MRI', '*_T2.png' );
Ycell = load_images(YpathCell );
for i =1 : length(Ycell)
	Ycell{i} = imresize(Ycell{i}, [NUMROWS, NUMCOLS], 'bicubic');
end
% 	figure; subplot(1,2,1); imshow(Xcell{1}); subplot(1,2,2); imshow(Ycell{1});
if length(Xcell) ~= length(Ycell)	
	error('Error: The number of X images is not equal to the number of Y images!');
end

% Normalize and truncate interpolated images in order to limit the pixel value between 0 and 1.
for i =1 : length(Xcell)
	Img = Xcell{i} ;
	Img = Img./max(abs(Img(:)));
	Img( Img > 1) = 1 ;
	Img( Img < 0) = 0 ;		
	Xcell{i} = Img ;

	Img = Ycell{i} ;
	Img = Img./max(abs(Img(:)));
	Img( Img > 1) = 1 ;
	Img( Img < 0) = 0 ;		
	Ycell{i} = Img ;		
end	
clear Img
DLMRIparams = [];

% load pre-learned dictionaries.
CDLname = 'CDL_T1T2_D64x512_WT0.5.mat';
load( CDLname );

n = paramsCDL.N;
K = paramsCDL.K;

% Trained dictionaries
Psi_cx = outputCDL.Psi_cx;
Psi_x = outputCDL.Psi_x;
Psi_cy = outputCDL.Psi_cy;
Psi_y = outputCDL.Psi_y;

Dict = [Psi_cx, Psi_x, zeros(n,K) ; ...
			Psi_cy, zeros(n,K), Psi_y]; 
								
% load masks
QA = load('./SamplingMasks/5fold_Size256/Q1_5fold_p3.mat');
QB = load('./SamplingMasks/5fold_Size256/Q2_5fold_p3.mat');
QA = ( QA.sampling_mask );
QA = ifftshift( QA );
QB = ( QB.sampling_mask );
QB = ifftshift( QB );

Q.Q1A = ( QA );
Q.Q1B = ( QB );
indexQA=find(QA==1); %Index the sampled locations in sampling mask
fold = round( (NUMROWS * NUMCOLS) / length(indexQA) );


% set noise power
sigmai.sigmaiA = 0 ; % Noise level for input image.
sigmai.sigmaiB = 0 ;
sigma.sigmaA = 0 ; % Noise level for k-space sampling.
sigma.sigmaB = 0 ;

% set parameters
numIter = 40; % number of outer iterations.
ErrTh_start = 0.03; % error threshold (RMSE) in the beginning, % RMSE = 10^(-PSNR_goal/20);
ErrTh_end = 0.01;  % error threshold at the end	
errNum = 40; % number of outer iterations with error threshold decreasing linearly from ErrTh_start to ErrTh_end.
errthArray = linspace(ErrTh_start, ErrTh_end, errNum); % linearly decreasing error thresholds.
errthArray = [errthArray, ErrTh_end*ones(1, numIter-errNum)];

DLMRIparams.num = numIter; % number of outer iterations.
DLMRIparams.n = n; % length of an atom, i.e., the number of pixels in an image patch.
DLMRIparams.K = K; % number of atoms in a dictionary.
DLMRIparams.N = 200*DLMRIparams.K;  % number of samples for dictionary training.
DLMRIparams.T0 = 10; % sparsity constraint, round((0.15)*DLMRIparams.n);
DLMRIparams.Lambda = 140; % a parameter for noisy case.
DLMRIparams.opt = 1; % If set to 1, CDL is done with fixed sparsity.
DLMRIparams.th = errthArray; % error thresholds for sparse coding during reconstruction.
DLMRIparams.MAX_ITER = 50; % number of inner iterations for dictionary learning
DLMRIparams.r = 1; % distance between two adjacent image patches.
DLMRIparams.fold = fold ; % subsampling fold computed according to the sampling mask.
DLMRIparams.XpathCell = XpathCell ; % image path
DLMRIparams.YpathCell = YpathCell ; % image path

DLMRIparams.lambda = [0.1, 0.05, 0.05] ; % lagrange multipliers for Lasso problem.
DLMRIparams.s = [DLMRIparams.T0, 1, 1] ; % sparsity constraints
DLMRIparams.weights = 0.5; % weights for guidance/reference info.
DLMRIparams.variance_Thresh = 0.1 ; % Discard those patch pairs with too small variance during training.
DLMRIparams.Psi_cx0 = outputCDL.Psi_cx(:,1:K); % initial dictionaries
DLMRIparams.Psi_x0 = outputCDL.Psi_x(:,1:K);
DLMRIparams.Psi_cy0 = outputCDL.Psi_cy(:,1:K);
DLMRIparams.Psi_y0 = outputCDL.Psi_y(:,1:K);
		

% % for debug only
% DLMRIparams.num = 5;
% DLMRIparams.MAX_ITER = 20;
% DLMRIparams.r = 1; %1


% file name
SIZE = ['_D',num2str( DLMRIparams.n ),'x',num2str( DLMRIparams.K )];
MaxIter = ['_Iter', num2str( DLMRIparams.num )];
% TrainSize = ['_T', num2str( DLMRIparams.N )]; % no need
current_date = date;
DATE = ['_Date',current_date];
S = ['_s', num2str( DLMRIparams.T0 )];
Fold = ['_', num2str( DLMRIparams.fold ), 'fold'];
PatchDist = ['_PatchDist', num2str( DLMRIparams.r )];
W = ['_W', num2str( DLMRIparams.weights )];
FILENAME = ['CDLMRI_T1T2', SIZE, Fold, S, MaxIter, W, DATE];


PSNR_all = [];
S_train_all = [];
S_test_all = [];
for i =1 : length(Xcell)

	fprintf('Processing image No. %d ... \n', i);
	
	Img.I1A = Xcell{i};
	Img.I1B = Ycell{i};
	
	%% The first stage uses maximum overlapping setting (stride=1) to perform CDLMRI
	% do CDLMRI
	[Iout1,param1] = CDLMRI(Img, Q, sigmai,sigma,DLMRIparams); 
% 	Iout1.A = uint8( 255.*abs(Iout1.A) ) ;
% 	Iout1.B = uint8( 255.*abs(Iout1.B) ) ;
	
	Result_cell{i}.Img = Img;
	Result_cell{i}.Iout1 = Iout1;
	Result_cell{i}.param1 = param1;
	
	PSNR_all(i, 1) = param1.PSNRA(end);
	PSNR_all(i, 2) = param1.PSNRB(end);
	S_train_all(i, :) = param1.S_mean_iter(end, :);
	S_test_all(i, :) = param1.S_mean_test_iter(end, :);
	
	
	Output.PSNR_all = PSNR_all;
	Output.S_train_all = S_train_all;
	Output.S_test_all = S_test_all;
	Output.PSNR_mean = mean(PSNR_all,1) ;
	Output.S_train_mean = mean(S_train_all,1) ;
	Output.S_test_mean = mean(S_test_all,1) ;

	% save data
	save([FILENAME,'.mat'], ...
		'DLMRIparams', 'sigmai', 'sigma', 'Q', 'Result_cell',  'Output' )

	%% The second stage uses smaller overlapping setting (e.g. stride=6) to refine the results. 
	% set parameters
	numIter = 10;
	ErrTh_start = 0.01; % error threshold (RMSE) in the beginning, % RMSE = 10^(-PSNR_goal/20);
	ErrTh_end = 0.001;  % error threshold at the end	
	errNum = 10;
	errthArray = linspace(ErrTh_start, ErrTh_end, errNum); % error thresholds.
	errthArray = [errthArray, ErrTh_end*ones(1, numIter-errNum)];
	
	
	DLMRIparams2 = DLMRIparams;
	DLMRIparams2.num = numIter; % number of outer iterations.
	DLMRIparams2.th = errthArray; % error thresholds for sparse coding during reconstruction.
	DLMRIparams2.r = 6; % distance between two adjacent image patches.
	DLMRIparams2.initImg = abs(Iout1.A);
	DLMRIparams2.initImgB = abs(Iout1.B);
	
	% do CDLMRI
	[Iout2,param2] = CDLMRI(Img, Q, sigmai,sigma, DLMRIparams2); 
	
	Result_cell{i}.Iout2 = Iout2;
	Result_cell{i}.param2 = param2;
	
	PSNR_all2(i, 1) = param2.PSNRA(end);
	PSNR_all2(i, 2) = param2.PSNRB(end);
	S_train_all2(i, :) = param2.S_mean_iter(end, :);
	S_test_all2(i, :) = param2.S_mean_test_iter(end, :);

	Output.PSNR_all2 = PSNR_all2;
	Output.S_train_all2 = S_train_all2;
	Output.S_test_all2 = S_test_all2;
	Output.PSNR_mean2 = mean(PSNR_all2,1) ;
	Output.S_train_mean2 = mean(S_train_all2,1) ;
	Output.S_test_mean2 = mean(S_test_all2,1) ;

	% save data
	save([FILENAME,'.mat'], ...
		'DLMRIparams', 'DLMRIparams2', 'sigmai', 'sigma', 'Q', 'Result_cell',  'Output' )

end

%% disp information for all images
for i = 1 : numel(Result_cell)
	info = sprintf('Image No. %d:\n', i) ;
	info_part1 = Result_cell{i}.param1.info_all ;
	info_part2 = Result_cell{i}.param2.info_all ;
	
	fprintf([ info, 'Stage 1: \n',  info_part1, 'Stage 2: \n', info_part2 ])
end



