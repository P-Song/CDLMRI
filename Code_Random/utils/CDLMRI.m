function [Iout1,param1] = CDLMRI(Img, Q, sigmai,sigma,DLMRIparams)
% =========================================================================
% Coupled Dictionary Learning for Multi-contrast MRI Reconstruction
% =========================================================================
% 
%The software is used to perform multi-contrast MRI reconstruction that jointly reconstructs two contrast, 
% e.g. T1-weighted and T2-weighted contrasts from their under-samled k-space measurements.
%Note that all input parameters need to be set prior to simulation. We provide some example settings for the input parameters below. However, the user is
%advised to choose optimal values for the parameters depending on the specific data or task at hand.
% 
% 
% Inputs -
%       Img.I1A and I1B : Input fully-sampled target and guidance MR Image (real valued, non-negative). 
%       Q1A, Q1B : Sampling Mask for 2D DFT data (Center of k-space corresponds to corners of mask Q1) with zeros at non-sampled locations and ones at sampled locations.
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



% Parameters initialization

Q1A = Q.Q1A ;
Q1B = Q.Q1B ;

I1A = Img.I1A;
I1B = Img.I1B;

sigmaiA = sigmai.sigmaiA ;
sigmaiB = sigmai.sigmaiB ;
sigmaA = sigma.sigmaA ;
sigmaB = sigma.sigmaB ;

sigma2A=sqrt((sigmaiA^2)+(sigmaA^2));  %Effective k-space noise level. 
sigma2B=sqrt((sigmaiB^2)+(sigmaB^2));  %Effective k-space noise level. 
%If sigma2A is very small (sigma2A<1e-10), then the sampled k-space locations are filled back without averaging in the reconstruction update step 
%of our algorithm.

num=DLMRIparams.num;  % iteration count
Lambda=DLMRIparams.Lambda; %Lambda parameter
La2A=(Lambda)/(sigma2A); % \nu weighting of paper
La2B=(Lambda)/(sigma2B); % \nu weighting of paper
r=DLMRIparams.r; %Overlap Stride
n=DLMRIparams.n; %patch size
K=DLMRIparams.K; %number of dictionary atoms
N=DLMRIparams.N; %number of training signals
T0=DLMRIparams.T0; %sparsity levels
th=DLMRIparams.th; %error threshold for patches - allows error of (th^2)*n per patch during sparse coding.
op=DLMRIparams.opt;  %Type of learning
MAX_ITER=DLMRIparams.MAX_ITER; %number of K-SVD iterations
weights = DLMRIparams.weights ;
variance_Thresh = DLMRIparams.variance_Thresh ; % Discard those patch pairs with too small variance during training.

s_c = T0;
s_x = DLMRIparams.s(2);
s_y = DLMRIparams.s(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN SIMULATION CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------for modal A-----------
I1A=double(I1A(:,:,1));
I1A=I1A/(max(max(I1A))); %Normalize input image - peak reconstruction pixel value restricted to 1 in this simulation.
[aa,bb]=size(I1A);     %Compute size of image

DZA=((sigmaA/sqrt(2))*(randn(aa,bb)+(0+1i)*randn(aa,bb)));  %simulating noise
I5A=fft2(I1A);          %FFT of input image
I5A=I5A+DZA;             %add measurement noise in k-space (simulated noise)

%-----------------for modal B-----------
I1B=double(I1B(:,:,1));
I1B=I1B/(max(max(I1B))); %Normalize input image - peak reconstruction pixel value restricted to 1 in this simulation.
% [aa,bb]=size(I1B);     %Compute size of image

DZB=((sigmaB/sqrt(2))*(randn(aa,bb)+(0+1i)*randn(aa,bb)));  %simulating noise
I5B=fft2(I1B);          %FFT of input image
I5B=I5B+DZB;             %add measurement noise in k-space (simulated noise)


%-----------------for modal A-----------
%Compute Input PSNR after adding simulated noise
IGA=abs(ifft2(I5A));
InputPSNRA=20*log10((sqrt(aa*bb))/norm(double(abs(IGA))-double(I1A),'fro'));
param1.InputPSNRA=InputPSNRA;

indexA=find(Q1A==1); %Index the sampled locations in sampling mask

I2A=(double(I5A)).*(Q1A);  %Apply mask in DFT domain
I11A=ifft2(I2A);          % Inverse FFT - gives zero-filled result

%initializing simulation metrics
itererrorA=zeros(num,1);
highfritererrorA=zeros(num,1);
PSNR1A=zeros(num,1);
RMSEA = zeros(num,1);
mssimA = zeros(num,1);

I11A=I11A./(max(max(abs(I11A))));  %image normalization

% use pre-recovered results as the initial image.
if isfield(DLMRIparams, 'initImg')
	I11A = DLMRIparams.initImg ;
	I11A = im2double(I11A);
% 	I11A=I11A./(max(abs(I11A(:))));  %image normalization
end

PSNR0A=20*log10(sqrt(aa*bb)/norm(double(abs(I11A))-double(I1A),'fro')); %PSNR of normalized zero-filled reconstruction
param1.PSNR0A=PSNR0A;

%-----------------for modal B-----------
%Compute Input PSNR after adding simulated noise
IGB=abs(ifft2(I5B));
InputPSNRB=20*log10((sqrt(aa*bb))/norm(double(abs(IGB))-double(I1B),'fro'));
param1.InputPSNRB=InputPSNRB;

indexB=find(Q1B==1); %Index the sampled locations in sampling mask

I2B=(double(I5B)).*(Q1B);  %Apply mask in DFT domain
I11B=ifft2(I2B);          % Inverse FFT - gives zero-filled result

%initializing simulation metrics
itererrorB = zeros(num,1);
highfritererrorB = zeros(num,1);
PSNR1B = zeros(num,1);
RMSEB = zeros(num,1);

I11B=I11B/(max(max(abs(I11B))));  %image normalization

% use pre-recovered results as the initial image.
if isfield(DLMRIparams, 'initImgB')
	I11B = DLMRIparams.initImgB ;
	I11B = im2double(I11B);
end

PSNR0B=20*log10(sqrt(aa*bb)/norm(double(abs(I11B))-double(I1B),'fro')); %PSNR of normalized zero-filled reconstruction
param1.PSNR0B=PSNR0B;

constr_array = []; % store constraints in each iteration.
S_mean_test_iter = [] ;
S_mean_iter = [] ;

TSTART = tic; % measure time.
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%CDLMRI iterations
for kp=1:num
    
    I11A=abs(I11A); % I11A=I11A/(max(max(I11A)));
    IiterA=I11A;
	
	I11B=abs(I11B); % I11B=I11B/(max(max(I11B)));
    IiterB=I11B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Creating image patches
	
	%-----------------for modal A-----------
    [blocksA,idxA] = my_im2col(I11A,[sqrt(n),sqrt(n)],r); 
	brA=mean(blocksA); %image patches
    br2A = (ones(n,1))*brA;
    TEA=blocksA-br2A;           %subtract means of patches
    [rowsA,colsA] = ind2sub(size(I11A)-sqrt(n)+1,idxA);
	
    %-----------------for modal B-----------
    [blocksB,idxB] = my_im2col(I11B,[sqrt(n),sqrt(n)],r); 
	brB=mean(blocksB); %image patches
    br2B = (ones(n,1))*brB;
    TEB=blocksB-br2B;           %subtract means of patches
    [rowsB,colsB] = ind2sub(size(I11B)-sqrt(n)+1,idxB);
	
	N2A=size(blocksA,2); %total number of overlapping image patches
	N2B=size(blocksB,2); %total number of overlapping image patches
	if N2A ~= N2B
		error('CDLMRI Error 1: The size of training data A is not equal to the size of training data B!');
	end
    de=randperm(N2A);
    
    %Check if specified number of training signals is less or greater than the available number of patches.
    if(N2A>N)
        N4=N;
    else
        N4=N2A;
	end

    YHA=TEA(:,de(1:N4));   %Training data - using random selection/subset of patches
    YHB=TEB(:,de(1:N4));   %Training data - using random selection/subset of patches
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %dictionary initialization : PCA + random training patches
	
	%-----------------for modal A-----------
    [UU,SS,VV]=svd(YHA*YHA');
    D0A=zeros(n,K);
    [hh,jj]=size(UU);
    D0A(:,1:jj)=UU;
    p1=randperm(N4);
    for py=jj+1:K
        D0A(:,py)=YHA(:,p1(py-jj));
    end
    paramA.initialDictionaryA=D0A;   %initial dictionary for K-SVD algorithm
	

%-----------------for modal B-----------
	[UU,SS,VV]=svd(YHB*YHB');
    D0B=zeros(n,K);
    [hh,jj]=size(UU);
    D0B(:,1:jj)=UU;
    p1=randperm(N4);
    for py=jj+1:K
        D0B(:,py)=YHB(:,p1(py-jj));
    end
    paramB.initialDictionaryB=D0B;   %initial dictionary for K-SVD algorithm
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     coupled dictionary learning algorithm 
    
    if(op==1)
		
		% retrain the dicts on the testing images. 
		% Discard those patch pairs with too small variance when construct traning dataset.
		Xnorm2 = sum(YHA.^2, 1);
		X_index = (Xnorm2 > variance_Thresh);
		Ynorm2 = sum(YHB.^2, 1);
		Y_index = (Ynorm2 > variance_Thresh);
		XY_index = X_index | Y_index;
		TrainData.X = YHA(:, XY_index);
		TrainData.Y = YHB(:, XY_index);
		
		errth = th(kp);
		DLMRIparams.errth = errth ;	 
		[outputCDL] = CDL_MRI_ODL_OMP(TrainData, DLMRIparams);	 
		% Trained dictionaries
		Psi_cx = outputCDL.Psi_cx;
		Psi_x = outputCDL.Psi_x;
		Psi_cy = outputCDL.Psi_cy;
		Psi_y = outputCDL.Psi_y;
		infoCDL = outputCDL(1).infoCDL;

		GammaAB = zeros(3*K, size(YHA,2)) ;
		
		TrainErr(:, kp) = mean( outputCDL(1).Res(end, :) ) ;  % store train error
		Ratio(kp) = {[]} ; 
		Dist(kp) = {[]} ; 
		
		ind1 = 1: K;
		ind2 = (K+1) : 2*K;
		ind3 = (2*K +1) : 3*K;

		ind1_row = 1: n;
		ind2_row = (n + 1) : 2*n;
	
		Zc = GammaAB(ind1,:);
		Zx = GammaAB(ind2,:);
		Zy = GammaAB(ind3,:);	

		S_array(1, : ) = sum(abs(Zc)>0.001);
		S_array(2, : ) = sum(abs(Zx)>0.001);
		S_array(3, : ) = sum(abs(Zy)>0.001);
		S_mean(1,1) = mean(S_array(1, : ) );
		S_mean(1,2) = mean(S_array(2, : ) );
		S_mean(1,3) = mean(S_array(3, : ) );
		
		S_mean_iter = [S_mean_iter; full(S_mean)] ;
    end
    
    %DLerror=(norm(YHA -(DA*output.CoefMatrix),'fro'))^2  %dictionary fitting error - uncomment to monitor.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Computing sparse representations of all patches and summing up the patch approximations
	
	%-----------------for both modal A and B-----------
    WeightA= zeros(aa,bb);
    IMoutA = zeros(aa,bb); 
	bbb=sqrt(n);
	
	WeightB= zeros(aa,bb);
	IMoutB = zeros(aa,bb); 
	
	S_mean_test = [] ;

    for jj = 1:10000:size(blocksA,2)
		
        jumpSize = min(jj+10000-1,size(blocksA,2)); 
		ZZA=TEA(:,jj:jumpSize);
		ZZB=TEB(:,jj:jumpSize);
		
		ZZAB = [ZZA ; ZZB];

		gain = 1; % 1.2; %1.15;
		D = [Psi_cx; Psi_cy];
		normMat = diag(1./sqrt(sum(D.^2)) ) ; D = D*normMat; 
		Psi_cx2 = Psi_cx*normMat; Psi_cy2 = Psi_cy*normMat; 
% 		D = D./ repmat(sqrt(sum(D.^2)) , [size(D ,1), 1]) ; 
		D0  =D;		

		data = [ZZA ; weights.*ZZB];
		errth = th(kp);

		% Use OMP in spams
		paramOMP.L=s_c; % not more than xx non-zeros coefficients
		paramOMP.eps= (sqrt(n) * errth * gain)^2; 
		paramOMP.numThreads = -1; % number of threads
		Zc=mexOMP(data, D, paramOMP);		

		X_lowR = ZZA - Psi_cx2*Zc;
		Y_highR = weights.*ZZB - Psi_cy2*Zc;

		paramOMP.L=s_x; 
		Zx = mexOMP(X_lowR, Psi_x, paramOMP); 
		paramOMP.L=s_y; 
		Zy = mexOMP(Y_highR, Psi_y, paramOMP); 
		
		CoefsAB = [Zc; Zx; Zy];
	
		ZZA_noM = Psi_cx2 * Zc + Psi_x * Zx;  % Recovery without mean;
		ZZB_noM = Psi_cy2 * Zc + Psi_y * Zy;  % Recovery without mean;
		ZZB_noM = ZZB_noM./weights;
		
        ZZA= ZZA_noM + (ones(size(blocksA,1),1) * brA(jj:jumpSize)); %sparse approximations of patches
        ZZB= ZZB_noM + (ones(size(blocksB,1),1) * brB(jj:jumpSize)); %sparse approximations of patches

		
		S_array_test = [];
		S_array_test(1, : ) = sum(abs(Zc)>0.001);
		S_array_test(2, : ) = sum(abs(Zx)>0.001);
		S_array_test(3, : ) = sum(abs(Zy)>0.001);
		
		s_c_batch = mean(S_array_test(1, : ) );
		s_x_batch = mean(S_array_test(2, : ) );
		s_y_batch = mean(S_array_test(3, : ) );
		
		S_mean_test = [S_mean_test ; [ s_c_batch, s_x_batch, s_y_batch] ] ;
		
        %summing up patch approximations
        for i  = jj:jumpSize
            col = colsA(i); 
			row = rowsA(i);
            block =reshape(ZZA(:,i-jj+1),[bbb,bbb]);
            IMoutA(row:row+bbb-1,col:col+bbb-1)=IMoutA(row:row+bbb-1,col:col+bbb-1)+block;
            WeightA(row:row+bbb-1,col:col+bbb-1)=WeightA(row:row+bbb-1,col:col+bbb-1)+ones(bbb);
		end
		
		 for i  = jj:jumpSize
            col = colsB(i); 
			row = rowsB(i);
            block =reshape(ZZB(:,i-jj+1),[bbb,bbb]);
            IMoutB(row:row+bbb-1,col:col+bbb-1)=IMoutB(row:row+bbb-1,col:col+bbb-1)+block;
            WeightB(row:row+bbb-1,col:col+bbb-1)=WeightB(row:row+bbb-1,col:col+bbb-1)+ones(bbb);
		 end

    end
	
    I3nA=IMoutA./WeightA;  %patch-averaged result
    innA= abs(I3nA)>1;
	I3nA(innA)=1;
	
	I3nB=IMoutB./WeightB;  %patch-averaged result
    innB= abs(I3nB)>1;
	I3nB(innB)=1;	
	
	S_mean_test_iter(kp, :) = mean(S_mean_test, 1) ;
	constr_array(kp, :) = [errth, s_c] ;
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%-----------------for modal A-----------
    I2A=fft2(I3nA);  %Move from image domain to k-space
    
    %K-space update formula
    if(sigma2A<1e-10)
        I2A(indexA)=I5A(indexA);
    else
        I2A(indexA)= (1/(1+(La2A)))*(I2A(indexA) + (La2A)*I5A(indexA));
    end
    
    I11A=ifft2(I2A);   %Use Inverse FFT to get back to image domain
    inn2A= abs(I11A)>1;
	I11A(inn2A)=1;
	
	%-----------------for modal B-----------
    I2B=fft2(I3nB);  %Move from image domain to k-space
    
    %K-space update formula
    if(sigma2B<1e-10)
        I2B(indexB)=I5B(indexB);
    else
        I2B(indexB)= (1/(1+(La2B)))*(I2B(indexB) + (La2B)*I5B(indexB));
    end
    
    I11B=ifft2(I2B);   %Use Inverse FFT to get back to image domain
    inn2B= abs(I11B)>1;
	I11B(inn2B)=1;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Compute various performance metrics
	
	%-----------------for modal A-----------
    itererrorA(kp)= norm(abs(IiterA) - abs(I11A),'fro');
    highfritererrorA(kp)=norm(imfilter(abs(I11A),fspecial('log',15,1.5)) - imfilter(I1A,fspecial('log',15,1.5)),'fro');
    PSNR1A(kp)=20*log10(sqrt(aa*bb)/norm(double(abs(I11A))-double(I1A),'fro'));
	RMSEA(kp) = norm(double(abs(I11A))-double(I1A),'fro')/sqrt(aa*bb) ;
   
	%-----------------for modal B-----------
	itererrorB(kp)= norm(abs(IiterB) - abs(I11B),'fro');
    highfritererrorB(kp)=norm(imfilter(abs(I11B),fspecial('log',15,1.5)) - imfilter(I1B,fspecial('log',15,1.5)),'fro');
    PSNR1B(kp)=20*log10(sqrt(aa*bb)/norm(double(abs(I11B))-double(I1B),'fro'));
	RMSEB(kp) = norm(double(abs(I11B))-double(I1B),'fro')/sqrt(aa*bb) ;
    
	% 	SSIM
	[mssimA(kp), ssim_mapA] = ssim(abs(I11A), I1A) ;
	
	% store all the information
	info_str = ['OutIter = %d; PSNR(T1,T2) = (%.3f, %.3f); RMSE_Test(T1,T2) = (%.4f, %.4f);\n'] ;
	info_num = [(1:kp)' , PSNR1A(1:kp), PSNR1B(1:kp), RMSEA(1:kp), RMSEB(1:kp)]' ;
	info_all = sprintf(info_str , info_num) ;
	disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	disp(info_all) ;
	disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

	if mod(kp,5) == 0
		figure(101)
		imagesc(abs([I11A, I11B])); colormap gray
		title('Target T1 / Guidance T2')
		pause(0.2)
	end
	
	% save data
	save(['tempCDLMRI','.mat'], ...
		'DLMRIparams', 'I11A', 'PSNR1A', 'itererrorA' , ...
		 'I11B', 'PSNR1B', 'itererrorB' , 'TrainErr', 'info_all')

end
time = toc(TSTART);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

Iout1.A=abs(I11A);
param1.PSNRA=PSNR1A;
param1.HFENA=highfritererrorA;
param1.itererrorA=itererrorA;
param1.outputCDL=outputCDL;

Iout1.B=abs(I11B);
param1.PSNRB=PSNR1B;
param1.HFENB=highfritererrorB;
param1.itererrorB=itererrorB;
param1.TrainErr=TrainErr;

param1.S_mean_test_iter=S_mean_test_iter;
param1.S_mean_iter=S_mean_iter;
param1.info_all=info_all;
param1.constr_array=constr_array;

param1.RMSEA=RMSEA;
param1.RMSEB=RMSEB;
param1.time=time;
param1.mssimA=mssimA;

end



