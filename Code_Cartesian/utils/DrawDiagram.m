% 18/12/2017 generate zero-filled MRI images for T1 and T2. Draw successive
% images
% 30/10/2017 generate zero-filled MRI images for T1
% close all;
% addpath('./Dicts');
% addpath('./utils')
% %% Extract a few atoms and use them to synthesize patches to demonstrate sparse coding
% % load trained dicts
% load('G:\Work\multi-modality super-resolution\CDLMRI_VS_DLMRI\Code_T1T2_EnforceCom_Cartesian\Dicts\CDL_T1T2_D64x512_WT0.5_LamC0.1_Iter1000_T53856_Date13-Sep-2017.mat')
% 
% blocksize = [8,8];
% [~, index] = sort( sum(Psi_cx.^2) );
% Psi_cx = Psi_cx(:,index);
% Psi_cy = Psi_cy(:,index);
% 
% [~, index] = sort( sum(Psi_x.^2) );
% Psi_x = Psi_x(:,index);
% 
% [~, index] = sort( sum(Psi_y.^2) );
% Psi_y = Psi_y(:,index);	
% 
% % extract a subset of atoms from common dicts
% Row = 23;
% start = Row * 4 + 4; 
% width = 3;
% BasicRange = start: (start+width-1);
% Psi_cx2= Psi_cx(:,[BasicRange, BasicRange+Row, BasicRange+Row*2]);
% Psi_cy2= Psi_cy(:,[BasicRange, BasicRange+Row, BasicRange+Row*2]);
% 
% % extract a subset of atoms from unique dicts
% start = Row * 6 + 4; 
% width = 3;
% BasicRange = start: (start+width-1);
% Psi_x2= Psi_x(:,[BasicRange, BasicRange+Row, BasicRange+Row*2]);
% Psi_y2= Psi_y(:,[BasicRange, BasicRange+Row, BasicRange+Row*2]);
% 
% 
% % synthesize patches use specified atoms
% CXatom1 = 7; CXatom2 = 9; Xatom1 = 9; 
% CYatom1 = 7; CYatom2 = 9; Yatom1 = 2; 
% 
% d1 = 1 * Psi_cx2(:,CXatom1) + 0.5 * Psi_cx2(:,CXatom2) - 0.2 * Psi_x2(:,Xatom1);
% d2 = 1 * Psi_cy2(:,CYatom1) + 0.5 * Psi_cy2(:,CYatom2) + 0.1 * Psi_y2(:,Yatom1);
% 
% % reshape the vectorized patches to blocks
% N = size(Psi_cx, 1);
% K = size(Psi_cx2, 2);
% 
% d1s = col2im(d1,blocksize,blocksize*1,'distinct');
% d2s = col2im(d2,blocksize,blocksize*1,'distinct');
% 
% % show dicts
% addpath('./utils')
% addpath('../export_fig')
% plot_enable = 1;
% FigFontName = 'Times New Roman';
% FigFontSize = 10;
% SaveFig = 0 ; % save figure or not
% FigFormatCell = { '.png', '.fig'}; %.eps
% 
% blocksize = [sqrt(N), sqrt(N)];
% 	
% dictimg1 = SMALL_showdict(Psi_cx2,blocksize,	round(sqrt( K )), round(sqrt( K )) ,'lines','highcontrast');  
% dictimg2 = SMALL_showdict(Psi_x2,blocksize, round(sqrt( K )), round(sqrt( K )) ,'lines','highcontrast');  
% dictimg3 = SMALL_showdict(Psi_cy2,blocksize, round(sqrt( K )), round(sqrt( K )) ,'lines','highcontrast');  
% dictimg4 = SMALL_showdict(Psi_y2,blocksize, round(sqrt( K )), round(sqrt( K )) ,'lines','highcontrast');
% 
% if plot_enable
% 	
% %--------------------------------------------------------
% 	% Display dictionary atoms as image patch
% 	Hcf = figure ;
% 	Hca = gca ;
% 
% 	subplot(2,2,1)
% 	imagesc(dictimg1);colormap(gray);axis off; axis image;
% 	title('Dict Psi\_{cx}')
% 
% 	subplot(2,2,2)
% 	imagesc(dictimg2);colormap(gray);axis off; axis image;
% 	title('Dict Psi\_{x}')
% 
% 	subplot(2,2,3)
% 	imagesc(dictimg3);colormap(gray);axis off; axis image;
% 	title('Dict Psi\_{cy}')
% 
% 	subplot(2,2,4)
% 	imagesc(dictimg4);colormap(gray);axis off; axis image;
% 	title('Dict Psi\_{y}')
% 	
% 	if SaveFig
% 		FigName = ['DictsSubset'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 	
% 
% %--------------------------------------------------------
% 	% display synthesized patches	
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(d1s); colormap gray;  axis off;
% 	set(gcf, 'position', [100,400,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['PatchX'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 	
% 	Hcf = figure ;	Hca = gca ;
% 	imagesc(d2s);colormap gray; axis off;
% 	set(gcf, 'position', [100, 100,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 		
% 	if SaveFig
% 		FigName = ['PatchY'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 	
% %--------------------------------------------------------
% 	% display specified atoms for common dicts
% 	Hcf = figure ; 	Hca = gca ; 
% 	imagesc(col2im(Psi_cx2(:,CXatom1),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [300,400,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['CXatom1'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(col2im(Psi_cx2(:,CXatom2),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [500,400,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['CXatom2'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 	
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(col2im(Psi_x2(:,Xatom1),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [700,400,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['Xatom1'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 
% 	% display specified atoms for unique dicts
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(col2im(Psi_cy2(:,CYatom1),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [300,100,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['CYatom1'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(col2im(Psi_cy2(:,CYatom2),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [500,100,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['CYatom2'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 	
% 	Hcf = figure ; 	Hca = gca ;
% 	imagesc(col2im(Psi_y2(:,Yatom1),blocksize,blocksize,'distinct') ); colormap gray;  axis off;
% 	set(gcf, 'position', [700,100,200,200])
% 	set(gca, 'position', [0,0,1,1])
% 	set(gcf, 'Color', 'w'); % make the background to be white.
% 	pause(0.1)
% 	
% 	if SaveFig
% 		FigName = ['Yatom1'];
% 		for i = 1 : numel(FigFormatCell)
% 			export_fig(Hcf, [FigName, FigFormatCell{i}]) % export_fig(figure_handle, filename);
% 		end
% 		pause(0.1)
% 	end
% 
% end
% 
% disp('Done!')
% close all


%% show zero-filled image and spectral


addpath('../TestImages_MRI');
addpath('./utils')

clear;

% load and adjust true images
NUMROWS = 256;
NUMCOLS = 256;

% Construct testing dataset.
XpathCell = glob('../TestImages_MRI', '*_T1.png' );
Xcell = load_images( XpathCell );
for i =1 : length(Xcell)
	Xcell{i} = imresize(Xcell{i}, [NUMROWS, NUMCOLS], 'bicubic');
end

YpathCell = glob('../TestImages_MRI', '*_T2.png' );
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
CDLname = 'CDL_T1T2_D64x512_WT0.5_LamC0.1_Iter1000_T53856_Date13-Sep-2017.mat';
load( CDLname );

n = paramsCDL.N;
K = 512;

% Trained dictionaries
Psi_cx = outputCDL.Psi_cx;
Psi_x = outputCDL.Psi_x;
Psi_cy = outputCDL.Psi_cy;
Psi_y = outputCDL.Psi_y;


% load masks
QA = imread('./CartesianMasks/Cart64/Q1.png');
QB = imread('./CartesianMasks/Cart64/Q2.png');
QA = (QA>0.1);
QB = (QB>0.1);
% QB = ones(size(QA)); % Full sampling for side information

QA = ifftshift( QA );
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

%%
for i = [4]
% for i =1 : length(Xcell)

	fprintf('Processing image No. %d ... \n', i);
	
	Img.I1A = Xcell{i};
	Img.I1B = Ycell{i};

	
	Q1A = Q.Q1A ;
	Q1B = Q.Q1B ;
	I1A = Img.I1A;
	I1B = Img.I1B;	
	sigmaiA = sigmai.sigmaiA ;
	sigmaiB = sigmai.sigmaiB ;
	sigmaA = sigma.sigmaA ;
	sigmaB = sigma.sigmaB ;

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

	I11A=I11A./(max(max(abs(I11A))));  %image normalization

	PSNR0A=20*log10(sqrt(aa*bb)/norm(double(abs(I11A))-double(I1A),'fro')); %PSNR of normalized zero-filled reconstruction


	%-----------------for modal B-----------
	%Compute Input PSNR after adding simulated noise
	IGB=abs(ifft2(I5B));
	InputPSNRB=20*log10((sqrt(aa*bb))/norm(double(abs(IGB))-double(I1B),'fro'));
	param1.InputPSNRB=InputPSNRB;

	indexB=find(Q1B==1); %Index the sampled locations in sampling mask

	I2B=(double(I5B)).*(Q1B);  %Apply mask in DFT domain
	I11B=ifft2(I2B);          % Inverse FFT - gives zero-filled result
	I11B=I11B/(max(max(abs(I11B))));  %image normalization

	PSNR0B=20*log10(sqrt(aa*bb)/norm(double(abs(I11B))-double(I1B),'fro')); %PSNR of normalized zero-filled reconstruction

	%% for A
	FigFormatCell = { '.png', '.fig'}; %.eps
	fm = 1;
	SaveFig = 0;
	
	% show the k-space mask, samples, zero-filled image, ground true image
	figure; imshow(abs(fftshift(Q1A)))
	set(gcf, 'position', [0,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_Mask'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end

	figure;  imshow(abs(I1A))
	set(gcf, 'position', [300,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_True'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end
	
	figure; imshow(abs(I11A))
	set(gcf, 'position', [600,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_ZeroFilled'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end
	
	figure; image(abs(fftshift(I5A))); colormap gray
	set(gcf, 'position', [300,300,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_TrueSpect'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
	end
		
	figure; image(abs(fftshift(I2A))); colormap gray
	set(gcf, 'position', [600,300,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_ZeroFilledSpect'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
	end
	
%% for B
	% show the k-space mask, samples, zero-filled image, ground true image
	figure; imshow(abs(fftshift(Q1B)))
	set(gcf, 'position', [0,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_Mask_B'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end
	
	figure;  imshow(abs(I1B))
	set(gcf, 'position', [300,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_True_B'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end
	
	figure; imshow(abs(I11B))
	set(gcf, 'position', [600,100,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_ZeroFilled_B'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
	end
	
	figure; image(abs(fftshift(I5B))); colormap gray
	set(gcf, 'position', [300,300,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_TrueSpect_B'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
	end
	
	figure; image(abs(fftshift(I2B))); colormap gray
	set(gcf, 'position', [600,300,256,256])
	set(gca, 'position', [0,0,1,1])
	FigName = ['Img4_ZeroFilledSpect_B'];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
	end
	
end

%% show selected images
ImIndex = 4;
Img_Cell = [Result_cell{ImIndex}.param1.Img_Cell, Result_cell{ImIndex}.param2.Img_Cell] ;
FigFormatCell = { '.png', '.fig'}; %.eps
fm = 1;

PosStart = [140, 50] ;  lineSize = 2; lineSize2 = 1; wSize = [50, 50]; % width
BoxColor = [1,1,0]; % yellow.
SaveFig = 1;

for i = [1, 2, 3, 4, 15, 25, 35, 50]
	kp = i;
	I1A = Result_cell{ImIndex}.Img.I1A ;
	I11A = Img_Cell{kp}.Iout1 ;
	I2A=fft2(I11A);
	
	% save figures for A
% 	figure(101); imshow(abs(I11A))
% 	set(gcf, 'position', [600,100,256,256])
% 	set(gca, 'position', [0,0,1,1])
% 	pause(0.2)
% 	FigName = ['Img_Iter', num2str(kp)];
% 	if SaveFig
% 		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
% 		pause(0.2)
% 	end
	
	imX = abs(I11A) ;	
	imX = im2double(imX);
	[imWhole,imPiece] = ShowFigInFig(imX, PosStart, wSize, lineSize, BoxColor);


	figure(101); clf
	imagesc(imWhole); hold on;
	set(gcf,'Position', [0,100,size(imWhole,2), size(imWhole,1)]); axis off
	set(gca, 'Position', [0,0,1,1]);
	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2),0.4]); 
	imagesc(imPiece); axis off
	colormap(gray)
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(kp)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2);

	if SaveFig
		FigName = ['Img_Iter', num2str(kp)];
		export_fig(gcf,  [FigName, FigFormatCell{fm}], '-nocrop', '-transparent')   % export_fig(figure_handle, filename, , '-nocrop', '-transparent');		
	end
	
	

	figure(102); clf
	image(abs(fftshift(I2A))); colormap gray
	set(gcf, 'position', [600,100,256,256])
	set(gca, 'position', [0,0,1,1])
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(kp)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2)
	FigName = ['ImgFFT_Iter', num2str(kp)];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
		pause(0.2)
	end

	% residual Truth - rec with SI
	imX =  abs(I11A);
	X_GT = abs(I1A);
	imX_resid = abs(X_GT - imX);

	[imWhole,imPiece] = ShowFigInFig(imX_resid, PosStart, wSize, lineSize, BoxColor);
	imWhole(:,:,2:3) = [];	imPiece(:,:,2:3) = []; 

	cmin = 0; cmax = 0.4; % min and max color value threshold;
	figure(103); clf
	imagesc(imWhole); hold on;
	% % show wider image with colorbar
% 	set(gcf,'Position', [300,100,size(imWhole,2)*1.2, size(imWhole,1)]); axis off
% 	set(gca, 'Position', [0,0,1/1.2,1]);
% 	colormap(jet); caxis([cmin, cmax]) ;
% 	c = colorbar ; c.LineWidth = 0.5; c.FontSize = 11; c.Position = [1/1.2, 0.02, 0.05, 0.96]; % [left, bottom, width, height]
% 	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2)/1.2, 0.4]); 
	
	colormap(jet); caxis([cmin, cmax]) ;
	set(gcf,'Position', [300,100,size(imWhole,2), size(imWhole,1)]); axis off
	set(gca, 'Position', [0,0,1,1]);
	set(gcf, 'Color', 'w'); % set figure background to be white
	
	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2), 0.4]); 
	imagesc(imPiece); axis off
	colormap(jet); caxis([cmin, cmax]) ;
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(kp)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2);

	if SaveFig
		FigName = ['ImgRes_Iter', num2str(kp)];
		export_fig(gcf,  [FigName, FigFormatCell{fm}], '-nocrop', '-transparent')   % export_fig(figure_handle, filename, , '-nocrop', '-transparent');		
	end
	
	pause(1)

end

%% show the zero-filled images;
ImIndex = 4;
FigFormatCell = { '.png', '.fig'}; %.eps
fm = 1;

PosStart = [140, 50] ;  lineSize = 2; lineSize2 = 1; wSize = [50, 50]; % width
BoxColor = [1,1,0]; % yellow.
SaveFig = 1;
for i = ImIndex

	I1A = Result_cell{ImIndex}.Img.I1A ;
	I1B = Result_cell{ImIndex}.Img.I1B ;
	Q1A = Q.Q1A ;
	Q1B = Q.Q1B ;
	sigmaiA = sigmai.sigmaiA ;
	sigmaiB = sigmai.sigmaiB ;
	sigmaA = sigma.sigmaA ;
	sigmaB = sigma.sigmaB ;

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

	I11A=I11A./(max(max(abs(I11A))));  %image normalization

	PSNR0A=20*log10(sqrt(aa*bb)/norm(double(abs(I11A))-double(I1A),'fro')); %PSNR of normalized zero-filled reconstruction


	%-----------------for modal B-----------
	%Compute Input PSNR after adding simulated noise
	IGB=abs(ifft2(I5B));
	InputPSNRB=20*log10((sqrt(aa*bb))/norm(double(abs(IGB))-double(I1B),'fro'));
	param1.InputPSNRB=InputPSNRB;

	indexB=find(Q1B==1); %Index the sampled locations in sampling mask

	I2B=(double(I5B)).*(Q1B);  %Apply mask in DFT domain
	I11B=ifft2(I2B);          % Inverse FFT - gives zero-filled result
	I11B=I11B/(max(max(abs(I11B))));  %image normalization

	PSNR0B=20*log10(sqrt(aa*bb)/norm(double(abs(I11B))-double(I1B),'fro')); %PSNR of normalized zero-filled reconstruction

%% for A
	% show the k-space samples, zero-filled image, residual
% 	figure; imshow(abs(I11A))
% 	set(gcf, 'position', [0,100,256,256])
% 	set(gca, 'position', [0,0,1,1])
% 	FigName = ['Img4_ZeroFilled'];
% 	if SaveFig
% 		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
% 	end
	
	imX = abs(I11A) ;	
	imX = im2double(imX);
	[imWhole,imPiece] = ShowFigInFig(imX, PosStart, wSize, lineSize, BoxColor);

	figure(101); clf
	imagesc(imWhole); hold on;
	set(gcf,'Position', [0,100,size(imWhole,2), size(imWhole,1)]); axis off
	set(gca, 'Position', [0,0,1,1]);

	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2),0.4]); 
	imagesc(imPiece); axis off
	colormap(gray)
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(0)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2);

	if SaveFig
		FigName = ['Img_Iter', num2str(0)];
		export_fig(gcf,  [FigName, FigFormatCell{fm}], '-nocrop', '-transparent')   % export_fig(figure_handle, filename, , '-nocrop', '-transparent');		
	end
	
	
	figure(102); clf
	image(abs(fftshift(I2A))); colormap gray
	set(gcf, 'position', [600,100,256,256])
	set(gca, 'position', [0,0,1,1])
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(0)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2)
	FigName = ['ImgFFT_Iter', num2str(0)];
	if SaveFig
		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
		pause(0.2)
	end

	% residual Truth - rec with SI
	imX =  abs(I11A);
	X_GT = abs(I1A);
	imX_resid = abs(X_GT - imX);

	[imWhole,imPiece] = ShowFigInFig(imX_resid, PosStart, wSize, lineSize, BoxColor);
	imWhole(:,:,2:3) = [];	imPiece(:,:,2:3) = []; 

	cmin = 0; cmax = 0.4; % min and max color value threshold;
	figure(103); clf
	imagesc(imWhole); hold on;
	% % show wider image with colorbar
% 	set(gcf,'Position', [300,100,size(imWhole,2)*1.2, size(imWhole,1)]); axis off
% 	set(gca, 'Position', [0,0,1/1.2,1]);
% 	colormap(jet); caxis([cmin, cmax]) ;
% 	c = colorbar ; c.LineWidth = 0.5; c.FontSize = 11; c.Position = [1/1.2, 0.02, 0.05, 0.96]; % [left, bottom, width, height]
% 	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2)/1.2, 0.4]); 

	colormap(jet); caxis([cmin, cmax]) ;
	set(gcf,'Position', [300,100,size(imWhole,2), size(imWhole,1)]); axis off
	set(gca, 'Position', [0,0,1,1]);
	set(gcf, 'Color', 'w'); % set figure background to be white
	axes('Position',[0,0,0.4*size(imWhole,1)/size(imWhole,2), 0.4]); 

	imagesc(imPiece); axis off
	colormap(jet); caxis([cmin, cmax]) ;
	annotation(gcf,'textbox', [0.8 0 0.1 0.1], 'Color',[1 1 0],	'String',{num2str(0)}, 'LineStyle','none', 'FontSize',16,	'FitBoxToText','off');
	pause(0.2);

	if SaveFig
		FigName = ['ImgRes_Iter', num2str(0)];
		export_fig(gcf,  [FigName, FigFormatCell{fm}], '-nocrop', '-transparent')   % export_fig(figure_handle, filename, , '-nocrop', '-transparent');		
	end
	
	pause(1)
	
% %% for B
% 	% show the k-space samples, zero-filled image, residual
% 	figure; imshow(abs(I11B))
% 	set(gcf, 'position', [0,100,256,256])
% 	set(gca, 'position', [0,0,1,1])
% 	FigName = ['Img4_ZeroFilled_B'];
% 	if SaveFig
% 		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);
% 	end
% 	
% 	figure; image(abs(fftshift(I2B))); colormap gray
% 	set(gcf, 'position', [600,100,256,256])
% 	set(gca, 'position', [0,0,1,1])
% 	FigName = ['Img4_ZeroFilledSpect_B'];
% 	if SaveFig
% 		export_fig(gcf, [FigName, FigFormatCell{fm}], '-nocrop', '-transparent') % export_fig(figure_handle, filename);	
% 	end
end













