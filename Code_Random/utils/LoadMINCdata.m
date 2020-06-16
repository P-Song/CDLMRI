
%%
addpath('../export_fig')
	% Load training images from MINC format data
	directoryX = '../TrainImages_MRI'; 
	directoryY = '../TrainImages_MRI'; 
	patternX = 't1_icbm_normal_1mm_pn0_rf0.mnc';
	patternY = 't2_icbm_normal_1mm_pn0_rf0.mnc';
	
	XpathCell = glob(directoryX, patternX );
	[T1,scaninfo] = loadminc(XpathCell{1});
	YpathCell = glob(directoryY, patternY );
	[T2,scaninfo] = loadminc(YpathCell{1});
	
	Channels = size(T1,3);
	SliceNo =[1:10: 141] ;
	
	h1 = figure;
	h2 = figure;
	FigFormat = ['.png'];
	SaveFig = 1;
	for i = SliceNo
		figure(h1)
		im = T1(:,:,i) ;
		amax = max(im(:)); amin = min(im(:));
		im = (im - amin)./(amax - amin) ; % normalization to 0 ~ 1;
		Xcell{i} = im;
		imshow(im); axis off; colormap gray; 
		set(h1,'position', [100,100, size(im,2), size(im,1)])
		set(gca, 'position', [0,0, 1, 1])
		pause(0.1);
		if SaveFig
			FileName = ['BrainWeb_T1_Slice', num2str(i) ];
			export_fig(h1, [FileName, FigFormat])  
		end	
		
		
		figure(h2)
		im = T2(:,:,i) ;
		amax = max(im(:)); amin = min(im(:));
		im = (im - amin)./(amax - amin) ; % normalization to 0 ~ 1;
		Ycell{i} = im ;
		imshow(im); axis off; colormap gray; 
		set(h2,'position', [400,100, size(im,2), size(im,1)])
		set(gca, 'position', [0,0, 1, 1])
		pause(0.1);
		if SaveFig
			FileName = ['BrainWeb_T2_Slice', num2str(i) ];
			export_fig(h2, [FileName, FigFormat])  
		end	
	
	end
	
	disp('done!')