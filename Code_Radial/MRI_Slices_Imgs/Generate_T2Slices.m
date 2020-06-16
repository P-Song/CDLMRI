%% transform T2 slices. From 512x512 to 256x256
% input T2Slices_SZ512.mat, 
% output T1Slices_SZ512.mat

load('T2Slices_SZ512.mat')

pos = [50, 490, 80, 430];
NUMROWS = 256; NUMCOLS = 256; 
T2Slices = single([]);
for i =1 : size(data,3)
	im = flipud(data(:,:, i)') ;
	im = single(im(pos(1): pos(2),  pos(3) : pos(4)));
    im = im./max(abs(im(:)));
	max(im(:)), min(im(:)),
	T2Slices(:,:,i) = imresize( im , [NUMROWS, NUMCOLS], 'bicubic'); 
end

figure;
for i =1 : 12
	imagesc(T2Slices(:,:, i))
	colormap gray
	set(gcf,'Position', [100,100, 256, 256]); axis off
	set(gca, 'Position', [0,0,1,1]);
	set(gcf, 'Color', 'w'); % set figure background to be white
	pause()
end


% figure;
% for i =1 : 12
% % 	subplot(3,4,i)
% 	imagesc(data(:,:, i))
% 	colormap gray
% 	set(gcf,'Position', [100,100, 256, 256]); axis off
% 	set(gca, 'Position', [0,0,1,1]);
% 	set(gcf, 'Color', 'w'); % set figure background to be white
% 	pause(0.2)
% end


