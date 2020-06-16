
function [imWhole,imPiece,varargout] = ShowFigInFig(im, PosStart, wSize, lineSize, BoxColor, varargin)
% Show figure in a figure. First, draw a box in the image and 
% also copy this piece of image in the box, at last, 
% show the cropped piece in the image,too.

%Input：
% im:	original image, RGB or gray.
% PosStart： Starting position for drawing the box, i.e., the coordinate of the top left point.
% wSize：   [wx, wy], the size of the box.
% lineSize： the width of the box line.
% BoxColor：    [r,  g,  b] , the color of the box
%----------------------------------------------------------------------
%Output：
% imWhole： the whole image with a box in it.
% imPiece: the piece of image in that box.

% Example
% [imWhole,imPiece] = ShowFigInFig(im, PosStart, wSize, lineSize,
% BoxColor);

% 23/11/2016 Code finished.

imWhole = im;
[m, n, z] = size(imWhole);

RecArea = [PosStart, wSize]; % the area in the box

imPiece = imcrop(imWhole,RecArea); % cropped image piece

PosStart_Piece = [1, 1]; % Starting position for drawing the box in the image piece.

wSize_Piece = size(imPiece(:,:,1))-[1, 1]; % size of the box.


imWhole = drawBox(imWhole,PosStart,wSize,lineSize, BoxColor );

imPiece = drawBox(imPiece,PosStart_Piece,wSize_Piece,1, BoxColor );


%%
show_figure = 0;
save_figure = 0;

if show_figure
	figure;
	imshow(imWhole); hold on;
	set(gcf,'Position', [100,100,ImWidth,ImHeigth]);
	set(gca, 'Position', [0,0,1,1]);

	axes('Position',[0,0,0.4,0.4]);
	imshow(imPiece)
	pause(0.2);
end


if save_figure
	addpath('../export_fig');
	FigFormat = ['.png'] ; % ['.jpg'] ; % ['.eps'];
% 	current_date = date;
	DATE = datestr(now,30); % 30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517 
	FileName = ['Im_' , DATE];
	export_fig(gcf, [FileName, FigFormat])  % export_fig(figure_handle, filename);		
end

end
		

function [ dest ] = drawBox( im, PosStart, wSize,  lineSize, color )
%简介：
% %将图像画上有颜色的框图，如果输入是灰度图，先转换为彩色图像，再画框图
% 图像矩阵
% 行向量方向  是  y
% 列向量方向  是  x
%----------------------------------------------------------------------
%输入：
% im：        原始图像，可以为灰度图，可为彩色图
% PosStart：         左上角坐标   [x1, y1]
% wSize：   框的大小      [wx, wy]
% lineSize： 线的宽度
% color：     线的颜色      [r,  g,  b] 
%----------------------------------------------------------------------
%输出：
% dest：           画好了的图像
%----------------------------------------------------------------------

%flag=1: 有缺口的框
%flag=2: 无缺口的框
flag = 2;


%判断输入参数个数
if nargin < 5
    color = [255 255 0];
end

if nargin < 4
    lineSize = 1;
end

if nargin < 3
    disp('输入参数不够 !!!');
    return;
end





%判断框的边界问题
[yA, xA, z] = size(im);
x1 = PosStart(1);
y1 = PosStart(2);
wx = wSize(1);
wy = wSize(2);
if  x1>xA || ...
        y1>yA||...
        (x1+wx)>xA||...
        (y1+wy)>yA

    disp('画的框将超过图像 !!!');
    return;
end

%如果是单通道的灰度图，转成3通道的图像
if 1==z
    dest(:, : ,1) = im;
    dest(:, : ,2) = im;
    dest(:, : ,3) = im;
else
    dest = im;
end

%开始画框图
for c = 1 : 3                 %3个通道，r，g，b分别画
    for dl = 1 : lineSize   %线的宽度，线条是向外面扩展的
        d = dl - 1;
        if  1==flag %有缺口的框
            dest(  y1-d ,            x1:(x1+wx) ,  c  ) =  color(c); %上方线条
            dest(  y1+wy+d ,     x1:(x1+wx) , c  ) =  color(c); %下方线条
            dest(  y1:(y1+wy) ,   x1-d ,           c  ) =  color(c); %左方线条
            dest(  y1:(y1+wy) ,   x1+wx+d ,    c  ) =  color(c); %左方线条
        elseif 2==flag %无缺口的框
            dest(  y1-d ,            (x1-d):(x1+wx+d) ,  c  ) =  color(c); %上方线条
            dest(  y1+wy+d ,    (x1-d):(x1+wx+d) ,  c  ) =  color(c); %下方线条
            dest(  (y1-d):(y1+wy+d) ,   x1-d ,           c  ) =  color(c); %左方线条
            dest(  (y1-d):(y1+wy+d) ,   x1+wx+d ,    c  ) =  color(c); %左方线条
        end
    end    
end %主循环尾


end %函数尾