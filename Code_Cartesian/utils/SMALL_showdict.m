function x = SMALL_showdict(D,sz,n,m,varargin)
%% SMALL_SHOWDICT Display a dictionary of image patches.
%  Reimplementation of showdict function from KSVD toolbox with Image
%  Processing toolbox dependecies removed
%
%  SMALL_SHOWDICT(D,SZ,N,M) displays the contents of the dictionary D, whos
%  columns are 2-D image patches (in column-major order). SZ = [SX SY] is
%  the size of the image patches. SHOWDICT displays the atoms on an N x M
%  grid. If there are more atoms in D then only the first N*M are
%  displayed.
%
%  SMALL_SHOWDICT(...,'lines') separates the dictionary atoms by black lines.
%  SMALL_SHOWDICT(...,'whitelines') separates the dictionary atoms by white
%  lines.
%
%  SMALL_SHOWDICT(...,'linewidth',W) when used with either 'lines' or
%  'whitelines' sets the width of the lines to W pixels (default=1).
%
%  SMALL_SHOWDICT(...,'highcontrast') increases the contrast of the figure by
%  normalizing the intensity values of each atom individually to the range
%  of [0,1] (the default behavior is to normalize the values of the entire
%  figure to [0,1] as one image). Note that in this way, the relative
%  intensities of the atoms are not maintained.
%
%  X = SMALL_SHOWDICT(...) returns a bitmat of the dictionary image without
%  displaying the figure.

%   Centre for Digital Music, Queen Mary, University of London.
%   This file copyright 2011 Ivan Damnjanovic.
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation; either version 2 of the
%   License, or (at your option) any later version.  See the file
%   COPYING included with this distribution for more information.
%   

% This function is from SMALLbox Version 2.1

% 25th May, 2012
% 
% Version 2.0 of SMALLbox is the release candidate that is distributed 
% for testing and bug fixing purposes. Please send all bugs, requests and suggestions to:
% 
% luis.figueira@soundsoftware.ac.uk
% ---------------------------------------------------------------------------
% 
% Copyright (2012): 	Luis Figueira, Ivan Damnjanovic, Matthew Davies
% 			Centre for Digital Music, 
% 			Queen Mary University of London
% 
% SMALLbox is distributed under the terms of the GNU General Public License 3
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% ---------------------------------------------------------------------------
% For more information on the SMALL Project, please visit the following website:
% 
% http://small-project.eu

if (size(D,2) < n*m)
  D = [D zeros(size(D,1),n*m-size(D,2))];
end


%%%  parse input arguments  %%%

linewidth = 1;
highcontrast = 0;
drawlines = 0;
linecolor = 0;

for i = 1:length(varargin)
  if (~ischar(varargin{i}))
    continue;
  end
  switch(varargin{i})
    case 'highcontrast'
      highcontrast = 1;
    case 'lines'
      drawlines = 1;
    case 'whitelines'
      drawlines = 1;
      linecolor = 1;
    case 'linewidth'
      linewidth = varargin{i+1};
  end
end



%%%  create dictionary image  %%%


if (drawlines)
  
  D = [D ; nan(sz(1)*linewidth,size(D,2))];
  sz(2) = sz(2)+linewidth;
  x = col2imstep(D(:,1:n*m),[n m].*sz, sz, sz);
  sz = [sz(2) sz(1)];
  D = im2colstep(x',sz, sz);
  D = [D ; nan(sz(1)*linewidth,size(D,2))];
  sz(2) = sz(2)+linewidth;
  x = col2imstep(D(:,1:n*m),[m n].*sz,sz,sz);
  x = x';
  x = x(1:end-linewidth,1:end-linewidth);
  
  if (highcontrast)
    for i = 0:n-1
      for j = 0:m-1
        x(i*sz(1)+1:i*sz(1)+sz(1)-linewidth, j*sz(2)+1:j*sz(2)+sz(2)-linewidth) = ...
          imnormalize(x(i*sz(1)+1:i*sz(1)+sz(1)-linewidth, j*sz(2)+1:j*sz(2)+sz(2)-linewidth));
      end
    end
  else
    x = imnormalize(x);
  end
  
  x(isnan(x)) = linecolor;
  
else
  
  x = col2imstep(D(:,1:n*m),[n m].*sz, sz, sz);
  
  if (highcontrast)
    for i = 0:n-1
      for j = 0:m-1
        x(i*sz(1)+1:i*sz(1)+sz(1), j*sz(2)+1:j*sz(2)+sz(2)) = ...
          imnormalize(x(i*sz(1)+1:i*sz(1)+sz(1), j*sz(2)+1:j*sz(2)+sz(2)));
      end
    end
  else
    x = imnormalize(x);
  end
end


if (nargout==0)
    imagesc(dictimg);colormap(gray);axis off; axis image; 
end
