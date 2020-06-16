function imgs = load_True_images(paths)
% 23/05/2016 load original true images.
imgs = cell(size(paths));
for i = 1:numel(paths)
    X = imread(paths{i});
	X = im2double(X); % to reduce memory usage
    imgs{i} = X;
end