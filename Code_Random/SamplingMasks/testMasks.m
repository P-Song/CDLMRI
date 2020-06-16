% compare different masks with different downsampling folds.

clear;

% QA = load('./10fold/Q1.mat');
% QB = load('./10fold/Q2.mat');

% QA = load('./20fold/Q1.mat');
% QB = load('./20fold/Q2.mat');
% 
% QA = load('./30fold/Q1.mat');
% QB = load('./30fold/Q2.mat');

% QA = load('./40fold/Q1.mat');
% QB = load('./40fold/Q2.mat');


QA = load('./40fold_Size128/Q1.mat');
QB = load('./40fold_Size128/Q2.mat');
% 

QA = ( QA.sampling_mask );
QB = ( QB.sampling_mask );
% QA = ifftshift( QA.sampling_mask );
% QB = ifftshift( QB.sampling_mask );

indexQA=find(QA==1); %Index the sampled locations in sampling mask
indexQB=find(QB==1); %Index the sampled locations in sampling mask


[m , n] = size(QA);
foldQA = (m * n) / length(indexQA)
foldQB = (m * n) / length(indexQB)

figure;
imshow(QA);
figure;
imshow(QB);

