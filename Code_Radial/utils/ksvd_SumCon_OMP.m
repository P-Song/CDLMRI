function [D,Gamma,err,gerr,varargout] = ksvd(params,varargin)
% 31/12/2016 Perform cleardict (clear dictionary) after computing
% RMSE.XY,etc. Otherwise, the RMSE.XY will fluctuate violently.
% 18/10/2016 Store training error of common dicts and unique dicts together.
% 30/09/2016 Store all the information in variable 'info_all' and pass out.
% 09/04/2016 Use OMP for sparse coding with the sum of sc, sx, sy as the sparsity constraint. 
% 16/11/2015 Use OMP for sparse coding. Use local dictionary update mode, that is, 
% divide the whole dictionary into 2 parts which only one of them is updated after the 
% sparse-coding while the other part is kept fixed, and then after a few iterations, 
% fix the updated part and update the other one after sparse-coding. 

% This is a coupled dictionary learning algorithm adapted from K-SVD to solve following problem.
% Given two training data matrices X: N x T and Y: M x T (same number of columns),
% finds a decomposition
%
%      X = Psi_cx*Z_c + Psi_x*Z_x
%      Y = Psi_cy*Z_c + Psi_y*Z_y
%
% where all the dictionaries Psi_cx, Psi, Psi_cy, and Psi_y have dimensions 
% N x sub_dictsize. The matrices Z_c, Z_x, and Z_y are sparse and contain the coefficients
% that represent X and Y in the corresponding dictionaries. The above
% decomposition is computed by solving
%
%   minimize 0.5*||X - Psi_cx*Z_c - Psi_x*Z_x||_F^2 + ||Y - Psi_cy*Z_c - Psi_y*Z_y||_F^2
%                 + lambda1*||Z_c||_1 + lambda2*||Z_x||_1 + lambda3*||Z_y||_1
%	where lambda1=lambda2=lambda3=1


%  [D,GAMMA] = KSVD(PARAMS) runs the K-SVD dictionary training algorithm on
%  the specified set of signals, returning the trained dictionary D with specified structure and the
%  signal representation matrix GAMMA, where
% D = [Psi_cx, Psi_x, zeros(N,sub_dictsize) ; ...
%		  Psi_cy, zeros(M,sub_dictsize), Psi_y];
% Gamma = [Z_c; Z_x; Z_y];
% params.data = [X;Y];
%

% In sparse-coding stage, modified OMP uses sparsity-based minimization for sparse-coding and the
% optimization problem is given by
%
%      min  ||data-D*Gamma||_F^2     
%    Gamma
%	   s.t.  ||Gamma_i(ind1)||_0 <= s_c, (i.e. ||Z_c||_0 <= s_c);
%			  ||Gamma_i(ind2)||_0 <= s_x, (i.e. ||Z_x||_0 <= s_x);
%			  ||Gamma_i(ind3)||_0 <= s_y, (i.e. ||Z_y||_0 <= s_y);
%  where data is the set of training signals, Gamma_i is the i-th column of
%  Gamma, and s_c, s_x and s_y are the target sparsities of Z_c, Z_x and Z_y, respectively.

% In dictionary update stage, D is decoupled into 3 independent sub
% dictionaries D1 = [Psi_cx; Psi_cy]; D2 = [Psi_x; zeros]; D3 = [zeros;Psi_y].
% When the updating atom belongs to D1, to classic k-svd. When the updating
% atom belongs to D2, only update the first half corresponding to Psi_x.
% When belongs to D3, only update the last half corresponding to Psi_y.

%
%   Summary of all fields in PARAMS:
%   --------------------------------
%     'data'                   training data
%     'Tdata' / 'Edata'        sparse-coding target
%     'initdict' / 'dictsize'  initial dictionary / dictionary size
%	'ColIndex' - a vector which includes column size of each dictionary.
%   'RowIndex' - a vector which includes row size of each dictionary.
% 	'Psi_cx0', 'Psi_x0', 'Psi_cy0', 'Psi_y0' - True dictionaries
%    'codemode' - Sparse-coding target mode. should be string 'sparsity'.

% This code has been moditied by :
% Pingfan Song
% Department of Electronic and Electrical Engineering
% University College London, Gower Street, UK
% pingfan.song.14@ucl.ac.uk
% 13/11/2015

%%
%KSVD K-SVD dictionary training.
%  [D,GAMMA] = KSVD(PARAMS) runs the K-SVD dictionary training algorithm on
%  the specified set of signals, returning the trained dictionary D and the
%  signal representation matrix GAMMA.
%
%  KSVD has two modes of operation: sparsity-based and error-based. For
%  sparsity-based minimization, the optimization problem is given by
%
%      min  |X-D*GAMMA|_F^2      s.t.  |Gamma_i|_0 <= T
%    D,Gamma
%
%  where X is the set of training signals, Gamma_i is the i-th column of
%  Gamma, and T is the target sparsity. For error-based minimization, the
%  optimization problem is given by
%
%      min  |Gamma|_0      s.t.  |X_i - D*Gamma_i|_2 <= EPSILON
%    D,Gamma
%
%  where X_i is the i-th training signal, and EPSILON is the target error.
%
%  [D,GAMMA,ERR] = KSVD(PARAMS) also returns the target function values
%  after each algorithm iteration. For sparsity-constrained minimization,
%  the returned values are given by
%
%    ERR(D,GAMMA) = RMSE(X,D*GAMMA) = sqrt( |X-D*GAMMA|_F^2 / numel(X) ) .
%
%  For error-constrained minimization, the returned values are given by
%
%    ERR(D,GAMMA) = mean{ |Gamma_i|_0 } = |Gamma|_0 / size(X,2) .
%
%  Error computation slightly increases function runtime.
%
%  [D,GAMMA,ERR,GERR] = KSVD(PARAMS) computes the target function values on
%  the specified set of test signals as well, usually for the purpose of
%  validation (testing the generalization of the dictionary). This requires
%  that the field 'testdata' be present in PARAMS (see below). The length
%  of ERR and GERR is identical.
%
%  [...] = KSVD(...,VERBOSE) where VERBOSE is a character string, specifies
%  messages to be printed during the training iterations. VERBOSE should
%  contain one or more of the characters 'i', 'r' and 't', each of which
%  corresponds to a certain piece of information:
%           
%    i - iteration number
%    r - number of replaced atoms
%    t - target function value (and its value on the test data if provided)
%
%  Specifying either 'r', 't' or both, also implies 'i' automatically. For
%  example, KSVD(PARAMS,'tr') prints the iteration number, number of
%  replaced atoms, and target function value, at the end of each iteration.
%  The default value for VERBOSE is 't'. Specifying VERBOSE='' invokes
%  silent mode, and cancels all messages.
%
%  [...] = KSVD(...,MSGDELTA) specifies additional messages to be printed
%  within each iteration. MSGDELTA should be a positive number representing
%  the interval in seconds between messages. A zero or negative value
%  indicates no such messages (default). Note that specifying VERBOSE=''
%  causes KSVD to run in silent mode, ignoring the value of MSGDELTA.
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'data' - Training data.
%      A matrix containing the training signals as its columns.
%
%    'Tdata' / 'Edata' - Sparse coding target.
%      Specifies the number of coefficients (Tdata) or the target error in
%      L2-norm (Edata) for coding each signal. If only one is present, that
%      value is used. If both are present, Tdata is used, unless the field
%      'codemode' is specified (below).
%
%    'initdict' / 'dictsize' - Initial dictionary / no. of atoms to train.
%      At least one of these two should be present in PARAMS.
%
%      'dictsize' specifies the number of dictionary atoms to train. If it
%      is specified without the parameter 'initdict', the dictionary is
%      initialized with dictsize randomly selected training signals.
%
%      'initdict' specifies the initial dictionary for the training. It
%      should be either a matrix of size NxL, where N=size(data,1), or an
%      index vector of length L, specifying the indices of the examples to
%      use as initial atoms. If 'dictsize' and 'initdict' are both present,
%      L must be >= dictsize, and in this case the dictionary is
%      initialized using the first dictsize columns from initdict. If only
%      'initdict' is specified, dictsize is set to L.
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'testdata' - Validation data.
%      If present, specifies data on which to compute generalization error.
%      Should be a matrix containing the validation signals as its columns.
%
%    'iternum' - Number of training iterations.
%      Specifies the number of K-SVD iterations to perform. If not
%      specified, the default is 10.
%
%    'memusage' - Memory usage.
%      This parameter controls memory usage of the function. 'memusage'
%      should be one of the strings 'high', 'normal' (default) or 'low'.
%      When set to 'high', the fastest implementation of OMP is used, which
%      involves precomputing both G=D'*D and DtX=D'*X. This increasese
%      speed but also requires a significant amount of memory. When set to
%      'normal', only the matrix G is precomputed, which requires much less
%      memory but slightly decreases performance. Finally, when set to
%      'low', neither matrix is precomputed. This should only be used when
%      the trained dictionary is highly redundant and memory resources are
%      very low, as this will dramatically increase runtime. See function
%      OMP for more details.
%
%    'codemode' - Sparse-coding target mode.
%      Specifies whether the 'Tdata' or 'Edata' fields should be used for
%      the sparse-coding stopping criterion. This is useful when both
%      fields are present in PARAMS. 'codemode' should be one of the
%      strings 'sparsity' or 'error'. If it is not present, and both fields
%      are specified, sparsity-based coding takes place.
%
%    'exact' - Exact K-SVD update.
%      Specifies whether the exact or approximate dictionary update
%      should be used. By default, the approximate computation is used,
%      which is significantly faster and requires less memory. Specifying a
%      nonzero value for 'exact' causes the exact computation to be used
%      instead, which slows down the method but provides slightly improved
%      results. The exact update uses SVD to solve the rank-1 minimization
%      problem, while the approximate upate performs alternate-optimization
%      to solve this problem.
%
%
%  Optional fields in PARAMS - advanced:
%  -------------------------------------
%
%    'maxatoms' - Maximal number of atoms in signal representation.
%      When error-based sparse coding is used, this parameter can be used
%      to specify a hard limit on the number of atoms in each signal
%      representation (see parameter 'maxatoms' in OMP2 for more details).
%
%    'muthresh' - Mutual incoherence threshold.
%      This parameter can be used to control the mutual incoherence of the
%      trained dictionary, and is typically between 0.9 and 1. At the end
%      of each iteration, the trained dictionary is "cleaned" by discarding
%      atoms with correlation > muthresh. The default value for muthresh is
%      0.99. Specifying a value of 1 or higher cancels this type of
%      cleaning completely. Note: the trained dictionary is not guaranteed
%      to have a mutual incoherence less than muthresh. However, a method
%      to track this is using the VERBOSE parameter to print the number of
%      replaced atoms each iteration; when this number drops near zero, it
%      is more likely that the mutual incoherence of the dictionary is
%      below muthresh.
%
%
%   Summary of all fields in PARAMS:
%   --------------------------------
%
%   Required:
%     'data'                   training data
%     'Tdata' / 'Edata'        sparse-coding target
%     'initdict' / 'dictsize'  initial dictionary / dictionary size
%
%   Optional (default values in parentheses):
%     'testdata'               validation data (none)
%     'iternum'                number of training iterations (10)
%     'memusage'               'low, 'normal' or 'high' ('normal')
%     'codemode'               'sparsity' or 'error' ('sparsity')
%     'exact'                  exact update instead of approximate (0)
%     'maxatoms'               max # of atoms in error sparse-coding (none)
%     'muthresh'               mutual incoherence threshold (0.99)
%
%
%  References:
%  [1] M. Aharon, M. Elad, and A.M. Bruckstein, "The K-SVD: An Algorithm
%      for Designing of Overcomplete Dictionaries for Sparse
%      Representation", the IEEE Trans. On Signal Processing, Vol. 54, no.
%      11, pp. 4311-4322, November 2006.
%  [2] M. Elad, R. Rubinstein, and M. Zibulevsky, "Efficient Implementation
%      of the K-SVD Algorithm using Batch Orthogonal Matching Pursuit",
%      Technical Report - CS, Technion, April 2008.
%
%  See also KSVDDENOISE, OMPDENOISE, OMP, OMP2.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009

% This code has been moditied by :
% Pingfan Song
% Department of Electronic and Electrical Engineering
% University College London, Gower Street, UK
% pingfan.song.14@ucl.ac.uk
% 13/11/2015


global CODE_SPARSITY CODE_ERROR codemode
global MEM_LOW MEM_NORMAL MEM_HIGH memusage
global ompfunc ompparams exactsvd
global ind1 ind2 ind3 ind1_row ind2_row

CODE_SPARSITY = 1;
CODE_ERROR = 2;

MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;


%%%%% parse input parameters %%%%%


data = params.data;

ompparams = {'checkdict','off'};

% coding mode %

if (isfield(params,'codemode'))
  switch lower(params.codemode)
    case 'sparsity'
      codemode = CODE_SPARSITY;
      thresh = params.Tdata;
    case 'error'
      codemode = CODE_ERROR;
      thresh = params.Edata;
    otherwise
      error('Invalid coding mode specified');
  end
elseif (isfield(params,'Tdata'))
  codemode = CODE_SPARSITY;
  thresh = params.Tdata;
elseif (isfield(params,'Edata'))
  codemode = CODE_ERROR;
  thresh = params.Edata;

else
  error('Data sparse-coding target not specified');
end


% max number of atoms %

if (codemode==CODE_ERROR && isfield(params,'maxatoms'))
  ompparams{end+1} = 'maxatoms';
  ompparams{end+1} = params.maxatoms;
end


% memory usage %

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% iteration count %

if (isfield(params,'iternum'))
  iternum = params.iternum;
else
  iternum = 10;
end


% omp function %

if (codemode == CODE_SPARSITY)
  ompfunc = @omp;
else
  ompfunc = @omp2;
end


% status messages %

printiter = 0;
printreplaced = 0;
printerr = 0;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'
      printiter = 1;
    case 'r'
      printiter = 1;
      printreplaced = 1;
    case 't'
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1; 
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% data norms %

XtX = []; XtXg = [];
if (codemode==CODE_ERROR && memusage==MEM_HIGH)
  XtX = colnorms_squared(data);
  if (testgen)
    XtXg = colnorms_squared(testdata);
  end
end


% mutual incoherence limit %

if (isfield(params,'muthresh'))
  muthresh = params.muthresh;
else
  muthresh = 0.99;
end
if (muthresh < 0)
  error('invalid muthresh value, must be non-negative');
end


% exact svd computation %

exactsvd = 0;
if (isfield(params,'exact') && params.exact~=0)
  exactsvd = 1;
end


% determine dictionary size %

if (isfield(params,'initdict'))
  if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
    dictsize = length(params.initdict);
  else
    dictsize = size(params.initdict,2);
  end
end
if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
  dictsize = params.dictsize;
end

if (size(data,2) < dictsize)
  error('Number of training signals is smaller than number of atoms to train');
end


% initialize the dictionary %

if (isfield(params,'initdict'))
  if (any(size(params.initdict)==1) && all(iswhole(params.initdict(:))))
    D = data(:,params.initdict(1:dictsize));
  else
    if (size(params.initdict,1)~=size(data,1) || size(params.initdict,2)<dictsize)
      error('Invalid initial dictionary');
    end
    D = params.initdict(:,1:dictsize);
  end
else
  data_ids = find(colnorms_squared(data) > 1e-6);   % ensure no zero data elements are chosen
  perm = randperm(length(data_ids));
  D = data(:,data_ids(perm(1:dictsize)));
end


% normalize the dictionary %

D = normcols(D);

err = zeros(1,iternum);
gerr = zeros(1,iternum);

if (codemode == CODE_SPARSITY)
  errstr = 'RMSE';
else
  errstr = 'mean atomnum';
end



%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%

% ColIndex is a vector which includes column size of each dictionary.
% RowIndex is a vector which includes row size of each dictionary.
ColIndex = params.ColIndex;
RowIndex = params.RowIndex;

ind1 = 1:ColIndex(1);
ind2 = (ColIndex(1)+1) : (ColIndex(1)+ColIndex(2));
ind3 = (ColIndex(1)+ColIndex(2)+1) : (ColIndex(1)+ColIndex(2)+ColIndex(3));

ind1_row = 1: RowIndex(1);
ind2_row = (RowIndex(1)+1): (RowIndex(1)+RowIndex(2));

% Check whether there exist true dictionaries
TrueDict_flag = 0 ;
if (isfield(params,'TrueDict'))
	TrueDict_flag = 1;
	TrueDict = params.TrueDict ;
	
	Psi_cx0 = TrueDict(ind1_row,ind1);
	Psi_x0 = TrueDict(ind1_row,ind2);
	Psi_cy0 = TrueDict(ind2_row,ind1);
	Psi_y0 = TrueDict(ind2_row,ind3);
end

TrueGamma_flag = 0 ;
if (isfield(params,'TrueGamma'))
	TrueGamma_flag = 1 ;
	Z_c0 = params.TrueGamma(ind1,:) ;
	Z_x0 = params.TrueGamma(ind2,:) ;
	Z_y0 = params.TrueGamma(ind3,:) ;
end

SparsityArray_flag = 0 ;
if (isfield(params,'sparsity'))
	SparsityArray_flag = 1 ;
% sparsity
	s_c = params.sparsity(1);
	s_x = params.sparsity(2); 
	s_y = params.sparsity(3);
end


% training dataset size
NumTrain = size(data,2);

Gamma = zeros(dictsize, NumTrain);

Main_Iternum = params.Main_Iternum;
Sub_Iternum = params.Sub_Iternum;

% initial
% dist_cx = nan ;
% dist_cy = nan ;
% dist_x = nan ; 
% dist_y = nan ; 
% 
% ratio_cx = 0;
% ratio_cy = 0;
% ratio_x = 0;
% ratio_y = 0;

dist = [];
ratio = [];

% dist_mc = [];
% ratio_mc = [];
% Conv_mc = [];
err_XY = [] ;

info_all = sprintf('Coupled Dictionary Learning ... \n');
iter_all = 0;

for main_iter = 1:Main_Iternum
	
	disp('---------------------------------------------------');
	disp(['cycle = ', num2str(main_iter), '/', num2str(Main_Iternum)]);	
	
	% store the prompt info.
	info_prompt = sprintf('---------------------------------------------------\n');
	info_prompt = sprintf(['%s','cycle = ', num2str(main_iter), '/', num2str(Main_Iternum)], info_prompt);
	info_all = sprintf('%s\n%s ', info_all, info_prompt);
	
	for dict_index = 1 : 2;
		info_prompt2 = [];
		if dict_index == 1;
			disp('Train common dictionaries.')
			info_prompt2 = sprintf(['%s','Train common dictionaries.'], info_prompt2);
		else
			disp('Train unique dictionaries.')
			info_prompt2 = sprintf(['%s','Train unique dictionaries.'], info_prompt2);
		end
		info_all = sprintf('%s\n%s ', info_all, info_prompt2);

	for iter = 1:Sub_Iternum
		
		iter_all = iter_all + 1 ;
  
		row_n = (main_iter-1)*(Sub_Iternum) + iter; % indicator of iterations;

		G = [];
		if (memusage >= MEM_NORMAL)
			G = D'*D;
		end


		  %%%%%  sparse coding  %%%%%

		Gamma = sparsecode(data,D,XtX,G,thresh); 
		Gam = full(Gamma);
	% 	Gamsum(1,:) = sort( sum(abs(Gam(ind1,:))>1e-7 , 1), 'descend');
	% 	Gamsum(2,:) = sort( sum(abs(Gam(ind2,:))>1e-7 , 1), 'descend' );
	% 	Gamsum(3,:) = sort( sum(abs(Gam(ind3,:))>1e-7 , 1), 'descend' );
		mean_s = zeros(3,1);
		mean_s(1) = sum(sum(abs(Gam(ind1,:))>1e-7 , 1))/ NumTrain;
		mean_s(2) = sum(sum(abs(Gam(ind2,:))>1e-7 , 1))/ NumTrain;
		mean_s(3) = sum(sum(abs(Gam(ind3,:))>1e-7 , 1))/ NumTrain;

		
		%%%%%  dictionary update  %%%%%

		replaced_atoms = zeros(1,dictsize);  % mark each atom replaced by optimize_atom

		unused_sigs = 1:size(data,2);  % tracks the signals that were used to replace "dead" atoms.
									 % makes sure the same signal is not selected twice

		p = randperm(dictsize);
		tid = timerinit('updating atoms', dictsize);

		for j = 1:dictsize
% 			if iter <= Sub_Iternum
			if dict_index == 1
				if sum(p(j)==ind1)	% if this combined atom is composed of an atom in dictionary Psi_cx and an atom in dictionary Psi_cy, do svd to update the whole combined atom.
					[D(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data,D,p(j),Gamma,unused_sigs,replaced_atoms);
					Gamma(p(j),data_indices) = gamma_j;
				end
			else
				if sum(p(j)==ind2)	% if this combined atom is composed of a atom in dictionary Psi_x and zeros, do svd to update only the part that belongs to Psi_x.
					data_ind2 = data(ind1_row,:); 
					D_ind2 = D(ind1_row,:);
					[D_ind2(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data_ind2,D_ind2,p(j),Gamma,unused_sigs,replaced_atoms);
					Gamma(p(j),data_indices) = gamma_j;
					D(ind1_row,:) = D_ind2;
				elseif  sum(p(j)==ind3) % tha is sum(p(j)==ind3)>0 % if this combined atom is composed of a atom in dictionary Psi_x and zeros, do svd to update only the part that belongs to Psi_x.
					data_ind3 = data(ind2_row,:); 
					D_ind3 = D(ind2_row,:);
					[D_ind3(:,p(j)),gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(data_ind3,D_ind3,p(j),Gamma,unused_sigs,replaced_atoms);
					Gamma(p(j),data_indices) = gamma_j;
					D(ind2_row,:) = D_ind3;
				end
			end

			if (msgdelta>0)
				timereta(tid, j, msgdelta);
			end
		
		end

		if (msgdelta>0)
			printf('updating atoms: iteration %d/%d', dictsize, dictsize);
		end


	  %%%%%  compute error  %%%%%

		if (comperr)
			err(iter) = compute_err(D,Gamma,data);
		end
		
		if (testgen)
			if (memusage >= MEM_NORMAL)
				G = D'*D;
			end
			GammaG = sparsecode(testdata,D,XtXg,G,thresh);
			gerr(iter) = compute_err(D,GammaG,testdata);
		end
   
  
   %%=================compute atoms recovery ratio======================

		Psi_cx = D(ind1_row,ind1);
		Psi_x = D(ind1_row,ind2);
		Psi_cy = D(ind2_row,ind1);
		Psi_y = D(ind2_row,ind3);
		
		Z_c = Gamma(ind1, :) ;
		Z_x = Gamma(ind2, :) ;
		Z_y = Gamma(ind3, :) ;
		
		if TrueDict_flag == 1
			
			[dist_cx,ratio_cx] = dictdist(Psi_cx,Psi_cx0);
			[dist_cy,ratio_cy] = dictdist(Psi_cy,Psi_cy0);
			[dist_x,ratio_x] = dictdist(Psi_x,Psi_x0);
			[dist_y,ratio_y] = dictdist(Psi_y,Psi_y0);
					
			dist.cx(iter_all,1) = dist_cx ;
			dist.cy(iter_all,1) = dist_cy ;
			dist.x(iter_all,1) = dist_x ;
			dist.y(iter_all,1) = dist_y;
			
			ratio.cx(iter_all,1) = ratio_cx ;
			ratio.cy(iter_all,1) = ratio_cy ;
			ratio.x(iter_all,1) = ratio_x ;
			ratio.y(iter_all,1) = ratio_y ;
			
		end
		
		Fnorm.X = sum(sum((data(ind1_row,:) - D(ind1_row,:)*Gamma).^2));
		Fnorm.Y = sum(sum((data(ind2_row,:) - D(ind2_row,:)*Gamma).^2));
		Fnorm.XY = Fnorm.X + Fnorm.Y ;
		RMSE.XY(iter_all, 1) = sqrt( Fnorm.XY / numel(data)) ;
		RMSE.X(iter_all, 1) = sqrt( Fnorm.X / numel(data(ind1_row,:)));
		RMSE.Y(iter_all, 1) = sqrt( Fnorm.Y / numel(data(ind2_row,:)));
		

		if TrueGamma_flag
			Fnorm.Xcom = sum(sum((Psi_cx0*Z_c0 - Psi_cx*Z_c).^2));
			Fnorm.Xuniq = sum(sum((Psi_x0*Z_x0 - Psi_x*Z_x).^2));
			Fnorm.Ycom = sum(sum((Psi_cy0*Z_c0 - Psi_cy*Z_c).^2));
			Fnorm.Yuniq = sum(sum((Psi_y0*Z_y0 - Psi_y*Z_y).^2));
			
			RMSE.Xcom(iter_all, 1) = sqrt( Fnorm.Xcom / numel(data(ind1_row,:)));
			RMSE.Xuniq(iter_all, 1) = sqrt( Fnorm.Xuniq / numel(data(ind1_row,:)));
			RMSE.Ycom(iter_all, 1) = sqrt( Fnorm.Ycom / numel(data(ind2_row,:)));
			RMSE.Yuniq(iter_all, 1) = sqrt( Fnorm.Yuniq / numel(data(ind2_row,:)));
		end		
		
		% store mean sparsity
		MeanS.sc(iter_all, 1) = mean_s(1) ;
		MeanS.sx(iter_all, 1) = mean_s(2) ;
		MeanS.sy(iter_all, 1) = mean_s(3) ;

		
	  %%%%%  clear dictionary  %%%%%

	  [D,cleared_atoms] = cleardict(D,Gamma,data,muthresh,unused_sigs,replaced_atoms);
  
	  
	  %%%%%  print info  %%%%%

% 		info = sprintf('Iteration %d / %d complete', row_n, iternum);
		info = sprintf('Iteration %d / %d complete', iter, Sub_Iternum);
		if (printerr)
			info = sprintf('%s, %s = %.4f', info, errstr, err(iter));
		end
		if (printgerr)
			info = sprintf('%s, test %s = %.4f', info, errstr, gerr(iter));
		end
		if (printreplaced)
			info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
		end

% 		info = sprintf('%s, ratio(cx,x,cy,y) = (%.2f%%, %.2f%%, %.2f%%, %.2f%%)', info, ratio_cx*100, ratio_x*100, ratio_cy*100, ratio_y*100);
		info = sprintf('%s, (sc,sx,sy) = (%.2f, %.2f, %.2f) ', info, mean_s(1), mean_s(2), mean_s(3));
		if TrueDict_flag == 1
			info = sprintf('%s, ratio(cx,x,cy,y) = (%.2f%%, %.2f%%, %.2f%%, %.2f%%)', info, ratio_cx*100, ratio_x*100, ratio_cy*100, ratio_y*100);
			info = sprintf('%s, dist(cx,x,cy,y) = (%.4f, %.4f, %.4f, %.4f);', info, dist_cx, dist_x, dist_cy, dist_y);
		end
		
		info_all = sprintf('%s\n%s ', info_all, info);
		
		if (printiter)
			disp(info);
			if (msgdelta>0), disp(' '); end
		end
   
	end
	end
	
end
disp('---------------------------------------------------');
info_all = sprintf(['%s\n', '==================================================',...
				'\nCoupled Dictionary Learning completed !! \n'],info_all);
err = RMSE ;

if TrueDict_flag == 1
	output.dist = dist;
	output.ratio = ratio;
end
output.info_all = info_all;
output.MeanS = MeanS;
varargout = {output};

% 	Gamma = full(Gamma); % Change sparse mode to full mode as required by paralell training

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            optimize_atom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [atom,gamma_j,data_indices,unused_sigs,replaced_atoms] = optimize_atom(X,D,j,Gamma,unused_sigs,replaced_atoms)

global exactsvd

% data samples which use the atom, and the corresponding nonzero
% coefficients in Gamma
[gamma_j, data_indices] = sprow(Gamma, j);

if (length(data_indices) < 1)
  maxsignals = 5000;
  perm = randperm(length(unused_sigs));
  perm = perm(1:min(maxsignals,end));
  E = sum((X(:,unused_sigs(perm)) - D*Gamma(:,unused_sigs(perm))).^2);
  [d,i] = max(E);
  atom = X(:,unused_sigs(perm(i)));
  atom = atom./norm(atom);
  gamma_j = zeros(size(gamma_j));
  unused_sigs = unused_sigs([1:perm(i)-1,perm(i)+1:end]);
  replaced_atoms(j) = 1;
  return;
end

smallGamma = Gamma(:,data_indices);
Dj = D(:,j);

if (exactsvd)

  [atom,s,gamma_j] = svds(X(:,data_indices) - D*smallGamma + Dj*gamma_j, 1);
  gamma_j = s*gamma_j;
  
else
  
  atom = collincomb(X,data_indices,gamma_j') - D*(smallGamma*gamma_j') + Dj*(gamma_j*gamma_j');
  atom = atom/norm(atom);
  gamma_j = rowlincomb(atom,X,1:size(X,1),data_indices) - (atom'*D)*smallGamma + (atom'*Dj)*gamma_j;

end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             sparsecode               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gamma = sparsecode(data,D,XtX,G,thresh)

global CODE_SPARSITY codemode
global MEM_HIGH memusage
global ompfunc ompparams

if (memusage < MEM_HIGH)
  Gamma = ompfunc(D,data,G,thresh,ompparams{:});
  
else  % memusage is high
  
  if (codemode == CODE_SPARSITY)
    Gamma = ompfunc(D'*data,G,thresh,ompparams{:});
    
  else
    Gamma = ompfunc(D'*data,XtX,G,thresh,ompparams{:});
  end
  
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             compute_err              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err = compute_err(D,Gamma,data)
  
global CODE_SPARSITY codemode

if (codemode == CODE_SPARSITY)
  err = sqrt(sum(reperror2(data,D,Gamma))/numel(data));
else
  err = nnz(Gamma)/size(data,2);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           cleardict                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D,cleared_atoms] = cleardict(D,Gamma,X,muthresh,unused_sigs,replaced_atoms)

global ind1 ind2 ind3 ind1_row ind2_row

use_thresh = 4;  % at least this number of samples must use the atom to be kept

dictsize = size(D,2);

% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1)-D*Gamma(:,blocks(i):blocks(i+1)-1)).^2);
end

cleared_atoms = 0;
usecount = sum(abs(Gamma)>1e-7, 2);

for j = 1:dictsize
	if sum(j==ind1)	% if this combined atom is composed of an atom in dictionary Psi_cx and an atom in dictionary Psi_cy, do svd to update the whole combined atom.
		% compute G(:,j)
		Gj = D'*D(:,j);
		Gj(j) = 0;

		% replace atom
		if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
			[y,i] = max(err(unused_sigs));
			D(:,j) = X(:,unused_sigs(i)) / norm(X(:,unused_sigs(i)));
			unused_sigs = unused_sigs([1:i-1,i+1:end]);
			cleared_atoms = cleared_atoms+1;
		end
	elseif sum(j==ind2)	% if this combined atom is composed of a atom in dictionary Psi_x and zeros, do svd to update only the part that belongs to Psi_x.
	
		X2 = X(ind1_row,:); 
		D2 = D(ind1_row,:);
		
		% compute G(:,j)
		Gj = D2'*D2(:,j);
		Gj(j) = 0;

		% replace atom
		if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
			[err_ordered, indices_ordered] = sort(err(unused_sigs),'descend');  % maybe sort function is more time-consuming than max.
			chosen_flag = 0;
			seq_num = 0;
			while ~chosen_flag && seq_num <= length(indices_ordered)
				seq_num = seq_num + 1;
				i = indices_ordered(seq_num);
				if sum(unused_sigs(i)==ind2) % if the max
					chosen_flag = 1;
				end
			end
		
			if ~chosen_flag  % if find no suitable data in ind2 for replacement of atom
				[y,i] = max(err(unused_sigs));
			end
			
			D2(:,j) = X2(:,unused_sigs(i)) / norm(X2(:,unused_sigs(i)));
			unused_sigs = unused_sigs([1:i-1,i+1:end]);
			cleared_atoms = cleared_atoms+1;
		end
		D(ind1_row,:) = D2;
	
	else  % sum(j==ind3) % tha is sum(j==ind3)>0 % if this combined atom is composed of a atom in dictionary Psi_x and zeros, do svd to update only the part that belongs to Psi_x.
		X2 = X(ind2_row,:); 
		D2 = D(ind2_row,:);
		
		% compute G(:,j)
		Gj = D2'*D2(:,j);
		Gj(j) = 0;

		% replace atom
		if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh) && ~replaced_atoms(j) )
			[err_ordered, indices_ordered] = sort(err(unused_sigs),'descend');  % maybe sort function is more time-consuming than max.
			chosen_flag = 0;
			seq_num = 0;
			while ~chosen_flag && seq_num <= length(indices_ordered)
				seq_num = seq_num + 1;
				i = indices_ordered(seq_num);
				if sum(unused_sigs(i)==ind3) % if the max
					chosen_flag = 1;
				end
			end
		
			if ~chosen_flag  % if find no suitable data in ind3 for replacement of atom
				[y,i] = max(err(unused_sigs));
			end
			
			D2(:,j) = X2(:,unused_sigs(i)) / norm(X2(:,unused_sigs(i)));
			unused_sigs = unused_sigs([1:i-1,i+1:end]);
			cleared_atoms = cleared_atoms+1;
		end
		D(ind2_row,:) = D2;
	end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end


function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  Y(blockids) = sum(X(:,blockids).^2);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Modified OMP function            % 20151029
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w] = ModifiedOMP(b, Theta, ind1, ind2, ind3, s_c, s_x, s_y)

% 14/10/2015 
% as function pinv() takes a lot of time, it is better to call 
% ompparams = {'checkdict','off', 'messages', -1}; 
% fastest : w = ompfunc(Theta'*b,Theta'*Theta,s,ompparams). 
% fast:	w = ompfunc(Theta,b,Theta'*Theta,s,ompparams{:});
% slow:	w = ompfunc(Theta,b,[],s,ompparams{:});
% refer to	Gamma = ompfunc(D'*data,G,thresh,ompparams{:});

% 16/09/2015
% [w] = ModifiedOMP(b, Theta, ind1, ind2, ind3, s_c, s_x, s_y)
% 
% Adapts the Orthogonal Matching Pursuit algorithm (OMP) to solve
%
%     minimize    ||b - Theta*w||_2^2                   (1)
%        w
%     subject to  card(w(ind1)) <= s_c
%						card(w(ind2)) <= s_x
%						card(w(ind3)) <= s_y,
%
% where w(ind) represents the components of the vector w indexed by the
% index set ind, and card(.) represents the cardinality of a vector. 
%
% Input:
%   - b: m x 1 vector of observations
%   - Theta: m x n matrix
%   - ind1~ind3: vectors of integers between 1 and n, such that the 
%                union of their numbers is the set {1,2,...,n}
%   - s_c~s_y: positive integers
%
% Output:
%   - w: (approximate) solution of (1)
% ====================================================================== 
% Check input

[m, n] = size(Theta);

if length(b) ~= m
  error('b should have length equal to the first dimension of Theta')
end

if sum(sort([ind1, ind2, ind3]) ~= 1:n) > 0
  error('Concatenation of ind1, ind2, and ind3 should equal the set {1,2,...,n}')
end

if mod(s_c, 1) || s_c <= 0 || mod(s_x, 1) || s_x <= 0 || mod(s_y, 1) || s_y <= 0 
  error('Numbers s_c, s_x, and s_y should be positive integers')
end
% ====================================================================== 

% ====================================================================== 
% Initialization

s_w = s_c + s_x + s_y; % Total cardinality

residual   = b;
Omega = [];           % Index of chosen columns
Theta_i      = [];           % List of selected Theta_i 

counter_ind1 = 0;          % Counts the number of used Theta_i in ind1 
counter_ind2 = 0;          % Counts the number of used Theta_i in ind2 
counter_ind3 = 0;          % Counts the number of used Theta_i in ind3 
% ====================================================================== 

% ====================================================================== 
% Algorithm
for i = 1 : s_w

	% Compute the absolute value of inner products and rearrange it in descending order
	inner_product = abs(residual'*Theta);
	[inner_product_ordered, indices_ordered] = sort(inner_product, 'descend');  

	% Determine the first index in indices_ordered that can be used

	chosen_index = [];
	iter = 0;       % Goes through indices_ordered

	while isempty(chosen_index)
		
		iter = iter + 1;
		candidate_index = indices_ordered(iter);

		% Find where candidate_index belongs
		if sum(candidate_index == ind1) && counter_ind1 < s_c % candidate_index is a scalar while ind1 is a vector !!??

			chosen_index = candidate_index; 
			counter_ind1 = counter_ind1 + 1;

		elseif sum(candidate_index == ind2) && counter_ind2 < s_x

			chosen_index = candidate_index; 
			counter_ind2 = counter_ind2 + 1;

		elseif sum(candidate_index == ind3) && counter_ind3 < s_y

			chosen_index = candidate_index; 
			counter_ind3 = counter_ind3 + 1;

		end	
	end

	Omega = [Omega, chosen_index];          
	Theta_i      = [Theta_i, Theta(:, chosen_index)];
	Theta(:, chosen_index) = zeros(m,1);
	w_i = pinv(Theta_i)*b;     
	residual = b - Theta_i*w_i;

end 
% -------------------------------------------------------------------------------------------------------------------
w = zeros(n,1);
w(Omega) = w_i;
end
