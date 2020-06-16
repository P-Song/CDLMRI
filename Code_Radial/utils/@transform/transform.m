function  res = transform(dims, mode, varargin)

%   Transform class that encapsulates the nufft code by Fessler or regular FFT
%   Handles time dimension data
%   Handles both Cartesian (FFT) and non-Cartesian (NUFFT) trajectories
%   Handles single channel and multi-coil data
%
%   NOTE:   Requires the nufft portion of Jeff Fessler's irt toolbox
%           See http://web.eecs.umich.edu/~fessler/irt/fessler.tgz
%   NOTE2:  Using old-style MATLAB classes because MATLAB OOP is horrible and
%           new-classdef classes are extremely slow with property access and
%           have very high overheads for method calls.
%
%   Last Modified:
%   Mark Chiew
%   July 2014
%
%   Class is modelled after the MCNUFFT class by
%   Li Feng & Ricardo Otazo, NYU, 2012
%
%   Global Required Inputs:
%       dims    =   [Nx, Ny, Nz, Nt] 4D vector of image dimensions
%       mode    =   "FFT" or "NUFFT"
%
%   Global Optional Inputs:
%       coils   =   [Nc, Nx, Ny, Nz] array of coil sensitivities 
%                   Defaults to single channel i.e. ones(Nx,Ny,Nz)
%       verbose =   0 or 1 toggle for extra verbose output
%                   Defaults to 0 (off)
%
%   This class can be called in FFT, NUFFT or TURBINE mode:
%
%   FFT Mode 
%       Optional:
%               mask    =   [Nx,Ny,Nz,Nt] logical mask of sampled points
%                           Defaults to full sampling
%               centric =   1 or 0, indicating whether the FFT is centric
%                           Defaults to 0 (non-centric FFT is faster)
%                           It is important to make sure the sampling mask is
%                           consistent with this parameter
%
%   NUFFT Mode
%       Required Inputs:
%               k       =   [Nsamp, Nt, 2] (2D) or [Nsamp, Nt, 3] (3D) 
%                           array of sampled k-space positions
%                           (in radians, normalised range -pi to pi)
%               wi      =   Density compensation weights
%       Optional:
%               Jd      =   Size of interpolation kernel
%                           Defaults to [6,6]
%               Kd      =   Size of upsampled grid
%                           Defaults to 200% of image size
%               shift   =   NUFFT shift factor
%                           Defaults to 50% of image size
%
%   TURBINE Mode (Special case, radial kx-ky, cartesian kz)
%       Required Inputs:
%               k       =   [Nsamp, Nt, 2] array of sampled k-space positions
%                           (in radians, normalised range -pi to pi)
%               wi      =   Density compensation weights
%       Optional:
%               Jd      =   Size of interpolation kernel
%                           Defaults to [6,6]
%               Kd      =   Size of upsampled grid
%                           Defaults to 200% of image size
%               shift   =   NUFFT shift factor
%                           Defaults to 50% of image size
%
%   Usage:
%           Forward transforms can be applied using the "*" or ".*" operators
%           or equivalently the "mtimes" or "times" functions
%           The input can be either the "Casorati" matrix (Nx*Ny*Nz, Nt) or
%           an n-d image array
%           The reverse transform can be accessed by transposing the transform
%           object 
%           The difference between "*" and ".*" is only significant for the 
%           reverse transform, and changes  the shape of the image ouptut:
%               "*" produces data that is in matrix form (Nx*Ny*Nz, Nt) whereas
%               ".*" produces data in n-d  array form (Nx, Ny, Nz, Nt) 
%
%   Example usage:
%   
%   img =   some_image_timeseries(64,64,64,512);
%   xfm =   transform(dims, 'FFT', 'mask', sampling_mask);
%   k   =   xfm*image;
%   k   =   xfm*reshape(image,[],512); %same as previous line
%
%   q_matrix    =   xfm'*k;
%   q_array     =   xfm'.*k;

%===========================================================

    %   Input parser object
    p   =   inputParser;

    %   Input validation functions
    modeValidator   =   @(x) sum(strcmp(x,{'FFT','NUFFT','TURBINE'})) == 1;
    coilValidator   =   @(x) isequal([size(x,2),size(x,3),size(x,4)],dims(1:3));
    toggleValidator =   @(x) x == 0 || x == 1;
    lengthValidator =   @(x) length(x) == 2 || length(x) == 3;

    %   Global options
    p.addRequired('dims', @(x) length(x) == 4);
    p.addRequired('mode', modeValidator);

    p.addParamValue('sensMat',  [],                     @(x) strcmp(class(x),'sensEncodingMatrix'));
    p.addParamValue('coils',    ones([1, dims(1:3)]),   coilValidator);
    p.addParamValue('verbose',  0,                      toggleValidator);
    p.addParamValue('norm',     1,                      toggleValidator);

    %   FFT-specific options
    p.addParamValue('mask',     logical(ones(dims)), @(x) numel(x) == prod(dims));
    p.addParamValue('centric',  0,          toggleValidator);

    %   NUFFT-specific options
    p.addParamValue('k',    [],         @(x) isequal(size(x,2),dims(4))); 
    p.addParamValue('wi',   [],         @(x) size(x,2) == dims(4)||isscalar);
    p.addParamValue('Jd',   [6,6,6],    lengthValidator);
    p.addParamValue('Kd',   2*dims(1:3),lengthValidator);
    p.addParamValue('shift',dims(1:3)/2,lengthValidator);

    p.parse(dims, mode, varargin{:});
    p   =   p.Results;

    %   Initialize properties (MATLAB OOP-restriction)
    res.adjoint =   0;
    res.Nd      =   p.dims(1:3);
    res.Nt      =   p.dims(4);
    res.mode    =   p.mode;
    res.centric =   p.centric;
    res.verbose =   p.verbose;
    if isempty(p.sensMat)
        res.S   =   sensEncodingMatrix(p.coils);
        res.ncoils  =   size(p.coils, 1);
    else
        res.S   =   p.sensMat;
        res.ncoils  =   res.S.ncoils;
    end

    res.mask    =   [];
    res.st      =   [];
    res.dataSize=   [];
    res.w       =   [];
    res.k       =   [];
    res.w_mean  =   [];
    res.st_mean =   [];
    res.norm    =   [];

    switch res.mode
        case 'FFT'

            res.mask    =   repmat(reshape(p.mask,[res.Nd(1:3),1,res.Nt]), 1,1,1,res.ncoils,1);

        case 'NUFFT'
            
            res.k   =   p.k;

            disp('Initialising NUFFT');
            for t = 1:res.Nt
                if res.Nd(3) == 1
                %   2D
                    res.st{t}   =   nufft_init(squeeze(p.k(:, t, 1:2)),...
                                               res.Nd(1:2),...
                                               p.Jd(1:2),...
                                               p.Kd(1:2),...
                                               p.shift(1:2));
                else    
                %   3D (use 'table' option to save memory)
                    res.st{t}   =   nufft_init(squeeze(p.k(:, t, :)),...
                                               res.Nd,...
                                               p.Jd,...
                                               p.Kd,...
                                               p.shift,...
                                               'table', 2^12);
                end
            end

            %   For computing temporal mean
            if res.Nd(3) == 1
                res.st_mean =   nufft_init(reshape(p.k(:, :, 1:2), [], 2),...
                                           res.Nd(1:2),...
                                           p.Jd(1:2),...
                                           p.Kd(1:2),...
                                           p.shift(1:2));
            else
                res.st_mean =   nufft_init(reshape(p.k, [], 3),...
                                           res.Nd,...
                                           p.Jd,...
                                           p.Kd,...
                                           p.shift,...
                                           'table', 2^12);
            end

            res.dataSize    =   size(p.k);
            if isempty(p.wi)
                disp('Generating Density Compensation Weights');
                %res.w   =   ones(size(p.k(:,:,1)));
                for t = 1:res.Nt
                    res.w(:,t)  =   ones(size(p.k,1),1);
                    for ii = 1:20
                        tmp =   res.st{t}.p*(res.st{t}.p'*res.w(:,t));
                        res.w(:,t)  =   res.w(:,t)./real(tmp);
                    end
                end
                res.w_mean  =   ones(size(res.st_mean.p,1),1);
                pp      =   res.st_mean.p;
                for ii  =   1:5
                    tmp =   pp*(pp'*res.w_mean);
                    res.w_mean   =   res.w_mean./real(tmp);
                end
                res.w_mean  =   res.st_mean.sn(end/2,end/2)^(-2)/prod(res.st_mean.Kd)*res.w_mean;
            elseif isscalar(p.wi)
                res.w       =   repmat(p.wi,size(p.k(:,:,1)));
                res.w_mean  =   res.w(:);
            else
                res.w       =   p.wi;
                res.w_mean  =   p.wi(:);
            end
            
            if p.norm
                %res.norm        =   0.5/numel(res.w(:,1));
                res.norm        =   res.st{1}.sn(end/2,end/2)^(-2)/prod(res.st{1}.Kd);
            else
                res.norm        =   1;
            end

        case 'TURBINE'

            res.k   =   p.k;

            disp('Initialising NUFFT');
            for t = 1:res.Nt
                res.st{t}   =   nufft_init(squeeze(p.k(:, t, 1:2)),...
                                           res.Nd(1:2),...
                                           p.Jd(1:2),...
                                           p.Kd(1:2),...
                                           p.shift(1:2));
            end

            res.st_mean     =   nufft_init(reshape(p.k(:,:,1:2), [], 2),...
                                           res.Nd(1:2),...
                                           p.Jd(1:2),...
                                           p.Kd(1:2),...
                                           p.shift(1:2));

            res.dataSize    =   size(p.k);

            if isempty(p.wi)
                disp('Generating Density Compensation Weights');
                %res.w   =   ones(size(p.k(:,:,1)));
                for t = 1:res.Nt
                    res.w(:,t)  =   ones(size(p.k,1),1);
                    for ii = 1:20
                        tmp =   res.st{t}.p*(res.st{t}.p'*res.w(:,t));
                        res.w(:,t)  =   res.w(:,t)./real(tmp);
                    end
                end
                res.w       =   repmat(res.w, [1,1,res.Nd(3)]);
                res.w_mean  =   ones(size(res.st_mean.p,1),1);
                pp      =   res.st_mean.p;
                for ii  =   1:5
                    tmp =   pp*(pp'*res.w_mean);
                    res.w_mean   =   res.w_mean./real(tmp);
                end
                res.w_mean  =   res.st_mean.sn(end/2,end/2)^(-2)/prod(res.st_mean.Kd)*res.w_mean;
                res.w_mean  =   repmat(res.w_mean,[1,1,res.Nd(3)]);
            else
                res.w       =   repmat(p.wi,[1,1,res.Nd(3)]);
                res.w_mean  =   repmat(p.wi(:),[1,1,res.Nd(3)]);
            end

            if p.norm
                res.norm        =   0.5/numel(res.w(:,1));
            else
                res.norm        =   1;
            end


    end
    
    res =   class(res, 'transform');
end
