function res = mtimes(a,b)

switch a.mode
case 'FFT'
    if a.adjoint
    %   Inverse FFT
        res =   zeros([a.Nd, a.ncoils, a.Nt]);
        res(a.mask)  =   b(:);
        res =   ifftfn(a, res, 1:3);
        res =   reshape(res, [prod(a.Nd), a.ncoils, a.Nt]);
        res =   a.S'*res;
    else
    %   Forward FFT and sampling
        b   =   a.S*b;
        b   =   fftfn(a, b, 1:3);
        res =   b(a.mask);
    end
case 'NUFFT'
    if a.adjoint
    %  Reverse NUFFT to generate Images on Cartesian grid
        res =   zeros([prod(a.Nd), a.ncoils, a.Nt]);

        for t = 1:a.Nt
            res(:,:,t)  =   sqrt(a.norm)*reshape(nufft_adj(repmat(sqrt(a.w(:,t)),1,a.ncoils).*b(:,:,t), a.st{t}), [], a.ncoils,1);
        end

        res =   a.S'*res;
    else
    %  Forward NUFFT to generate k-space data on non-uniform points
        res =   zeros(a.dataSize(1), a.ncoils, a.Nt);
        b   =   a.S*b;

        for t = 1:a.Nt
            res(:,:,t)  =   sqrt(a.norm)*reshape(nufft(b(:,:,:,:,t), a.st{t}),[],a.ncoils).*repmat(sqrt(a.w(:,t)),1,a.ncoils);
        end
    end
case 'TURBINE'
    if a.adjoint
    %   Reverse transform by doing 2D inverse NUFFT for each kz and then iFFTz
        res =   zeros([a.Nd, a.ncoils, a.Nt]);

        %   NUFFT along kx-,ky-directions
        for c = 1:a.ncoils, for t = 1:a.Nt
            if a.verbose
                printf('Coil: %02d, Time: %04d', c, t);
            end
            res(:,:,:,c,t)  =   a.norm*nufft_adj(squeeze(a.w(:,t,:)).*shiftdim(b(c,:,:,t)), a.st{t});

        end, end

        res =   a.S'*reshape(res, [], a.ncoils, a.Nt);

    else
    %   Forward transform by performing FFTz and then 2D NUFFT on (x,y,kz) data
        res =   zeros(a.ncoils, a.dataSize(1), a.Nd(3), a.Nt);
        b   =   reshape(a.S*b, [a.Nd, a.ncoils, a.Nt]);

        %   NUFFT along x-,y-directions
        for c = 1:a.ncoils, for t = 1:a.Nt
            if a.verbose
                printf('Coil: %02d, Time: %04d', c, t);
            end
            res(c,:,:,t)    =   nufft(squeeze(b(:,:,:,c,t)), a.st{t});
        end, end
    end
case 'MB-FFT'
    if a.adjoint
    %   Inverse FFT
    else
    %   Forward FFT and sampling
    end
end
