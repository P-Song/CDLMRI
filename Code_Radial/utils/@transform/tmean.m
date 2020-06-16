function res = tmean(a,b)

%   Computes temporal mean of (under)sampled k-t data and returns an image
%   Input should be shaped exactly like the output of mtimes

switch a.mode
case 'FFT'
    res     =   zeros([prod(a.Nd)*a.ncoils, a.Nt]);
    mask    =   false([prod(a.Nd)*a.ncoils, a.Nt]);

    res(a.mask)  =   b;
    mask(a.mask) =   true;

    res     =   reshape(sum(res, 2)./(sum(mask, 2)+eps), [a.Nd, a.ncoils]);
    res     =   reshape(ifftfn(a, res, 1:3), [prod(a.Nd), a.ncoils]);
    res     =   a.S'*res;

case 'NUFFT'
    res     =   zeros([prod(a.Nd), a.ncoils]);
    b       =   reshape(permute(b,[1,3,2]),[], a.ncoils);


    res =   sqrt(a.norm)*nufft_adj(repmat(sqrt(a.w_mean),1,a.ncoils).*b, a.st_mean);
    res =   a.S'*reshape(res, [], a.ncoils);

case 'TURBINE'
    res =   zeros([a.Nd, a.ncoils]);
    b   =   reshape(permute(b,[1,2,4,3]), a.ncoils, [], a.Nd(3));

    w       =   ones(size(a.st_mean.p,1),1);
    for ii  =   1:20
        tmp =   a.st_mean.p*(a.st_mean.p'*w);
        w   =   w./real(tmp);
    end
    w       =   repmat(w,[1,a.Nd(3)]);
    w       =   a.st_mean.sn(end/2,end/2)^(-2)/prod(a.st_mean.Kd)*w;

    %   NUFFT along kx-,ky-directions
    for c = 1:a.ncoils
        res(:,:,:,c)  =   nufft_adj(w.*squeeze(b(c,:,:)), a.st_mean);
        res(:,:,:,c)  =   res(:,:,:,c)*a.norm/a.Nt;
    end

    res =   a.S'*reshape(res, [], a.ncoils);

end
