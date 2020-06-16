function res = weight(a,b)

%   Pre-weight k-space data to accommodate weighting scheme used in transform
%   Only in NUFFT or TURBINE modes
switch a.mode
case 'FFT'
    res =   b;
case 'NUFFT'
    res =   b.*repmat(reshape(sqrt(a.w),[],1,a.Nt),[1,a.ncoils,1]);
case 'TURBINE'
    res =   b;
    for c = 1:a.ncoils
        res(c,:,:,:)  =   shiftdim(res(c,:,:,:)).*permute(sqrt(a.w),[1,3,2]);
    end
end
