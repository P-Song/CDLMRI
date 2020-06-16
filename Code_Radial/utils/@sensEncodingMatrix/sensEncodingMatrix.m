classdef sensEncodingMatrix
%   Generates sparse stripe-diagonal sensitivity encoding matrix
%
%   Last Modified:
%   Mark Chiew
%   June 2014
%
%   Input:
%           sens    -   ND array of complex coil sensitivities
%                       Assumes 1st dimension is coils
%                       [Nc, Nx, (Ny), (Nz)]

properties (SetAccess = private, GetAccess = private)

    adjoint =   0;
    dims    =   [];
    coils   =   [];
    ccoils  =   [];
    sos     =   [];
end

properties (SetAccess = private, GetAccess = public)
    ncoils  =   0;
end

methods
function res = sensEncodingMatrix(coils)

    res.dims    =   size(coils);
    res.dims(4) =   size(coils,4);
    res.coils   =   permute(reshape(coils,res.dims(1),res.dims(2),res.dims(3),size(coils,4)),[2,3,4,1]);
    res.ccoils  =   permute(reshape(conj(coils),size(coils,1),[]),[2,1]);
    res.sos     =   reshape(sum(conj(coils).*coils,1),[],1);

    res.ncoils  =   size(coils,1);

end
end
end
