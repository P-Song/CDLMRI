function res = times(a,b)

if a.adjoint
%   Inverse (NU)FFT
    res =   reshape(mtimes(a,b), [a.Nd, a.Nt]);
else
%   Forward (NU)FFT and sampling
    res =   mtimes(a,b);
end
