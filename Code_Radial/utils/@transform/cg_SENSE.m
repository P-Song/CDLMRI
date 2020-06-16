function est = cg_SENSE(xfm, d, tol, iters, est)

%   Performs iterative conjugate gradient SENSE recon
%   Input d should be shaped like the output of mtimes

if nargin == 2
    tol     =   1E-4;
    iters   =   10;
    est     =   zeros(prod(xfm.Nd),xfm.Nt);
end

switch xfm.mode
case 'FFT'
    disp('Not implemented yet');
    est =   [];

case {'NUFFT','TURBINE'}
    a   =   0.05;
    b   =   0.1;
    err =   norm(d(:));

    k   =   0;
    g   =   2*(xfm'*(xfm*est-d));
    dx  =   -g;

    while norm(g(:)) > tol*err && k < iters
        t   =   1;
        while norm(reshape(d-xfm*(est+t*dx),[],1))^2 > norm(reshape(d-xfm*est,[],1))^2+a*t*real(reshape(g,[],1)'*reshape(dx,[],1))
            t   =   b*t;
        end
        est =   est+t*dx;
        g2  =   2*(xfm'*(xfm*est-d));
        dx  =   -g2+dx*(norm(g2(:))/norm(g(:)))^2;
        k   =   k+1;
        g   =   g2;
        disp(sprintf('Iter: %03i, Grad: %f', k, norm(g(:))/err));
    end

    est =   reshape(est,[xfm.Nd, xfm.Nt]);

case 'TURBINE'
    disp('Not implemented yet');
    est =   [];

end
