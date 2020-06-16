function [k2, wi, varargout] = gen_radial(incr, N_r, N_proj, N_meas, mod_angle, cont, shift, r)

%   If incr = 0, means Golden Ratio-based projection spacing

if ~isa(incr, 'function_handle')
    if incr == 0
        ang_fn  =   @(x) (x-1)*360/(1+sqrt(5));
    else
        ang_fn  =   @(x) incr*(x-1);
    end
else
    ang_fn  =   incr;
end

%   Generate radial spacing information
dr  =   1/N_r;
if nargin < 8 
    r   =   (-0.5:dr:(N_r-1)*dr-0.5);
end
if nargin < 7
   shift        =   0.5/N_r; 
end
if nargin < 5
    mod_angle   =   180;
    cont        =   0;
end

if cont ~= 0
    N_proj  =   N_proj*N_meas;
    N_meas  =   1;
end

%   Generate angles
angles  =   ang_fn(1:N_proj);
angles  =   repmat(angles, [1, N_meas]);


%r   =   linspace(-0.5,0.5,N_r);
r   =   r+shift;

%   Initialise output arrays
k   =   zeros(N_r, N_proj*N_meas);

%   Populate output arrays
for i = 1:N_proj*N_meas
   k(:,i)   =   (2*pi*r)*exp(+1j*pi*angles(i)/180);
end

%   Generate density compensation weighting 
%   Weight every point by |r| (dr, dthetha const) 
wi  =   abs(k);

k2  =   zeros(N_r, N_proj*N_meas, 2);
k2(:,:,1)   =   real(k);
k2(:,:,2)   =   imag(k);

if nargout == 3
    varargout{1}    =   angles;
else
    varargout{1}    =   angles;

    q   =   sort(mod(angles,180))*pi/180;
    E   =   ((q(2)-(q(end)-pi))/2).^2;
    for i = 2:length(q)-1
        E   =   E + ((q(i+1)-q(i-1))/2).^2;
    end
    E   =   E + ((q(1)-(q(end-1)-pi))/2).^2;
    varargout{2}    =   1/sqrt(E);
end
