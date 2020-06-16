function est = sliding_window(xfm, d, Nt, t0, window)

%   xfm is an object of the transform class
%   d is measured data (output dimensions equal to size(xfm*imdata)
%   Nt is the total number of projections in the data
%   t0 are the centre points for each window
%   window is the odd window size

switch xfm.mode
case {'FFT','TURBINE'}
    disp('Not implemented yet');
    est =   [];
case {'NUFFT'}

    dims    =   xfm.Nd;
    dims(4) =   Nt;

    nc      =   xfm.ncoils;
    sMat    =   xfm.S;

    k       =   reshape(xfm.k,[],Nt,2);
    d       =   reshape(permute(d,[2,1,3]),nc, [],Nt);

    for i = 1:length(t0)
        disp(sprintf('Time Point %i',i));
        idx = max(1,t0(i)-(window-1)/2):min(Nt,t0(i)+(window-1)/2);
        dd  =   permute(reshape(d(:,:,idx),nc,[],1),[2,1,3]);
        xx  =   transform([dims(1:2) 1 1],'NUFFT','k',reshape(k(:,idx,:),[],1,2),'sensMat',sMat);
        est0=   tmean(xx,dd);
        xx  =   transform([dims(1:2) 1 1],'NUFFT','k',reshape(k(:,idx,:),[],1,2),'sensMat',sMat,'wi',1);
        est(:,:,:,i)    =   cg_SENSE(xx,dd, 1E-3, 10,est0);
    end

end
