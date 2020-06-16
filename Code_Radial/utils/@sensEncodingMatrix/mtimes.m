function res = mtimes(a,b)

if a.adjoint
    nc  =   a.dims(1);
    nt  =   size(b,3);
    res =   zeros(size(b,1),1,nt);
    for c = 1:a.dims(1)
        res =   res + b(:,c,:).*repmat(a.ccoils(:,c),[1,1,nt]);
    end
    res =   reshape(res,[],nt)./repmat(a.sos,[1,nt]);
else
    b   =   reshape(b,a.dims(2),a.dims(3),a.dims(4),1,[]);
    nc  =   a.dims(1);
    nt  =   size(b,5);
    res =   zeros([a.dims(2:end),nc,nt]);
    for c = 1:nc
        res(:,:,:,c,:)  =   b.*repmat(a.coils(:,:,:,c),[1,1,1,1,nt]);
    end
end
