function b = ifftfn(a,b,dims)

if a.centric
    for i = dims
        b   =   fftshift(ifft(ifftshift(b, i), [], i), i);
    end
else
    for i = dims
        b   =   ifft(b, [], i);
    end
end
