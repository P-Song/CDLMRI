function b = fftfn(a,b,dims)

if a.centric
    for i = dims
        b   =   fftshift(fft(ifftshift(b, i), [], i), i);
    end
else
    for i = dims
        b   =   fft(b, [], i);
    end
end
