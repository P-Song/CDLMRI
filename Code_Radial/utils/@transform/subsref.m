function b = subsref(xfm, s)

switch s.type
case '.'
    switch s.subs
    case 'k'
        b = xfm.k;
    case 'dims'
        b = [xfm.Nd xfm.Nt];
    case 'S'
        b = xfm.S;
    case 'mode'
        b = xfm.mode;
    case 'ncoils'
        b = xfm.ncoils;
    case 'dataSize'
        b = [xfm.dataSize(1) xfm.ncoils xfm.dataSize(2)];
    case 'msize'
        b = [prod(xfm.Nd) xfm.Nt];
    case 'w'
        b = xfm.w;
    end
end
