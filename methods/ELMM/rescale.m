function out = rescale(in,a)

m = min(in(:));
M = max(in(:));
out = a*(in-m)/(M-m);
if a == 1
    out = double(out);
end
if a == 255
    out = uint8(out);
end