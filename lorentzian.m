
function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
% g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
g = ones(size(x)) ./ ( 4* (x-position).*(x-position) + width.^2);
% g = g./trapz(x,g,2);
end