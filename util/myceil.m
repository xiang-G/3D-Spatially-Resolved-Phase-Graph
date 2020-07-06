function out=myceil(in)
% abs(in)<1     out 0
% abs(in)>1     out +/- ceil(in)
in(abs(in)<1)=0;
out=sign(in).*ceil(abs(in));
end