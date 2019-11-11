function out=myceil(in)
% abs(in)<1     out 0
% abs(in)>1     out +/- ceil(in)
if(abs(in))<1
    out=0;
else
    out=sign(in)*ceil(abs(in));
end
end