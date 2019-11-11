function out=myfloor(in)
% abs(in)<1     out 0
% abs(in)>1     out +/- ceil(in)
if(abs(in))<1
    out=0;
else
    out=sign(in)*floor(abs(in));
end
end