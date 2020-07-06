function out=myfloor(in)
% abs(in)<1     out 0
% abs(in)>1     out +/- floor(in)
in(abs(in)<1)=0;
out=sign(in).*floor(abs(in));
end