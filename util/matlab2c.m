% to convert matlab structure to C job file.
% function success=matlab2c(seq, filename)
function success=matlab2c(seq, filename)
f=fopen(filename,'w');

fprintf(f,'Numberofintervals=%d\n',seq.npulse);
fprintf(f,'Timeintervals/ms=\n');
fprintf(f,'%g\n',seq.time);
fprintf(f,'Appliedgradientsx/(mT/m)=\n');
fprintf(f,'%g\n',seq.gradx);
fprintf(f,'Appliedgradientsy/(mT/m)=\n');
fprintf(f,'%g\n',seq.grady);
fprintf(f,'Appliedgradientsz/(mT/m)=\n');
fprintf(f,'%g\n',seq.gradz);
fprintf(f,'Flipangles/grad=\n');
fprintf(f,'%g\n',seq.angle);
fprintf(f,'Flipaxis/grad=\n');
fprintf(f,'%g\n',seq.axes);
fprintf(f,'Diffusioncoeff/(mcm^2/ms)=%g\n',seq.diff);
fprintf(f,'RelaxationtimeT1/ms=%g\n',seq.T1);
fprintf(f,'RelaxationtimeT2/ms=%g\n',seq.T2);
fprintf(f,'Offresonance/Hz=%g\n',seq.Omega);
fprintf(f,'Tolerance_k=%g\n',seq.ktolerance);
fprintf(f,'Tolerance_k_phys=%g\n',seq.ktolerance_phys);


fprintf(f,'Velocityx/(mm/s)=\n');
if(isfield(seq,'velocityx'))
 fprintf(f,'%g\n',seq.velocityx);
else
 velox=zeros(seq.npulse,1);
 fprintf(f,'%g\n',velox);
end
fprintf(f,'Velocityy/(mm/s)=\n');
if(isfield(seq,'velocityy'))
 fprintf(f,'%g\n',seq.velocityy);
else
 veloy=zeros(seq.npulse,1);
 fprintf(f,'%g\n',veloy);
end
fprintf(f,'Velocityz/(mm/s)=\n');
if(isfield(seq,'velocityz'))
 fprintf(f,'%g\n',seq.velocityz);
else
 veloz=zeros(seq.npulse,1);
 fprintf(f,'%g\n',veloz);
end
fprintf(f,'nOutput=\n');
if(isfield(seq,'nOutput'))
 fprintf(f,'%g\n',seq.nOutput);
else
 nOutput=0;
 fprintf(f,'%g\n',nOutput);
end
if(isfield(seq,'echo'))
fprintf(f,'Echointervals/ms=\n');
fprintf(f,'%g\n',seq.echo);    
end

fclose(f);
success=1;
