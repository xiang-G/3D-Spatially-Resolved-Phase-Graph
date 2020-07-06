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
fprintf(f,'RelaxationtimeT1/ms=%g\n',seq.T1);
fprintf(f,'RelaxationtimeT2/ms=%g\n',seq.T2);
fprintf(f,'Offresonance/Hz=%g\n',seq.Omega);
fprintf(f,'Tolerance_k=%g\n',seq.ktolerance);
fprintf(f,'Tolerance_k_phys=%g\n',seq.ktolerance_phys);
input_diff(f,seq);
input_velocity(f,seq);
fprintf(f,'nOutput=\n');
if(isfield(seq,'nOutput'))
    fprintf(f,'%g\n',seq.nOutput);
else
    fprintf(f,'%g\n',0);
end
if(isfield(seq,'echo'))
    fprintf(f,'Echointervals/ms=\n');
    fprintf(f,'%g\n',seq.echo);
end
fclose(f);
success=1;
end

function input_diff(f,seq)
if isfield(seq,'diff')
    if seq.diff>=0
        fprintf(f,'Diffusioncoeff/(mcm^2/ms)=%g\n',seq.diff);
    else
        error('Diffusioncoeff should not be negative');
    end
elseif isfield(seq,'diffTensor')
    if all(seq.diffTensor(:)>=0) && isequal(reshape(seq.diffTensor,3,3),transp(reshape(seq.diffTensor,3,3)))
        fprintf(f,'DiffusionTensorcoeff/(mcm^2/ms)=\n');
        fprintf(f,'%g\n',seq.diffTensor);
    else
        error('Diffusioncoeff should not be negative; or non-symmetric');
    end
end
end

function input_velocity(f,seq)
fprintf(f,'Velocityx/(mm/s)=\n');
if(isfield(seq,'velocityx'))
    fprintf(f,'%g\n',seq.velocityx);
else
    fprintf(f,'%g\n',zeros(seq.npulse,1));
end
fprintf(f,'Velocityy/(mm/s)=\n');
if(isfield(seq,'velocityy'))
    fprintf(f,'%g\n',seq.velocityy);
else
    fprintf(f,'%g\n',zeros(seq.npulse,1));
end
fprintf(f,'Velocityz/(mm/s)=\n');
if(isfield(seq,'velocityz'))
    fprintf(f,'%g\n',seq.velocityz);
else
    fprintf(f,'%g\n',zeros(seq.npulse,1));
end
end
