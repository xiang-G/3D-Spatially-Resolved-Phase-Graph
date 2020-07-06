% function [Kgrid3D, Bgrid3D]=prep_3dkmap(slogData,dK)
% in:                   slogData; output from KmapSignal2 function
%                       dK;       input post-gridding bins [1,3]
% out:                  Kgrid3D;  3D(K-vector grids) + 1D(time)
%                       Bgrid3D;  3D(Background grids) + 1D(time)
%
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

function [Kgrid3D, Bgrid3D]=prep_3dkmap(slogData,dK)
global Mat;
global Kmax2Pi;     %Skope of 3D-K space

%%
NE=length(slogData.NKperEcho)-1;
NKperEcho=slogData.NKperEcho;
kvector = slogData.kvector./2/pi*1e6; %1/m
signal = abs(slogData.signal(:,1)+1i*slogData.signal(:,2));

if Kmax2Pi
    Kmax=[1,1,1]*Kmax2Pi*1e3; %1/m
else
    Kmax=max(abs(kvector(1:end,:))); %1/m
end
disp(['Kmax for 3dkmap ',int2str(Kmax)]);
Kgrid3D=zeros([2*ceil(Kmax./dK)+1,NE]);
for necho=1:NE
    Kgrid=zeros(2*ceil(Kmax./dK)+1);
    shift=abs(ceil(Kmax./dK)+1);
    for i=NKperEcho(necho)+1:NKperEcho(necho+1)
        bins = myfloor(kvector(i,:)./dK)+shift;
        if all(bins>0) && all(bins<2*ceil(Kmax./dK)+1)
            Kgrid(bins(1),bins(2),bins(3))=Kgrid(bins(1),bins(2),bins(3))+signal(i);
        end
    end
    Kgrid3D(:,:,:,necho) = Kgrid;
end

%% receive zone
Bgrid3D = zeros(size(Kgrid3D));
bg1_st = ceil((size(Kgrid3D,1)-Mat)/2);
bg2_st = ceil((size(Kgrid3D,2)-Mat)/2);
bg3_st = ceil((size(Kgrid3D,3)-Mat)/2);
if Mat<=1
    for necho=1:NE
        Bgrid3D(bg1_st+1,bg2_st+1,bg3_st+1,necho)=1;
    end
else
    for necho=1:NE
        Bgrid3D(bg1_st:Mat/2:bg1_st+Mat,bg2_st:Mat/2:bg2_st+Mat,bg3_st:Mat/2:bg3_st+Mat,necho)=1;
    end
end

end