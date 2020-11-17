%% IMG_DFT
function slogData=IMG_DFT(varargin)
slogData = varargin{1};
global FOV;
global Mat;
global TDImage;
res = FOV/Mat;
%% simulate signal on which position
if nargin <2 || isempty(varargin{2})
        FOVx=-FOV/2+res:res:FOV/2;
        FOVy=-FOV/2+res:res:FOV/2;
        FOVz=-FOV/2+res:res:FOV/2;
else
    FOVp= varargin{2};
    if size(FOVp,1)==1
        FOVp=FOVp';
        FOVx=FOVp(1,:);FOVy=FOVp(2,:);FOVz=FOVp(3,:);
    else
        FOVx=FOVp(1,:);FOVy=FOVp(2,:);FOVz=FOVp(3,:);
    end
end
N=length(FOVx);
NE=length(slogData.NKperEcho)-1;
NKperEcho=slogData.NKperEcho;
%% signal indensity (x dim without deltaG, y dim with deltaG)
if ~TDImage
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]);
    pos = zeros(3,N);
    posx=pos;posx(1,:)=FOVx;%posx(2,:)=FOVy;posx(3,:)=0;   
    sigx = zeros(N,NE); 
    for i=1:NE
        [sigx(:,i)]=DFT(slogData,posx*1e3,res*1e3,NKperEcho(i:i+1)); 
    end
    %write back
    slogData.sigxy = (sigx);
    slogData.sig3d = abs(sigx);
    slogData.pha3d = angle(sigx);
else
    [posxy,posyz,posxz]= deal(zeros(3,N*N));
    [sigxy,sigyz,sigxz]= deal(zeros(N*N,NE));
    disp(['K_filter = ',num2str(2*pi/(res*1e3)/2)]); % don't foget to divid by 2

    %% 2D Image (x-y)
    posxy(1,:)=repmat(FOVx,[1,N]);
    posxy(2,:)=col(repmat(FOVy',[1,N])');
    posxy(3,:)=0;
    for i=1:NE
        [sigxy(:,i)]=DFT(slogData,posxy*1e3,res*1e3,NKperEcho(i:i+1)); 
    end
    ImgPerExy = abs(reshape(sigxy,[N,N,NE]));%
    PhasePerExy = angle(reshape(sigxy,[N,N,NE]));%
    
    %% 2D Image (y-z)
    posyz(1,:)=0;
    posyz(2,:)=repmat(FOVy,[1,N]);
    posyz(3,:)=col(repmat(FOVz',[1,N])');
    for i=1:NE
        [sigyz(:,i)]=DFT(slogData,posyz*1e3,res*1e3,NKperEcho(i:i+1)); 
    end
    ImgPerEyz = abs(reshape(sigyz,[N,N,NE]));%
    PhasePerEyz = angle(reshape(sigyz,[N,N,NE]));%
    
    %% 2D Image (x-z)
    posxz(1,:)=repmat(FOVx,[1,N]);
    posxz(2,:)=0;
    posxz(3,:)=col(repmat(FOVz',[1,N])');
    for i=1:NE
        [sigxz(:,i)]=DFT(slogData,posxz*1e3,res*1e3,NKperEcho(i:i+1)); 
    end
    ImgPerExz = abs(reshape(sigxz,[N,N,NE]));%
    PhasePerExz = angle(reshape(sigxz,[N,N,NE]));%
    
    %%
    slogData.sig3d = cat(4,ImgPerExy,ImgPerEyz,ImgPerExz);
    slogData.pha3d = cat(4,PhasePerExy,PhasePerEyz,PhasePerExz);
end
end

function [s]=DFT(slogData,pos,resseq,NKlength) %for KmapSignal
global TDImage;
if slogData.KmapSignal == 1
    kvec = slogData.kvector(1+NKlength(1):NKlength(2)-1,:);
    sk = slogData.signal(1+NKlength(1):NKlength(2)-1,1)+1i*slogData.signal(1+NKlength(1):NKlength(2)-1,2);
elseif slogData.KmapSignal == 2
    kvec = slogData.kvector(1+NKlength(1):NKlength(2),:);
    sk = slogData.signal(1+NKlength(1):NKlength(2),1)+1i*slogData.signal(1+NKlength(1):NKlength(2),2);
end
if isempty(TDImage) || TDImage
    kfilter=2*pi/resseq/2;
    kvecF = kvec((abs(kvec(:,1))<=kfilter)&(abs(kvec(:,2))<=kfilter)&(abs(kvec(:,3))<=kfilter),:);
    skF = sk((abs(kvec(:,1))<=kfilter)&(abs(kvec(:,2))<=kfilter)&(abs(kvec(:,3))<=kfilter),:);    
    s = skF'*exp(-1i*kvecF*pos);
else
    s = sk'*exp(-1i*kvec*pos);
end
end
