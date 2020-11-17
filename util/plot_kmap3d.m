% function kmap3d(slogData,nth,timestamp)
% in:                   slogData; output from KmapSignal2 function
%                            nth; nth echo
%                      Kmax4plot; Kmap boundary (1/m)
% out:
%
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

function kmap3d(slogData,nth,Kmax4plot,timestamp)
global Mat;
global FOV;
plot3proj=true;
if nargin<4
    timestamp=[];
end

NKperEcho = slogData.NKperEcho;
kvector =slogData.kvector./2/pi*1e6;
signal = abs(slogData.signal(:,1)+1i*slogData.signal(:,2));

k=kvector(NKperEcho(nth)+1:NKperEcho(nth+1),:);
s=signal(NKperEcho(nth)+1:NKperEcho(nth+1));

%k sampling area %maybe a more efficent way?
boundary=[-1,0,1];cnt=0;
ksa=zeros(9*length(boundary),length(boundary));
for a=1:3
    for b=1:3
        for c=1:3
            cnt=cnt+1;
            ksa(cnt,:)=Mat/FOV/2*1e3.*[boundary(a),boundary(b),boundary(c)];
        end
    end
end

if ~isempty(k)
    scatter3(k(:,1),k(:,2),k(:,3),[],s,'*');hold on
end
if plot3proj
    kc=[];sc=[];
    for i=1:size(k,1)
        if (all(abs(k(i,:))<=Mat/FOV/2*1e3))
            kc=cat(1,kc,k(i,:));
            sc=cat(1,sc,s(i));
        end
    end
    if ~isempty(kc)
        scatter3(kc(:,1),kc(:,2),-Kmax4plot+0*kc(:,3),[],sc,'s');hold on
        scatter3(-Kmax4plot+0*kc(:,1),kc(:,2),kc(:,3),[],sc,'^');hold on
        scatter3(kc(:,1),Kmax4plot+0*kc(:,2),kc(:,3),[],sc,'p');hold on
    end
end
if Mat==1
    scatter3(ksa(:,1),ksa(:,2),ksa(:,3),[],ones(1,size(ksa,1)),'ko');
else
    scatter3(ksa(:,1),ksa(:,2),ksa(:,3),60,ones(1,size(ksa,1)),'k.');
end
hold off

if ~plot3proj
    legend('Sampled area','3D k-vectors');
    set(legend,...
        'Position',[0.00630073579088022 0.535565113347637 0.143292680201007 0.104422601640078]);
end
if plot3proj&& ~isempty(timestamp)
    title(['Time ',int2str(timestamp)]);
    if timestamp==3
        legend('Sampled area','3D k-vectors','Proj. XY','Proj. YZ','Proj. XZ','Location','northwest');
        set(legend,...
            'Position',[0.447810587423252 0.423250434344067 0.142276420068692 0.171627366100058]);
    end
end
axis(Kmax4plot.*[-1,1,-1,1,-1,1]);
colormap(winter);
colorbar;%caxis([0,0.5]);
set(gca,'FontSize',12);
xlabel('K_x');ylabel('K_y');zlabel('K_z [1/m]');
view(60,30)
end
