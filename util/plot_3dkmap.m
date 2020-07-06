% function plot_3dkmap(Kgrid3D,Bgrid3D,timestamp)
% in:                   Kgrid3D; 3D(K-vector grids) + 1D(time)
%                       Bgrid3D; 3D(Background grids) + 1D(time)
%                       timestamp; Nth image(optimal)
% out:
%
% author: xiang gao
% xiang.gao@uniklinik-freiburg.de

function plot_3dkmap(Kgrid3D,Bgrid3D,timestamp)
global Kmax2Pi;
global Mat;
global TDImage;
global plot3proj;
if isempty(plot3proj)
    plot3proj=true;
end
if isempty(Bgrid3D)
    TDImage = false;
end
if nargin<3
    timestamp=[];
end
values=[];
coords=[];
count=0;
for i=1:size(Kgrid3D,1)
    for j=1:size(Kgrid3D,2)
        for k=1:size(Kgrid3D,3)
            if Kgrid3D(i,j,k)>0
                count=count+1;
                values(count)=Kgrid3D(i,j,k);
                coords=cat(1,coords,[i j k]);
            end
        end
    end
end

if TDImage
    valuesB=[];
    coordsB=[];
    count=0;
    for i=1:size(Bgrid3D,1)
        for j=1:size(Bgrid3D,2)
            for k=1:size(Bgrid3D,3)
                if Bgrid3D(i,j,k)>0
                    count=count+1;
                    valuesB(count)=Bgrid3D(i,j,k);
                    coordsB=cat(1,coordsB,[i j k]);
                end
            end
        end
    end
    
    coordsInBox = [];
    valuesInBox = [];
    count=0;
    for i=1:size(coords,1)
        if (all(coordsB(1,:)<coords(i,:) & coords(i,:)<coordsB(end,:)))
            % we don't count the kvecotrs who sits on the gird(to be consistant
            % with image domain who has higher precision than k-gird)
            count=count+1;
            valuesInBox(count)=values(i);
            coordsInBox=cat(1,coordsInBox,coords(i,:));
        end
    end
    
    % figure;
    if Mat==1
        scatter3(coordsB(:,1),coordsB(:,2),coordsB(:,3),[],valuesB','ko');
    else
        scatter3(coordsB(:,1),coordsB(:,2),coordsB(:,3),60,valuesB','k.');
    end
    hold on
    if ~isempty(coords)
        scatter3(coords(:,1),coords(:,2),coords(:,3),[],values','*');
    end
    if coordsInBox & plot3proj
        hold on
        scatter3(coordsInBox(:,1),coordsInBox(:,2),0*coordsInBox(:,3),[],valuesInBox','s');
        hold on
        scatter3(0*coordsInBox(:,1),coordsInBox(:,2),coordsInBox(:,3),[],valuesInBox','^');
        hold on
        scatter3(coordsInBox(:,1),size(Kgrid3D,2)+0*coordsInBox(:,2),coordsInBox(:,3),[],valuesInBox','p');
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
    colormap(winter);
    colorbar;
    zlabel('K_z [1/m]');set(gca,'FontSize',12);
    view(60,30)
    xlabel('K_x');ylabel('K_y');
    xticks(0:(size(Bgrid3D,1)+1)/2:size(Bgrid3D,1)+1)
    yticks(0:(size(Bgrid3D,2)+1)/2:size(Bgrid3D,2)+1)
    zticks(0:(size(Bgrid3D,3)+1)/2:size(Bgrid3D,3)+1)
    xticklabels(Ticklabel(Kmax2Pi))
    yticklabels(Ticklabel(Kmax2Pi))
    zticklabels(Ticklabel(Kmax2Pi))
    axis([0,size(Bgrid3D,1)+1,0,size(Bgrid3D,2)+1,0,size(Bgrid3D,3)+1])
else
    if ~isempty(coords)
        scatter3(coords(:,1),coords(:,2),coords(:,3),[],values','*');
    end
    colormap(winter);
    colorbar;zlabel('k [1/m]')
    view(90,0)
end
end

function string=Ticklabel(factor)

% string = [{mat2str(-1000*factor)},{mat2str(-500*factor)},...
%     {mat2str(0*factor)},{mat2str(500*factor)},{mat2str(1000*factor)}];
string = [{mat2str(-1000*factor)},{mat2str(0*factor)},{mat2str(1000*factor)}];
end