% load matrix3d.mat

function plot_3dkmap(Kgrid3D,bgmatrix,timestamp)
global Kmax2Pi;
global Mat;
global TDImage;
global plot3proj;
if isempty(plot3proj)
    plot3proj=true;
end
if isempty(bgmatrix)
    TDImage = false;
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
    for i=1:size(bgmatrix,1)
        for j=1:size(bgmatrix,2)
            for k=1:size(bgmatrix,3)
                if bgmatrix(i,j,k)>0
                    count=count+1;
                    valuesB(count)=bgmatrix(i,j,k);
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
            % we don't count the kvecotrs who sits on the gird( to be consistant
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
        scatter3(coordsB(:,1),coordsB(:,2),coordsB(:,3),[],valuesB','k.');
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
    if isnan(timestamp)
        legend('Receive Zone','Whole K-vectors','Proj. XY','Proj. YZ','Proj. XZ','Location','northwest');
        set(legend,...
            'Position',[0.0016174430593756 0.478469762050102 0.142276420068692 0.121468923186178]);
    else
        title(['3D K-Space Time ',int2str(timestamp)]);
        if timestamp==3
            legend('Receive Zone','Whole K-vectors','Proj. XY','Proj. YZ','Proj. XZ','Location','northwest');
            set(legend,...
                'Position',[-0.0186528272108947 0.43782627308726 0.142276420068692 0.159448814497689]);
                %'Position',[0.0016174430593756 0.478469762050102 0.142276420068692 0.121468923186178]);
        end
    end
    colormap(winter);
    colorbar;zlabel('k [1/m]')
    view(60,30)
    xticks(0:(size(bgmatrix,1)+1)/4:size(bgmatrix,1)+1)
    yticks(0:(size(bgmatrix,2)+1)/4:size(bgmatrix,2)+1)
    zticks(0:(size(bgmatrix,3)+1)/4:size(bgmatrix,3)+1)
    xticklabels(Ticklabel(Kmax2Pi))
    yticklabels(Ticklabel(Kmax2Pi))
    zticklabels(Ticklabel(Kmax2Pi))
    axis([0,size(bgmatrix,1)+1,0,size(bgmatrix,2)+1,0,size(bgmatrix,3)+1])
else
    %     figure;
    if ~isempty(coords)
        scatter3(coords(:,1),coords(:,2),coords(:,3),[],values','*');
    end
    colormap(winter);
    colorbar;zlabel('k [1/m]')
    view(90,0)
end
end

function string=Ticklabel(factor)

string = [{mat2str(-1000*factor)},{mat2str(-500*factor)},...
    {mat2str(0*factor)},{mat2str(500*factor)},{mat2str(1000*factor)}];

end