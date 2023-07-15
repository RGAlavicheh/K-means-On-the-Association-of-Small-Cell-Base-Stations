% Simulate Matern hard-core point processes (Type I).
% Author: H. Paul Keeler, 2019.
% Website: hpaulkeeler.com
% Repository: github.com/hpaulkeeler/posts
% For more details, see the post:
% hpaulkeeler.com/simulating-matern-hard-core-point-processes

%%%%%%%%%%%%%%%%%%%%  about this code %%%%%%%%%%%%%%%
% I have edited the original code to align with my specific objectives.(Reza Ghasemi Alavicheh)[r.ghasemi.reza@gmail.com]
% This code implements both the Matérn I pattern generation and the k-means clustering algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X,C] = materni(cluster,figure_status)
xMin=0;
xMax=4;
yMin=0;
yMax=4;


lambda=2;%density 
r=0.25;%radius |0.25km = 250m

%Extended windows
rExt=r; %extension parameter -- use core radius
xMinExt=xMin-rExt;
xMaxExt=xMax+rExt;
yMinExt=yMin-rExt;
yMaxExt=yMax+rExt;
%dimensions
xDeltaExt=xMaxExt-xMinExt;
yDeltaExt=yMaxExt-yMinExt;
area_Ext=xDeltaExt*yDeltaExt; %extended area

%generate Poisson point 
numbPointsExt=poissrnd(area_Ext*lambda,1,1);%Poisson number
%x and y coordinates 
xxPoissonExt=xMinExt+xDeltaExt*rand(numbPointsExt,1);
yyPoissonExt=yMinExt+yDeltaExt*rand(numbPointsExt,1);

%remove points if outside the window
booleWindow=((xxPoissonExt>=xMin)&(xxPoissonExt<=xMax)&(yyPoissonExt>=yMin)&(yyPoissonExt<=yMax));

%retain points inside simulation window
xxPoisson=xxPoissonExt(booleWindow);
yyPoisson=yyPoissonExt(booleWindow);

numbPoints=length(xxPoisson); %number of Poisson points in window
%create random marks for ages

% START Removing
%% false is logical zeros
booleRemoveI=false(numbPoints,1);%Index for removing points -- Matern I

for ii=1:numbPoints
    distTemp=hypot(xxPoisson(ii)-xxPoissonExt,...
        yyPoisson(ii)-yyPoissonExt);  %distances to other points

    booleInDisk=(distTemp<r)&(distTemp>0);%check if inside disk

    %Matern I
    booleRemoveI(ii)=any(booleInDisk);

end

%Matérn I
booleKeepI=~(booleRemoveI);
xxMaternI=xxPoisson(booleKeepI);
yyMaternI=yyPoisson(booleKeepI);


if figure_status == "on"
    %% Plotting
    markerSize=30; %marker size of markers colors
    figure
    plot(xxPoisson,yyPoisson,'ko','MarkerSize',markerSize/4);
    hold on;
    plot(xxMaternI,yyMaternI,'rx','MarkerSize',markerSize/2);
    legend('PPP','PPP (Matern I)');
    xlabel 'X(Km)';
    ylabel 'Y(Km)';
    grid on;
    hold on;
end
%% %%%%%%%%%%% k-means %%%%%%%%%%

X = zeros(length(xxMaternI),2);
X(:,1) = xxMaternI ;
X(:,2) = yyMaternI ;

opts = statset('Display','final');
[idx,C] = kmeans(X,cluster,'Replicates',5,'Options',opts);

if figure_status == "on"
    %% 2D plot for 4 clusters
    markerSize=30; %marker size of markers colors
    figure
    plot(xxPoisson,yyPoisson,'ko','MarkerSize',markerSize/4);
    hold on;
    plot(xxMaternI,yyMaternI,'rx','MarkerSize',markerSize/2);
    plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',16)
    hold on
    plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',16)
    plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',16)
    plot(X(idx==4,1),X(idx==4,2),'c.','MarkerSize',16)
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
    legend('PPP (removed points)','PPP (Matern i)','cluster 1','cluster 2','cluster 3','cluster 4',...
        'UAVs')
    title 'Cluster Assignments and Centroids'
    xlabel 'X(Km)';
    ylabel 'Y(Km)';
    hold off;
    grid on;

    %% 3D plot for 4 clusters
    figure
    temp1 = X(idx==1,1);temp2 = X(idx==2,1);temp3 = X(idx==3,1);temp4 = X(idx==4,1);

    plot3(X(idx==1,1),X(idx==1,2),zeros(length(temp1),1),'r^','MarkerSize',10,'MarkerFaceColor',"r","MarkerEdgeColor","k")
    hold on
    plot3(C(1,1),C(1,2),0.3,'xr','MarkerSize',15,'LineWidth',3)

    plot3(X(idx==2,1),X(idx==2,2),zeros(length(temp2),1),'b^','MarkerSize',10,'MarkerFaceColor',"b","MarkerEdgeColor","k")
    plot3(C(2,1),C(2,2),0.3,'xb','MarkerSize',15,'LineWidth',3)

    plot3(X(idx==3,1),X(idx==3,2),zeros(length(temp3),1),'g^','MarkerSize',10,'MarkerFaceColor',"g","MarkerEdgeColor","k")
    plot3(C(3,1),C(3,2),0.3,'xg','MarkerSize',15,'LineWidth',3)

    plot3(X(idx==4,1),X(idx==4,2),zeros(length(temp4),1),'c^','MarkerSize',10,'MarkerFaceColor',"c","MarkerEdgeColor","k")
    plot3(C(4,1),C(4,2),0.3,'xc','MarkerSize',15,'LineWidth',3)
    legend('SCBSs of cluster1','UAV 1','SCBSs of cluster2','UAV 2','SCBSs of cluster3','UAV 3','SCBSs of cluster4','UAV 4')
    title 'Cluster Assignments and Centroids'
    xlabel 'X(Km)';
    ylabel 'Y(Km)';
    zlabel 'Z(Km)';
    hold off;
    grid on;
    %% % another figure
    x1 = min(X(:,1)):0.01:max(X(:,1));
    x2 = min(X(:,2)):0.01:max(X(:,2));
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot
    idx2Region = kmeans(XGrid,4,'MaxIter',1,'Start',C);
    figure;
    gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
        [0,0.32,0.75;0.75,0,0.75;0.75,0.75,0;0,0.4,0.15],'..');
    hold on;
    plot(X(:,1),X(:,2),'k*','MarkerSize',5);
    title 'coverage region of each UAV';
    xlabel 'X(Km)';
    ylabel 'Y(Km)';
    legend('Region 1','Region 2','Region 3','Region 3','SCBSs','Location','SouthEast');
    hold off;
end
end
