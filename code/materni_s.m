% Simulate Matern hard-core point processes (Type I).
% Author: H. Paul Keeler, 2019.
% Website: hpaulkeeler.com
% Repository: github.com/hpaulkeeler/posts
% For more details, see the post:
% hpaulkeeler.com/simulating-matern-hard-core-point-processes/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I have edited the original code to align with my specific objectives.(Reza Ghasemi Alavicheh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = materni_s(s_min)
numbSim=1000; %number of simulations default 1000
%Simulation window parameters
xMin=0;
xMax=4;
yMin=0;
yMax=4;

%Parameters for the parent and daughter point processes
lambdaPoisson=2;%density of parent Poisson point process default = 2
radiusCore=(s_min)/1000;%radius of hard core |0.25km = 250m

%Extended simulation windows parameters
rExt=radiusCore; %extension parameter -- use core radius
xMinExt=xMin-rExt;
xMaxExt=xMax+rExt;
yMinExt=yMin-rExt;
yMaxExt=yMax+rExt;
%rectangle dimensions
xDeltaExt=xMaxExt-xMinExt;
yDeltaExt=yMaxExt-yMinExt;
areaTotalExt=xDeltaExt*yDeltaExt; %area of extended rectangle

%%%START Simulations START%%%%
%initialize arrays for collecting statistics
numbPointsAll=zeros(numbSim,1); %number of Poisson points
numbPointsAll_I=zeros(numbSim,1); %number of Matern I points
%numbPointsAll_II=zeros(numbSim,1); %number of Matern II points
%loop through for each simulation
for ss=1:numbSim
    
    %Simulate Poisson point process for the parents
    numbPointsExt=poissrnd(areaTotalExt*lambdaPoisson,1,1);%Poisson number
    %x and y coordinates of Poisson points for the parent
    xxPoissonExt=xMinExt+xDeltaExt*rand(numbPointsExt,1);
    yyPoissonExt=yMinExt+yDeltaExt*rand(numbPointsExt,1);
    
    %thin points if outside the simulation window
    booleWindow=((xxPoissonExt>=xMin)&(xxPoissonExt<=xMax)...
        &(yyPoissonExt>=yMin)&(yyPoissonExt<=yMax));
    indexWindow=find(booleWindow);
    %retain points inside simulation window
    xxPoisson=xxPoissonExt(booleWindow);
    yyPoisson=yyPoissonExt(booleWindow);
    
    numbPoints=length(xxPoisson); %number of Poisson points in window
    %create random marks for ages
    markAge=rand(numbPointsExt,1);
    
    %%%START Removing/thinning points START%%%
    booleRemoveI=false(numbPoints,1);%Index for removing points -- Matern I
%     booleKeepII=false(numbPoints,1);%Index for keeping points -- Matern II
    for ii=1:numbPoints
        distTemp=hypot(xxPoisson(ii)-xxPoissonExt,...
            yyPoisson(ii)-yyPoissonExt);  %distances to other points
        
        booleInDisk=(distTemp<radiusCore)&(distTemp>0);%check if inside disk
        
        %Matern I
        booleRemoveI(ii)=any(booleInDisk);
        
        %Matern II
        %keep the younger points
%         booleKeepII(ii)=all(markAge(indexWindow(ii))<markAge(booleInDisk));
        %Note: if markAge(booleInDisk) is empty, keepBooleII(ii)=1.
    end
    %%%END Removing/thinning points END%%%
    
    %Remove/keep points to generate Matern hard-core processes
    %Matérn I
    booleKeepI=~(booleRemoveI);
    xxMaternI=xxPoisson(booleKeepI);
    yyMaternI=yyPoisson(booleKeepI);
    

    
    %Update statistics
    numbPointsAll(ss)=numbPoints;
    numbPointsAll_I(ss)=length(xxMaternI);

end
%%%END Simulations END%%%%

%Compare statistics
diskArea=pi*radiusCore^2; %area of disk
areaWindow=(xMax-xMin)*(yMax-yMin); %area of simulation window

%underlying Poisson point process
lambdaExact=lambdaPoisson;
lambdaEmp=mean(numbPointsAll/areaWindow);
X = zeros(length(xxMaternI),2);
X(:,1) = xxMaternI ; 
X(:,2) = yyMaternI ; 
C = X;

end
