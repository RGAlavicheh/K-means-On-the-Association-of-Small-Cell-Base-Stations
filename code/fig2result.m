% Author: Reza Ghasemi Alavicheh [r.ghasemi.reza@gmail.com]

clc;clear
%% parameters
iter =300;
L =1:1:17;
R_B = 1400;
W = 200;
power = 1.3;
cluster = 4;
[X,~] = materni(cluster,"off");
if length(X) >29 || length(X)<27
    clc;
    clear;
    fig2result
    return
end
pd = makedist('Nakagami','mu',1,'omega',1);
h = random(pd,1,2);
%%
for conventional = 1:1:2
    
    data_a = zeros(length(L),1);
    data_b = zeros(length(L),1);
    if conventional == 1
%         [~,C] = materni(cluster,"off"); %kmeans in matern function
        LW = 1;
        [idx,C] = kmeans(X,cluster);
        C_kmeans = C;
    else
        [s_min] = s_computer(X);
        [C] = test(s_min);
        C = C(1:4,:);
        LW =1;
        C_conventional = C;
    end

    scenario =[1 2];
    for i = 1:1:length(L)
        [probOFasso,average_num,sum_rate1] = func_s1(iter,L(i),R_B,W,power,X,h,C,cluster);
        data_a(i,1) = average_num;
        data_b(i,1) = probOFasso;  

    end

    LS = ["-" "--"];
    if conventional == 1
        figure
    else
        color_v = ['b','m','c'];
        marker_v = ['.','.','.'];

    end
    if conventional == 2
        color_p = color_v(1);
        marker_p = marker_v(1);
    else
        color_p = 'r';
        marker_p = 'x' ;
    end
    plot(L,data_a,marker_p,color=color_p,Marker=marker_p,LineWidth=LW,MarkerEdgeColor='k',LineStyle=LS(conventional));
    ylim([1 30]);
    xlim([1 15]);
    grid on;
    hold on;

    for i = 1:1:length(L)
        [probOFasso,average_num,sum_rate2] = func_s2(iter,L(i),R_B,W,power,X,h,C,cluster);
        data_a(i,1) = average_num;
        data_b(i,1) = probOFasso;

    end
    if conventional == 2
        color_p = color_v(2);
        marker_p = marker_v(2);

    else
        color_p = 'k';
        marker_p = 'v' ;
    end
    plot(L,data_a,Color=color_p,marker=marker_p ,MarkerEdgeColor='k',LineWidth=LW,LineStyle=LS(conventional));
    ylim([1 30]);
    xlim([1 15]);
    grid on;
    hold on;




    for i = 1:1:length(L)
        [probOFasso,average_num,sum_rate3] = func_s3(iter,L(i),R_B,W,power,X,h,C,cluster);

        data_a(i,1) = average_num;
        data_b(i,1) = probOFasso;

    end
    if conventional == 2
        color_p = color_v(3);
        marker_p = marker_v(3);

    else
        color_p = 'g';
        marker_p = '+';

    end
    plot(L,data_a,Color=color_p,marker=marker_p,MarkerEdgeColor='k',LineWidth=LW,LineStyle=LS(conventional));
    ylim([1 30]);
    xlim([1 15]);
    grid on;
    hold on;
    legend("S1 k-means","S2 k-means","S3 k-means", "S1 conventional","S2 conventional","S3 conventional");
    xlabel('L (number of links)','FontSize',12)
    ylabel("Avg. number of association")

end
compare_func(iter,11,R_B,W,power,cluster,X,C_kmeans,C_conventional,h)
