function Clustering_iterInfoAnalyzer(iterInfo,ifPlot,ifSave,saveFileName,saveDir)
%CLUSTERING_ITERINFOANALYZER 此处显示有关此函数的摘要
%   此处显示详细说明


convergeIter = iterInfo.convergeIter;
iterData = iterInfo.iterData;
% iterData(1,1) = 0;

colorMap = myColorMap("BR1");
color1 = colorMap(ceil(size(colorMap,1)*0.25),:);
color2 = colorMap(ceil(size(colorMap,1)*0.75),:);

if convergeIter
    iterData = iterData(1:convergeIter,:);
end


if ifPlot || ifSave
    f = figure(Units="inches");
    f.Position(3:4) = [3.5 2.5];

    x = 1 : size(iterData,1);
    hold on
    yyaxis left
    plot(x,iterData(:,1),"Marker",".","LineStyle","-","Color",color1)
    ylabel("Inter Iteration Similarity",VerticalAlignment="baseline")
    set(gca,'ycolor',color1);
    plot(convergeIter,iterData(convergeIter,1),"Marker","o","LineStyle","-","Color",color1)

    yyaxis right
    plot(x,iterData(:,2),"Marker","x","LineStyle",":","Color",color2)
    ylabel("Cluster Number",Rotation=-90,VerticalAlignment="baseline")
    set(gca,'ycolor',color2);

    hold off


    grid on
    xlabel("Iteration",VerticalAlignment="top")

    set(gca,'tickdir','out')
    totalIter = length(iterData(:,1));
    xlim([1-0.025*totalIter, 1.025 * totalIter])
    box on
    if ifSave
        %exportgraphics(f, saveDir + "iterInfo_" + saveFileName + ".png", "Resolution","500")
        exportgraphics(f, saveDir + "iterInfo_" + saveFileName + ".svg", "ContentType","vector");
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = figure(Units="inches");
    f.Position(3:4) = [3.5 2.5];
    hold on
    plot(x,iterData(:,3),"Marker","diamond","LineStyle","-","Color",'k')
    ylabel("Outlier Count", VerticalAlignment="baseline")
    set(gca,'ycolor','k');
    plot(convergeIter,iterData(convergeIter,3),"Marker","o","LineStyle","-","Color",'k')
    hold off
    grid on
    xlabel("Iteration",VerticalAlignment="top")

    set(gca,'tickdir','out')
    totalIter = length(iterData(:,1));
    xlim([1-0.025*totalIter, 1.025 * totalIter])
    box on
    if ifSave
        %exportgraphics(f, saveDir + "iterInfo_outlierCount_" + saveFileName + ".png", "Resolution","500")
        exportgraphics(f, saveDir + "iterInfo_outlierCount_" + saveFileName + ".svg", "ContentType","vector");
    end

end


end

