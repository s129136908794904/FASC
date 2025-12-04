function colorMap = sanitizeBR1ColorMap(colorMap)
%SANITIZEBR1COLORMAP remove extremes/white from BR1 palette
if size(colorMap,1) > 2
    colorMap = colorMap(2:end-1,:);
end
whiteRows = all(colorMap >= 0.999,2);
colorMap = colorMap(~whiteRows,:);
if isempty(colorMap)
    colorMap = myColorMap("BR1");
end
end
