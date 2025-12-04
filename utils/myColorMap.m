function colorMap = myColorMap(colorMapName)
%MYCOLORMAP 此处显示有关此函数的摘要
%   此处显示详细说明
switch colorMapName
    case "BR1"
        colorMap = [
            51  51 255;
            51 102 255;
            102 153 255;
            153 204 255;
            204 236 255;
            255 255 255;
            255 204 204;
            255 153 153;
            255 124 128;
            255  80  80;
            204   0   0;
            ]/255;
    case "BR2"
        colorMap = [
            '#00BBFF'
            '#0EB1F1'
            '#1CA6E3'
            '#2B9CD5'
            '#3991C6'
            '#4787B8'
            '#557DAA'
            '#63729C'
            '#71688E'
            '#805E80'
            '#8E5371'
            '#9C4963'
            '#AA3E55'
            '#B83447'
            '#C62A39'
            '#D41F2B'
            '#E3151C'
            '#F10A0E'
            '#FF0000'
            ];
end







end


