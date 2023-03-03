function avgImg = mov2DAvg(inputImg, windXY)

kernel = ones([windXY(1), windXY(2)])./(windXY(1)*windXY(2));
avgImg = imfilter(inputImg, kernel, 'same');

end
