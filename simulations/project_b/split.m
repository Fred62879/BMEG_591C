function RawData_Split = split(RawData)
% repeat project a processing w/o hanning windowing

    % split spectra window
    Fraction = 2;
    numPoints = size(RawData, 1);
    numBands = 2*Fraction - 1;
    winWidth = round(numPoints/Fraction);
    scale = sqrt(sum(hanning(numPoints).^2)/sum(hanning(winWidth).^2));

    for I = 1:numBands
        wins(:,I) = circshift(cat(1,hanning(winWidth),...
                                  zeros(numPoints-winWidth,1)),...
                              round((I-1)*numPoints/(numBands+1)))*scale;
    end
    
    splitWins = permute(shiftdim(wins,-1),[2,1,3]);
    RawData_Split = bsxfun(@times, RawData, splitWins);
end