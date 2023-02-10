function rawData_dispComp = compDisPhase(rawData, maxDispOrders, dispCoeffs)
    ScanPts = size(rawData, 1);
    LinePerFrame = size(rawData, 2);
    klinear = linspace(-1,1,ScanPts);
    kaxis = repmat(klinear',1,LinePerFrame);

    rawData_dispComp = rawData;

    for i = 1:maxDispOrders - 1
        rawData_dispComp = rawData_dispComp .* exp(1j .* (dispCoeffs(i) * (kaxis .^ (i+1))));
    end
end