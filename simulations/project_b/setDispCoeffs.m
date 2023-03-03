function dispCoeffs = setDispCoeffs(rawData, depthROI, maxDispOrders, coeffRange)
    dispCoeffs = zeros(1, maxDispOrders-1);
    for i = 1:length(dispCoeffs)
        dispCoeffRng = [coeffRange, -1 * coeffRange];
        costs = zeros(size(dispCoeffRng));

        for j = 1:length(dispCoeffRng)
            dispCoeffs(i) = dispCoeffRng(j);
            costs(j) = calCostFun(rawData, depthROI, maxDispOrders, dispCoeffs);
        end

        for k = 1:20
            [~,idx] = sort(costs);
            dispCoeff_min1 = dispCoeffRng(idx(1));
            dispCoeff_min2 = dispCoeffRng(idx(2));
            dispCoeff_new = (dispCoeff_min1 + dispCoeff_min2)/2;
            dispCoeffRng = [dispCoeffRng, dispCoeff_new];
            dispCoeffs(i) = dispCoeff_new;
            cost_new = calCostFun(rawData, depthROI, maxDispOrders, dispCoeffs);
            costs = [costs, cost_new];
        end

        [~,argmin] = min(costs);
        dispCoeffs(i) = dispCoeffRng(argmin);
    end
end

function cost = calCostFun(rawData, depthROI, maxDispOrders, dispCoeffs)
    rawData_dspComp = compDisPhase(rawData, maxDispOrders, dispCoeffs);

    OCT = (abs(fft(rawData_dspComp))).^2;
    OCT_ROI = OCT(depthROI(1) : depthROI(2),:);
    OCT_norm = OCT_ROI ./ sum(OCT_ROI(:));
    entropy = -1*((OCT_norm) .* log(OCT_norm));
    cost = sum(entropy(:));
end