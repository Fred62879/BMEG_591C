function volume_correction()
    % fast b scan
    for i = 1:size(prodc_OCT_ROI, 3)
        plot()
        pause(0.1)
    end

    % slow b scan (very bumpy)
    plot( abs(squeeze(scanROI(:,250,:))) )

    % motion correction
    yShift_axial


    figure()
    plot(yShift_axisl_tile())


    % fit a curve to the data
    x = [1:numLines]';
    cx = polyFit();

    % apply the fitted tilt estimation

    % tilde correction
end