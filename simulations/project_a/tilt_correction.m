function OCT_tcorr = tilt_correction(OCT_mcorr, axialshift)
    numLines = size(OCT_mcorr, 2);
    disp(axialshift(1,1));

    for I = 1:numLines
        OCT_tcorr(:,I,:) = circshift(OCT_mcorr(:,I,:), [round(axialshift(1,I)), 0]);
    end
end