function OCT_mcorr = motion_correction(OCT_Data)
    % Depth_ROI = [A B];
    Depth_ROI = [1 260];
    numFrames = size(OCT_Data, 3);
    OCT_ROI = OCT_Data(Depth_ROI(1):Depth_ROI(2),:,:);

    usfac = 1; % upsampling factor

    axialShift = zeros(numFrames, 1);

    OCT_mcorr = OCT_ROI;

    for I = 1:numFrames
        [output, ~] = dftregistration(...
            fft2( 20 .* log10(abs(OCT_ROI(:,:,round(numFrames)./2)))),...
            fft2(20 .* log10(abs(OCT_ROI(:,:,I)))),...
            usfac);

        axialShift(I) = round(output(3));

        OCT_mcorr(:,:,I) = circshift(OCT_ROI(:,:,I), [output(3), 0]);
    end
end
