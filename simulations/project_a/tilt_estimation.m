function axialshift = tilt_estimation(OCT_mcorr)
    numLines = size(OCT_mcorr, 2);
    usfac = 1;

    for I = 1:numLines
        [output,~] = dftregistration( fft2(squeeze(OCT_mcorr(:,round(numLines./2),:) )), ...
                                      fft2(squeeze(OCT_mcorr(:,I,:))), usfac );

        axialshift(I) = round(output(3));
    end
end