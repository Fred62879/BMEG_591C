function tile_correction()
    numLines = size(OCT_mcorr, 2);
    msfac = 1;

    for I = 1:numines
        [output,~] = dftregistration( fft2(squeeze( OCT_mcorr(:,round(numLines./2),:) )), ...
                                      fft2(squeeze(OCT_mcorr(:,I,:))), usfac );

        axialshift(I) = round(output(3));

        OCT_mcorr(:,I,:) = circshift(OCT_mcorr(:,I,:), [output(3),0]);
    end
end