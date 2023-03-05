function [xShift_axial, yShift_axial, cplxOCT_mcorr_local] = local_motion_correction(...
    procd_OCT_BM_ROI, xShift_axial, yShift_axial, estimate_shift)

    cplxOCT_ROI = procd_OCT_BM_ROI;
    numFrames = size(cplxOCT_ROI, 3);
    %Depth_ROI = [40 399];
    BM = 2;
    usfac = 1;

    %cplxOCT_ROI = cplxOCT(Depth_ROI(1):Depth_ROI(2),:,:);
    cplxOCT_mcorr_local = cplxOCT_ROI;

    for I = 1:BM:numFrames
        if estimate_shift
            [output, ~] = dftregistration(fft2(20.*log10(abs(cplxOCT_ROI(:,:,I)))),...
                                          fft2(20.*log10(abs(cplxOCT_ROI(:,:,I+1)))), usfac);
            xShift_axial(I) = round(output(4));
            yShift_axial(I) = round(output(3));
        end

        % cplxOCT_mcorr_local(:,:,I+1) = circshift(cplxOCT_ROI(:,:,I+1),[output(3), output(4)]);
        cplxOCT_mcorr_local(:,:,I+1) = circshift(cplxOCT_ROI(:,:,I+1),...
            [yShift_axial(I), xShift_axial(I)]);
    end

end