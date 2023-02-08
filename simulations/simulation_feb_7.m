

function simulation_feb_7

    function rawData_rescaled = reSampling_CalSig(fftData, cal_RawData)
        ifftData = ifft(fftData);
        rawData_real = real(ifftData);
        rawData_imag = image(ifftData);
        rawData_rescaled = zeros(size(fftData));

        idx = (1:length(cal_RawData))';

        for i = 1:size(fftData,2)
            phase = unwrap(angle(cal_RawData(:,i)));
            phase_norm = phase - min(phase(:));
            phase_norm = phase_norm ./ max(phase_norm(:));

            p = fit(phase_norm, idx, 'cubicinterp');

            phase_lin = linspace(0, 1, length(cal_RawData))';

            idx_linK = p(phase_lin);

            rawData_rescaled(:,i) = (interp1(...
                [1:size(rawData_real,1)]', rawData_real(:,i), idx_linK,'spline'))...
                + 1j.*(interp1([1:size(rawData_image,1)]',rawData_imag(:,i),idx_linK,'spline'))
        end
    end

    function rawData_PhaseComp = compPhaseShift(ref_fftData_1D, fftData_2D, depthIdx)
        ref_Ascan = ref_fftData_1D;

        cplxConjX = fftData_2D...
            .* repmat(conj(ref_Ascan), [1 size(fftData_2D,2)]);

        calSigDepth = depthIdx;

        phaseSlope = (angle(cplxConjX(calSigDepth,L)) /calSigDepth)...
            .* linspace(1,length(ref_Ascan), length(ref_Ascan))';

        for i = 1:size(fftData_2D, 2)
            fftData_PhaseComp(:,i) = fftData_2D(:,i)...
                .* exp(-1j * phaseSlope(:,i));
        end

        rawData_PhaseComp = ifft(fftData_PhaseComp);
    end


    function rawData_FPNSub = FPNSubtracton(rawData)
        rawData_FPNSub = rawData...
            - (repmat(median(real(rawData, 2)), [1,size(rawData,2)]) ...
            +1j .* repmat(median(imag(rawData), 2), [1, size(rawData, 2)]));
    end

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

    function motion_correction()
        Depth_ROI = [A B];
        numFrames = size(OCT_Data, 3);
        OCT_ROI = OCT_Data(Depth_ROI(1):Depth_ROI(2),:,:);

        usfac = 1;

        axialShift = zeros(numFrames, 1);

        OCT_mcorr = OCT_ROIl
        for I = 1LnumFrames
            [output, ~] = dftregistration(fft2( 20 .* log10(abs(OCT_ROI(:,:,round(numFrames)./2))) )),...
                fft2(20 .* log10(abs(OCT_ROI(:,:,I))), usfac);

            axialShift = round(output(3));

            OCT_mcorr(:,:,I) = circshift(OCT_ROI(:,:,I), [output(3), 0]);
        end
    end

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

    function volume_processing(rawOCT)
        for frame_num = 1:size(rawOCT, 3)
            raw_data = rawOCT(:,:,frame_num);

            procd_data(:,:,frame_num) = FFTData(1:size(raw_data, 1)/2,:);
        end
    end

    fn = '/media/fred/working_drive/2022W2/bmeg591c/simulations/data/RawOCT.mat';
    %load(fn)
    raw_data = load(fn);
    rawOCT = getfield(raw_data,'rawOCT'); %_array');

    ref_rawData = rawOCT(:,:,250);
    ref_FFTData = fft(hilbert(ref_rawData));

    figure(); plot(abs(ref_FFTData(1:end/2,:))); ylim([0 4e6])
%     figure(); plot(abs(ref_FFTData(1:end/2,1))); ylim([0 4e6])
    figure(); imagesc( imadjust(mat2gray(abs(ref_FFTData(1:end/2,:))))); colormap(gray)
    figure(); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_FFTData(1:end/2,:)))))); colormap(gray);

    calSigMin = [29,40];
    winFunc = zeros(size(ref_FFTData));
    winFunc(calSigMin(1):calSigMin(2),:) = 1;

    cal_FFTData = ref_FFTData .* winFunc;
    cal_RawData = ifft(cal_FFTData);

    figure(); plot(real(cal_RawData()));
%     figure(); plot(real(cal_RawData(1:end/2,1)));
%     figure(); imagesc( imadjust(mat2gray(abs(cal_RawData(1:end/2,:))))); colormap(gray)
%     figure(); imagesc( imadjust(mat2gray(20 .* log10(...
%         abs(cal_RawData(1:end/2,:)))))); colormap(gray);

end