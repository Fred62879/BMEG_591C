

function simulation_feb_7

    function rawData_rescaled = reSampling_CalSig(fftData, cal_RawData):
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

    fn = 'RawOCT.mat';
    load(fn)
    %raw_data = load(fn);
    %raw_data = getfield(raw_data,'rawData'); %_array');

    ref_rawData = rawOCT(:,:,250);
    ref_FFTData = fft(hilbert(ref_rawData));

    calSigMin = [29,40];
    winFunc = zeros(size(ref_FFTData));
    windFunc(calSigMin(1):calSigMin(2),:) = 1;

    cal_FFTData = ref_FFTData .* winFunc;
    cal_RawData = ifft(cal_FFTData);

    figure(); plot(abs(ref_FFTData(1:end/2,:)));
    figure(); imagesc( imadjust(mat2gray(abs(ref_FFTData(1:end/2,:))))); colormap(gray)
    figure(); imagesc( imadjust(mat2gray(20 .* log10(...
        abs(ref_FFTData(1:end/2,:)))))); colormap(gray);

    figure(); plot(real(cal_RawData));

end