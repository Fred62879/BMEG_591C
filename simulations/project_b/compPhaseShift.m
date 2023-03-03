function rawData_PhaseComp = compPhaseShift(ref_fftData_1D, fftData_2D, depthIdx)
    ref_Ascan = ref_fftData_1D;

    cplxConjX = fftData_2D...
        .* repmat(conj(ref_Ascan), [1 size(fftData_2D,2)]);

    calSigDepth = depthIdx;

    phaseSlope = (angle(cplxConjX(calSigDepth,:)) /calSigDepth)...
        .* linspace(1,length(ref_Ascan), length(ref_Ascan))';

    for i = 1:size(fftData_2D, 2)
        fftData_PhaseComp(:,i) = fftData_2D(:,i)...
            .* exp(-1j * phaseSlope(:,i));
    end

    rawData_PhaseComp = ifft(fftData_PhaseComp);
end
