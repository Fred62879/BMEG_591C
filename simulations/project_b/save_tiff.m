function save_tiff(data, fname)

%     scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(data(1,:,:))))));
    scaled = imadjust(mat2gray(abs(squeeze(data(1,:,:)))));
    imwrite(scaled,fname);

    for i = 2:size(data,1)
%         scaled = imadjust(mat2gray(20 .* log10(abs(squeeze(data(i,:,:))))));
        scaled = imadjust(mat2gray(abs(squeeze(data(i,:,:)))));
        imwrite(scaled,fname,"WriteMode","append")
    end
end