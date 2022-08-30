function [outImage] = subSampling(image, subSamplingFactor)
[rows cols] = size(image);
outImage = image(1:subSamplingFactor:rows,1:subSamplingFactor:cols);
end