function [outImage] = upSampling(image, upSamplingFactor)
[rows cols] = size(image);
newRows = rows*upSamplingFactor;
newCols = cols*upSamplingFactor;
rowStart = 1;
for rowsIndex=1:upSamplingFactor:newRows
    colStart = 1;
    for columnIndex=1:upSamplingFactor:newCols
        outImage(rowsIndex:rowsIndex+upSamplingFactor-1,columnIndex:columnIndex+upSamplingFactor-1) = image(rowStart,colStart);
    colStart = colStart + 1;
    end
    rowStart = rowStart + 1;
end
