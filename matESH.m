function [mat] = matESH(vec)
%matES matrixizes a vector from genEigH
mat = flipud(reshape(vec, [round(sqrt(length(vec))) round(sqrt(length(vec)))])');
end

