function [vec] = vecES(mat)
%vecES vectorizes a matrix for use in ESdtNL code
vec = reshape(flipud(mat)',[numel(mat) 1]);
end

