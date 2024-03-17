function [val] = isemptystruct(S)
%ISEMPTYSTRUCT 

val = numel(fieldnames(S)) == 0;

end

