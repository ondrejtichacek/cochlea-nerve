function [ ] = disp_title( STRING, BORDER, WIDTH )
%DISP_TITLE Prints formated title on screen

if nargin < 3
    WIDTH = 80;
end

if nargin < 2
    BORDER = '=';
end

n = WIDTH - length(STRING) - 2; % 2 for spaces

m = (n - mod(n,2)) / 2;

BORDER_STRING = repmat(BORDER,1,m+1);

fprintf('\n%s\n', [ BORDER_STRING(1:m), ' ', STRING, ' ', BORDER_STRING(1:m+mod(n,2)) ])

end

