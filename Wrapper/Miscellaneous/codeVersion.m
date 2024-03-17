function [ version_string, VERSION, REVISION ] = codeVersion( )
%CODEVERSION 

VERSION = 1;
REVISION = 0;

version_string = sprintf('%d.%03d',VERSION, REVISION);


end

