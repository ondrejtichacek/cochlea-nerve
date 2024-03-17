function [ str ] = sprintfBytes( formatSpec, varargin )
%SPRINTFBYTES 

[ val, unit ] = bestBytes( varargin{:} );

str = sprintf(formatSpec, val, unit);

end

