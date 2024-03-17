function [ varargout ] = copy_all( varargin )
%COPY_ALL 

varargout{numel(varargin)} = [];

for i = 1:numel(varargin)
    
    varargout{i} = copy(varargin{i});
    
end


end

