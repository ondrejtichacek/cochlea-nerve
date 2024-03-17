function [ varargout ] = loadv( FILE, varargin )
%LOADV 

varargout{numel(varargin)} = [];

switch class(FILE)
    case 'char'
        for i = 1:numel(varargin)
            tmp = load(FILE, varargin{i});
            varargout{i} = tmp.(varargin{i});
        end
    case 'matlab.io.MatFile'
        for i = 1:numel(varargin)
            varargout{i} = FILE.(varargin{i});
        end
    otherwise
        error('Unsupported FILE class %s', class(FILE))
end

end

