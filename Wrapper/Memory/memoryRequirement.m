function [ variableSize, VarProp, Prop ] = memoryRequirement( FILE, variableName )
%MEMORYREQUIREMENT 

if ischar( FILE ) % FILE is file name
    Prop = whos('-file',FILE);
else % FILE is file handle
    Prop = whos(FILE);
end


VarProp = [];
variableSize = 0;

if nargin < 2

    variableSize = sum([Prop.bytes]);
    
else
    
    for i = 1:length(Prop)
        if strcmp(Prop(i).name, variableName)
            VarProp = Prop(i);
            continue;
        end
    end

    if ~isempty(VarProp)
        variableSize = VarProp.bytes;
    else
        warning('specified variable not found in given file');
    end
end

end

