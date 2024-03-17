function [ PARAM, n, par_info ] = ParameterProduct( names, values )
%PARAMETERPRODUCT 
% Example:
% [Conf,n] = ParameterProduct( 'Amplitude', 0:10:100, 'Frequency', 1e3:1e3:1e4);
% Conf(sub2ind(n,1,3))
% 
% ans = 
% 
%   struct with fields:
% 
%     Amplitude: 0
%     Frequency: 3000

arguments (Repeating)
    names (1,:) char
    values {isvector(values)}
end

for i = 1:numel(values)
    if ischar(values{i})
        warning('Assuming %s is a single value, not a range.', values{i});
        values{i} = string(values{i});
    end
end

par_info = struct();
par_info.order = names;
for i = 1:numel(values)
    par_info.range.(names{i}) = values{i};
end

N = numel(names);

for i = 1:N
    assert(numel(values{i}) > 0, 'Range %s seems to be empty', names{i});
end

n = cellfun(@numel, values, 'UniformOutput', true);
sub = cell(1,numel(n));

S = cell(2, N);

[S{1,:}] = names{:};

for i = 1:N
    S{2,i} = cell(1,prod(n));
end

par_info.sub = zeros(prod(n),N);

for j = 1:prod(n)
    [sub{:}] = ind2sub(n,j);
    par_info.sub(j,:) = [sub{:}];
    for i = 1:N
        S{2,i}{j}= values{i}(sub{i});
    end
end

PARAM = struct(S{:});


end