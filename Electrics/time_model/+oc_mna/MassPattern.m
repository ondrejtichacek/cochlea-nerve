function mass_pattern = MassPattern(Ic, channels, Numstacks, spec_size)
% MASSPATTERN Calculates sparsity pattern of Mass matrix
%
% d/dy M(t,y) ?!= 0
%
if all(structfun(@isempty, Ic))
    
    % mass_pattern = [];
    % return

    n = Numstacks*spec_size;

    mass_pattern = spalloc(n,n,10*Numstacks);
    
else
    n = Numstacks*spec_size;

    mass_pattern = spalloc(n,n,10*Numstacks);

    deps = fieldnames(Ic);
    for i = 1:numel(deps)
        dep = deps{i};
%         if any(strcmp(dep, {'vihc_ss', 'vohc_ss'}))
%             continue
%         end
        if ~isempty(Ic.(dep))            
            for j = 1:numel(Ic.(dep))
                II = Ic.(dep){j} ~= 0;

                mp = kron(spdiags(ones(Numstacks,1),0,Numstacks,Numstacks), II);

                mass_pattern = mass_pattern | mp;
            end
        end
    end

end

if ~isempty(channels)
    
    V = struct( ...
    'IHC', -ones(Numstacks,1), ...
    'OHC', -ones(Numstacks,1), ...
    'IHC_ss', -ones(Numstacks,1), ...
    'OHC_ss', -ones(Numstacks,1));
    
    [~, ~, ~, ~, ~, ~, MM_V] = ChannelspOpen(V, channels, Numstacks);

    [a,b] = size(mass_pattern);
    [aa,bb] = size(MM_V);

    MM_V = MM_V ~= 0;

    mass_pattern =  [ mass_pattern, zeros(a,bb);
                       zeros(aa,b),     MM_V ];
end


mass_pattern = sparse(mass_pattern);


end