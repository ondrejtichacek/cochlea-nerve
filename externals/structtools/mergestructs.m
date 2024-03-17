function S = mergestructs(S, T, args)
%MERGESTRUCTS Combines fields of both input structures into a single structure
%   Output structure will contain all fields and their corresponding values of
%   the two input structures. If a field is present in both input structures,
%   the field in S is used, and a warning is produced. This warning can be
%   turned off using warning('off', "structtools:commonfield")
%
%   Example:
%   S.one = 1;
%   T.two = 2;
%   mergestructs(S,T)
%
%   ans = 
%   struct with fields:
%     one: 1
%     two: 2
%
%   Inputs:
%   S       first structure to merge
%   T       second structure to merge
%
%   Output:
%   S       combination structure

    arguments
        S(1,1) struct
        T(1,1) struct
        args.nowarn (1,1) logical = false
    end

    f = fieldnames(T);
    for i = 1:length(f)
        if isfield(S, f{i})
            if args.nowarn == false
                warning("structtools:commonfield",...
                    "Field '%s' exists in both structures.", f{i});
            end
            continue
        end
        S.(f{i}) = T.(f{i});
    end
end