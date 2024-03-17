function S = flattenstruct(S, delim)
%FLATTENSTRUCT Converts a structure with substructures to a single-level struct
%   All fields in S of type struct will be removed and their fields added to the
%   parent struct, S. The new field names will be "fieldname_subfieldname". This
%   process is recursive.
%
%   Example:
%   S.field1 = "parent";
%   S.substruct.subfield1 = "child1"
%   S.substruct.subfield2 = "child2"
%   S.substruct.subsubstruct.subsubfield1 = "grandchild";
%   flattenstruct(S)
%
%   ans = 
%   struct with fields:
%                                field1: "parent"
%                   substruct_subfield1: "child1"
%                   substruct_subfield2: "child2"
%   substruct_subsubstruct_subsubfield1: "grandchild"
%
%   Inputs:
%   S       structure to flatten
%   delim   delimiter to use between subtructure and subfield name (Default '_')
%
%   Output:
%   S       flattened structure

    arguments
        S(1,1) struct
        delim(1,:) char = '_'
    end
    
    f = fieldnames(S);
    for i = 1:length(f)
        if ~isstruct(S.(f{i}))
            continue
        end
        B = flattenstruct(S.(f{i}));
        fB = fieldnames(B);
        for j = 1:length(fB)
            S.([(f{i}) delim (fB{j})]) = B.(fB{j});
        end
        S = rmfield(S, f{i});
    end
end
