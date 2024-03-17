function mustBeScalar(A)
%MUSTBESCALAR Validate that value is scalar or issue error

    if numel(A) > 1 % assume empty is scalar
        error('Value must be scalar.');
    end
end


