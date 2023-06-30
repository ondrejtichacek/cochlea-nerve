function [ varargout ] = BMdispl( t, t_mech, mechfs, varargin )
%BMDISPL

tind = 1 + round((t-t_mech(1))*mechfs);

if tind <= 0 || tind > numel(t_mech)
    error('error in computed time index tind = %d, should be from 0 to %d', tind, numel(t_mech))
end

if abs(t_mech(tind) - t) > 1.5/mechfs % factor 1.5 to account for numerics
    warning('Time distance from the mechanical time sample is too large: %g\n\tshould be max. %g', abs(t_mech(tind) - t), 1/mechfs)
end

varargout = cell(numel(varargin), 1);

for i = 1:numel(varargin)
    loadedBM = varargin{i}; % it does not need to be BM, maybe other mech variable such as TM
    if tind > numel(t_mech)
        BM = loadedBM(end,:);
        warning('BMx size exceeded');
    else
        BM = loadedBM(tind,:);
    end
    varargout{i} = BM;
end

end
