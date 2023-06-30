function [ Current, sel ] = getCurrent( y, NAME, mnaopt )
%GETCurrent

numNode = mnaopt.circuit.num_node;
numV = mnaopt.circuit.num_v;
Numstacks = mnaopt.Numstacks;

%     for jj = 1:numV
%         res.Current(ts:te, 1:Numstacks, jj) = Y(:,((1:Numstacks)-1)*(numNode+numV)+numNode+jj)

if isempty(NAME)
    Current = zeros(size(y,1), Numstacks, numV);
    sel{numV} = [];
    for jj = 1:numV
        sel{jj} = ((1:Numstacks)-1)*(numNode+numV) + numNode + jj;
        Current(:, 1:Numstacks, jj) = y(:,sel{jj});
    end
else
    
    switch NAME
        case 'IHC'
            jj = find(strcmp('VI', {mnaopt.circuit.vsources.name}),1);
        case 'OHC'
            jj = find(strcmp('VO', {mnaopt.circuit.vsources.name}),1);
        case 'IHC_MET'
            jj = find(strcmp('VIapical', {mnaopt.circuit.vsources.name}),1);
        case 'OHC_MET'
            jj = find(strcmp('VOapical', {mnaopt.circuit.vsources.name}),1);
        case 'StV'
            jj = find(strcmp('VStV', {mnaopt.circuit.vsources.name}),1);
        otherwise
            error('Unknown option %s', NAME);
    end    
    
    sel = ((1:Numstacks)-1)*(numNode+numV)+numNode+jj;
    if isempty(y)
        Current = [];
    else
        Current(:,1:Numstacks) = y(:,sel);
    end
    
end

end