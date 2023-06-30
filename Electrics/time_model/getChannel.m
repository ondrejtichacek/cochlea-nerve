function [ Channel, sel, fun ] = getChannel( y, channel, mnaopt )
%GETChannel
arguments
    y (:,:) double
    channel (1,1) struct
    mnaopt (1,1) mnaOpt
end

Numstacks = mnaopt.Numstacks;
spec_size = mnaopt.circuit.num_v + mnaopt.circuit.num_node;

orders = [mnaopt.channels.order];

m = Numstacks * sum(orders(1:(channel.index-1)));

% located at the end of y with stride 2 (first eq. => -1)
jj = m + channel.order * (1:Numstacks) - 1;

sel = spec_size*Numstacks + jj;

if isempty(y)
    Channel = [];
else
    Channel(:,1:Numstacks) = y(:,sel);
end

fun = @(y) y(sel);

end