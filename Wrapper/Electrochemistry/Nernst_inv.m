function C = Nernst_inv(V, args)
%NERNST
arguments
    V
    args.charge = 1;
    args.temperature = 37 + 273.15; % [K]
    args.conc_in = [];
    args.conc_out = [];
end

z = args.charge;
T = args.temperature;

F = 96485.3328959;   % [C/mol]
R = 8.314459848;     % [J/mol/K]

if and(isempty(args.conc_in), isempty(args.conc_out)) || ...
   ~xor(isempty(args.conc_in), isempty(args.conc_out)) 
    error('exactly one of conc_in or conc_out has to be specified')
end

% V = R.*T./z./F.*log(Cout./Cin); % [J/mol/K*K/C*mol] = [J/C] = [V]
% exp(V.*z.*F./R./T) = Cout./Cin

if isempty(args.conc_in)
    Cout = args.conc_out;
    Cin = Cout./exp(V.*z.*F./R./T);
    C = Cin;
end

if isempty(args.conc_out)
    Cin = args.conc_in;
    Cout = Cin.*exp(V.*z.*F./R./T);
    C = Cout;
end
    
end

