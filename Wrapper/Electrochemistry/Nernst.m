function V = Nernst(Cin, Cout, args)
%NERNST
arguments
    Cin
    Cout
    args.charge = 1;
    args.temperature = 310; % [K] ~ 37C
end

z = args.charge;
T = args.temperature;

F = 96485.3328959;   % [C/mol]
R = 8.314459848;     % [J/mol/K]


V = R.*T./z./F.*log(Cout./Cin); % [J/mol/K*K/C*mol] = [J/C] = [V]

end