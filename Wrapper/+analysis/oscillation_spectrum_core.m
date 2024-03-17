function [ ppp, P, PP, args ] = oscillation_spectrum_core( PS, freq, FREQUENCY, args )
%OSCILLATION_SPECTRUM_CORE
arguments
    PS double
    freq double
    FREQUENCY (1,1) double
    args.delta (1,1) double = 0.02
    args.tol (1,1) double = 0.05
    args.maxf (1,1) double = 30000
end

maxf_ind = find(freq >= args.maxf, 1);
if isempty(maxf_ind)
    maxf_ind = numel(freq);
end

if isnan(FREQUENCY)
    f_plus = +inf;
    f_minus = 0;
else
    FAC = max([0.1,(1-log10(FREQUENCY/1000))]);

    f_plus = exp(log(FREQUENCY)*(1+args.delta*FAC));
    f_minus = exp(log(FREQUENCY)*(1-args.delta*FAC));
end

f_plus_ind = min([maxf_ind, find(freq >= f_plus, 1)]);
f_minus_ind = max([1, find(freq >= f_minus, 1)]);

assert(f_plus_ind > f_minus_ind)

% f = freq(1:maxf_ind);

P = PS(1:maxf_ind,:);

PP = P(f_minus_ind:f_plus_ind,:);
P = [P(1:f_minus_ind,:); P(f_plus_ind:end,:)];

P = P / max(PP(:));

Ptol = P(P >= args.tol);

ppp = 1e6*sum(Ptol(:))/numel(P);

args.f_plus = f_plus;
args.f_minus = f_minus;

end

