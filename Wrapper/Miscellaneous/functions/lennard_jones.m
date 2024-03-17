function [f] = lennard_jones(r, epsilon, r0, n)
%LENNARD_JONES 


f = epsilon * ((r0./r).^(2*n) - 2 * (r0./r).^n);

end

