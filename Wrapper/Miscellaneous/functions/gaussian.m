function [ y ] = gaussian(x, A, mu, sigma)

y = A .* exp(-(((x-mu)./sigma).^2)/2);

end

