function [out] = u0_prime_x

syms x;  % Declare x as a symbolic variable
out = diff(x);  % Calculate x^n

end