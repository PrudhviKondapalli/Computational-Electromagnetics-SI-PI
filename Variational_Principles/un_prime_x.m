function [out] = un_prime_x(n)

syms x;  % Declare x as a symbolic variable
out = diff((x-1)*x^n);  % Calculate x^n

end