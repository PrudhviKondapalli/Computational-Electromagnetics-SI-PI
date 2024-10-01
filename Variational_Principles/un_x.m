function [out] = un_x(n)

syms x;  % Declare x as a symbolic variable
out = (x-1)*x^n;  % Calculate x^n

end