function [out] = new_element(col,A,B)
%New element calculation 
%Here A is the last calculated latest element
%Here B is the newly calculated latest element
out = (((2^(col-1))*A) - B)/((2^(col-1))-1);
end