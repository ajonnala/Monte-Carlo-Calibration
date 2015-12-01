function [final] = mserror(base_exp, base_st, comp_exp, comp_st)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
base_i = 1;
comp_i = 1;

sum = 0;

while (comp_st(comp_i) <= base_st(base_i))
    sum = sum + (base_exp(base_i) - comp_exp(comp_i))^2;
    comp_i= comp_i + 1;
end
tol=1e-10;
for base_i = 2:length(base_st)
    while (comp_st(comp_i) <= base_st(base_i))
        d1 = (base_st(base_i) - comp_st(comp_i));
        d2 = comp_st(comp_i) - base_st(base_i-1);
        if (abs(d1 -d2) < tol) | (d1 < d2)
            sum = sum + (base_exp(base_i) - comp_exp(comp_i))^2;
        else
            sum = sum + (base_exp(base_i-1) - comp_exp(comp_i))^2;
        end
        if (comp_i == length(comp_st)) break;
        else
        comp_i= comp_i + 1;
        end
    end
end
final = sqrt(sum/length(base_st));
end

