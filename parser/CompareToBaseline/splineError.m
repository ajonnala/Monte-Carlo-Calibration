function [final] = splineError(base_exp, base_st, comp_exp, comp_st)
%   Spline interpolation with linear extrapolation

interp = interp1(comp_st,comp_exp,base_st(base_st>min(comp_st) & base_st<max(comp_st)),'spline','extrap'); %interpolate
extrap = interp1(comp_st,comp_exp,base_st(base_st<=min(comp_st) | base_st>=max(comp_st)),'linear','extrap'); %exterpolate
ind = (base_st<=min(comp_st) | base_st>=max(comp_st));

new_comp_exp = zeros(size(base_st,1), size(base_st,2));
int_j = 1;
ext_j = 1;
for j = 1:length(base_st)
    if (ind(j))
        new_comp_exp(j) = extrap(int_j);
        int_j = int_j + 1;
    else
        new_comp_exp(j) = interp(ext_j);
        ext_j = ext_j + 1;
    end
end

%plot(base_st,base_exp, 'o', comp_st, comp_exp, '*', base_st, new_comp_exp, '+')

sum = 0;

 for i= 1:length(base_st)
     sum = sum + (base_exp(i) - new_comp_exp(i))^2;
 end
 
final = sqrt(sum/length(base_st));

end

