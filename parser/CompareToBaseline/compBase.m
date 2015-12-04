function [sumError] = compBase(tensor_base, tensor_comp)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

a = tensor_base(1,1);
b = tensor_base(2,1);
base_exp = a{1,1};
base_st = b{1,1};

a = tensor_comp(1,1);
b = tensor_comp(2,1);
comp_exp = a{1,1};
comp_st = b{1,1};


error = ones(1,260);
sumError = 0;

for ts = 1:size(base_exp,1)
   error(ts) = splineError(base_exp(ts, :), base_st(ts, :), comp_exp(ts, :), comp_st(ts, :));
   sumError = sumError + splineError(base_exp(ts, :), base_st(ts, :), comp_exp(ts, :), comp_st(ts, :));
end

