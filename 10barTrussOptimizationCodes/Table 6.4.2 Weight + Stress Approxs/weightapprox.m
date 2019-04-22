function [gw] = weightapprox(x,x0,w,dwdA)
% This computes the weight approximation for each step of optimization;
gw = double(w + sum((x - x0)*dwdA'));
end