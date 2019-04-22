function [gw] = weight(x,L,rhow)
gw = rhow*x*L';
end