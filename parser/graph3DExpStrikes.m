function [] = graph3DExpStrikes(tensor, col, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


a = tensor(1,col);
b = tensor(2,col);
expe = a{1,1};
strikes = b{1,1};

surf(expe,strikes,'EdgeColor','none','FaceLighting','phong')
title(t)

