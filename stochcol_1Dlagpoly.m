function [lagpoly, weight]= stochcol_1Dlagpoly(node, nodes, varargin)
%STOCHCOL_1DLAGPOLY function handle and weight of 1D Lagrange polynomial
%
% [lagpoly, weight]= stochcol_1Dlagpoly(node,nodes,varargin)
%
%   Latest modification: AS; 25 November 2022
% Copyright (c) 2019 F. Xu

if length(nodes)==1
    if numel(varargin) == 0
        lagpoly = @(y) 1;
        weight = 2;
    else 
        lagpoly = @(y) 1;
        weight = 1; %\int_{-1}^{1} pdf dy = 1 by definition of pdf
    end
else
    nodes = setdiff(nodes,node);
    n = length(nodes);
    lagpoly = @(y) 1;
    for k = 1:n
        lagpoly = @(y) lagpoly(y).*(y-nodes(k))./(node-nodes(k));
    end
    if numel(varargin) == 0
        weight = integral(lagpoly,-1,1);
    else
        pdf = varargin{1};
        lagpoly = @(y) lagpoly(y).*pdf(y);
        weight = integral(lagpoly,-1,1);
    end
end
end
