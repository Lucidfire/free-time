%Function for simulating discrete random variates via an
%acceptance/rejection sampling procedure.
%   indices are the x values over which the pmf is defined, passed in as an
%   array. TODO: there should be a seperate n and indices passed in,
%   instead of using length(indices) as n.
%
%   pmf is a function handle containing the pmf to simulate variates from
%
%Copyright (c) 2021, Anthony Zinos
% All rights reserved.
% 
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. 
function result = discreteRejection(indices,pmf)
    n = length(indices);
    p = pmf(indices);
    q = (1/n)*ones(1,n);
    c = max(p./q);
    result = zeros(1,n);
    i=1;
    while(i<=n)
        u = rand(1);
        Y = indices(ceil(u*n));
        u = rand(1);
        if u<(p(Y)/(c*q(Y)))
            result(i) = Y;
            i = i+1;
        end
    end    
end