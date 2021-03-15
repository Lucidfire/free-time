%Function for simulating continuous random variates via an
%acceptance/rejection sampling procedure.
%   n is number of iid variates to generate
%
%   lims are the upper and lower bounds to which variates are restricted
%   expressed as [lower upper] TODO: add suport for infinite bounds and
%   make this arg optional
%
%   pdf is a function handle containing the pdf to simulate variates from
%
%Copyright (c) 2021, Anthony Zinos
% All rights reserved.
% 
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. 
function result = continuousRejection(n,lims,pdf)
    p = pdf(lims(1):1/1000:lims(2));
    c = max(p);
    result = zeros(1,n);
    i=1;
    while(i<=n)
        Y = lims(1)+rand*(lims(2)-lims(1));
        u = rand;
        if u<(pdf(Y)/(c))
            result(i) = Y;
            i = i+1;
        end
    end 
end