%Perform Monte Carlo integration on specified function func (passed as
%function handle) over domain given in matrix form:
%domain = [a1 b1;
%          a2 b2;
%          a3 b3]
%where ai,bi are the boundaries of integration in the ith dimension.
%
%samples specifies the number of random samples to use in calculating the
%integral. Higher samples leads to greater precision at the cost of
%computation time. Default value is 10^6.
%
%Copyright (c) 2021, Anthony Zinos
% All rights reserved.
% 
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. 

function integral = mcInt(func,domain,samples)

%Default values for optional arg 3
if(nargin<3)
    samples = 10^6;
end
dim = size(domain,1);

%Check formatting of domain
if(~isequal(size(domain),[dim,2]))
    disp('Error: domain must be dimensionsx2.');
    return;
end

%Check positivity and integrality of samples
if(~(samples==floor(samples))||samples<0)
    disp('Error: number of samples must be a positive whole number.');
    return;
end
%Check ordering and non equality of bounds
if(sum(domain(:,1)>domain(:,2))>0)
    disp('Error: upper bounds of integration are required to be strictly greater than lower bounds.');
    return;
end

Jacobians = cell(1,dim);
transforms = cell(1,dim);
u = rand(dim,samples);
for i = 1:dim
    if isinf(domain(i,1))
        a=0;
    else
        a = domain(i,1);
    end
    if isinf(domain(i,2))
        b=0;
    else
        b = domain(i,2);
    end
    if(isequal(domain(i,:),[a Inf]))
        Jacobians{i} = 1./((1-u(i,:)).^2);
        transforms{i} = makeAiCase(a);
    end
    if(isequal(domain(i,:),[-Inf Inf]))
        Jacobians{i} = 1./(u.*(1-u));
        transforms{i} = @IiCase;
    end
    if(isequal(domain(i,:),[-Inf b]))
        Jacobians{i} = 1./(u.^2);
        transforms{i} = makeIbCase(b);
    end
    if(isequal(domain(i,:),[a b]))
        Jacobians{i} = b-a;
        transforms{i} = makeAbCase(a,b);
    end
end
args = cell(1,dim);
for i = 1:dim
    t = transforms{i};
    u(i,:) = t(u(i,:));
    args{i} = u(i,:);
end
if nargin(func)>1
    y = func(args{:});
else
    y = func(u);
end
for i = 1:dim
   y = y.*Jacobians{i}; 
end
integral = mean(y,'all');
    function transform = makeAbCase(a,b)
        transform = @abCase;
        function result = abCase(x)
            result = a+(b-a).*x;
        end
    end
    function transform = makeAiCase(a)
        transform = @aiCase;
        function result = aiCase(x)
            result = a+(x./(1-x));
        end
    end
    function transform = makeIbCase(a)
        transform = @ibCase;
        function result = ibCase(x)
            result = b+1-(1./x);
        end
    end
    function result = IiCase(x)
        result = log(x./(1-x));
    end
end