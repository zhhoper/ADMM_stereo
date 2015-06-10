function r = myRandom(varargin)
%myRandom Draw a random number in the interval [a,b]
switch nargin
    case 2
        a = varargin{1};
        b = varargin{2};
        r = (b-a)*rand() + a;
    case 3
        a = varargin{1};
        b = varargin{2};
        n = varargin{3};
        r = (b-a)*rand(n) + a;
    otherwise
        r = rand();
end

end

