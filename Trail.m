classdef Trail
    properties (Access = public)
        pos;
        alg;
    end
    methods
        function obj = Trail(posT, algT)
            obj.pos = posT;
            obj.alg = zeros(1,3);
            obj.alg(algT) = 1;
        end
    end
end