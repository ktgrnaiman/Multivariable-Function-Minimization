classdef Minima
    properties (Access = public)
        pos; f; trail;
    end
    methods
        function obj = Minima(posT, fT, posI, alg)
            obj.pos = posT;
            obj.f = fT;
            obj.trail = Trail(posI, alg);
        end
    end
end