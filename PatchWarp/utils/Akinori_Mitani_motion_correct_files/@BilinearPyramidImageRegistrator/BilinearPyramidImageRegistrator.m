classdef BilinearPyramidImageRegistrator < PyramidImageRegistrator & BilinearImageRegistrator

    methods(Access = public)
        function obj = BilinearPyramidImageRegistrator(varargin)
            obj = obj@PyramidImageRegistrator(varargin{1:min(3,length(varargin))});
            obj = obj@BilinearImageRegistrator(varargin{1:min(2,length(varargin))});
        end
        function [d, c] = register_(obj,d0)
            [d_int] = register_@PyramidImageRegistrator(obj,round(d0));
            
            [d, c] = register_@BilinearImageRegistrator(obj,d_int,true);
        end
    end
end