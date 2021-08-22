classdef OpticalFlowImageRegistrator < BilinearImageRegistrator
    properties(Access = private)
        factor;
        downsized_target;
        downsized_shifted_x_target;
        downsized_shifted_y_target;
        ix;
        iy;
    end
    methods(Access = public)
        function obj = OpticalFlowImageRegistrator(target,central_part_to_use, factor) 
            if(nargin<3)
                factor = 4;
            end
            if(nargin<2)
                central_part_to_use = [];
            end
            obj = obj@BilinearImageRegistrator(target,central_part_to_use);
            
            obj.factor = factor;
            
            obj.downsized_target = obj.resize(obj.get_central(obj.target),obj.factor);
            obj.downsized_shifted_x_target = obj.resize(obj.get_central(obj.shift(obj.target,[1 0])),obj.factor);
            obj.downsized_shifted_y_target = obj.resize(obj.get_central(obj.shift(obj.target,[0 1])),obj.factor);
            obj.ix = obj.downsized_target - obj.downsized_shifted_x_target;
            obj.iy = obj.downsized_target - obj.downsized_shifted_y_target;
        end
        
        function [d] = register_(obj,source,d0)
            if(nargin<3 && any(d0))
                shifted_source = obj.shift(source,d0);
            else
                shifted_source = source;
            end
            
            downsized_source = obj.resize(obj.get_central(shifted_source),obj.factor);
            it = downsized_source - obj.downsized_target;
            
            d = [obj.ix(:) obj.iy(:)]\it(:);
        end
        
        function [d, shifted] = register(obj,source,d0)
            assert(ismatrix(source));
            assert(all(size(source) > size(obj.central_part_of_target_with_moment.image)));
            if(nargin<3)
                d0 = [0 0];
            end
            d = register_(obj,source,d0);
            if(nargout>1)
                shifted = obj.shift(source,d);
            end
        end
    end
    methods(Static)
        function resized = resize(source,factor)
            if(numel(factor)==1)
                factor = [factor, factor];
            end
            s =size(source);
            S = ceil(s./factor);
            if(~isequal(s,S.*factor))
                source = padarray(source,S.*factor-s,'replicate','post');
            end
                resized = permute(mean(mean(reshape(source, [factor(1) S(1) factor(2) S(2)]),1),3),[2 4 1 3]);
        end
    end
end