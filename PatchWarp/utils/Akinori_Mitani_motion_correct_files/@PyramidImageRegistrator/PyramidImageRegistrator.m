classdef PyramidImageRegistrator < ImageRegistrator
    properties(SetAccess = immutable, GetAccess = protected)
        imageRegistratorPyramid
        pyramid_depth
        pyramid_multiply
    end
    methods(Access = public)
        function obj = PyramidImageRegistrator(target,central_part_to_use,pyramid_depth)
            if(nargin < 3) 
                pyramid_depth = 3;
            end
            if(nargin < 2) 
                central_part_to_use = [];
            end
            
            obj = obj@ImageRegistrator(target,central_part_to_use);
            
            obj.pyramid_multiply = 2;
            obj.pyramid_depth = pyramid_depth;
            obj.imageRegistratorPyramid = cell(obj.pyramid_depth,1);
            current_target = obj.target;
            for i_p = 1:obj.pyramid_depth
                current_target = imresize(current_target,1/obj.pyramid_multiply,'box');
                obj.imageRegistratorPyramid{i_p} = ...
                    ImageRegistrator(current_target,central_part_to_use);
            end
        end
        function [d, c] = register_(obj,d0)
            source_pyramid = cell(obj.pyramid_depth,1);
            current_source = obj.correlationCalculator.source;
            for i_p = 1:obj.pyramid_depth
                current_source = imresize(current_source,0.5,'box');
                source_pyramid{i_p} = current_source;
            end
            
            d0_scaled = floor(d0/(obj.pyramid_multiply^obj.pyramid_depth));
            for i_p = obj.pyramid_depth:-1:1
                d0_scaled = obj.pyramid_multiply*obj.imageRegistratorPyramid{i_p}.register(...
                    source_pyramid{i_p},d0_scaled);
            end
            
            [d, c] = register_@ImageRegistrator(obj,d0_scaled);
        end
    end
    
end