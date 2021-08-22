classdef ImageRegistratorPyramid < handle
    methods(Access = public)
        function obj = ImageRegistratorPyramid(target,central_part_to_use,pyramid_depth,pyramid_proportion,ir_constructor)
            if(nargin < 3) 
                pyramid_depth = 3;
            end
            if(nargin < 3)
                pyramid_proportion = 0.25;
            end
            if(nargin < 5)
                ir_constructor = @ImageRegistrator;
            end
            assert(pyramid_depth>0);
            
            obj.ir_constructor = ir_constructor;
            obj.pyramid_depth = pyramid_depth;
            obj.pyramid_proportion = pyramid_proportion;
            
            obj.ir_pyramid = cell(obj.pyramid_depth,1);
            current_target = target;
            obj.ir_pyramid{1} = ir_constructor(current_target,central_part_to_use);
            for i_p = 2:obj.pyramid_depth
                current_target = imresize(current_target,pyramid_proportion,'box');
                obj.ir_pyramid{i_p} = ir_constructor(current_target,central_part_to_use);
            end
        end
    end
    methods(Access = public)
        function [d] = register(obj,source,d0)
            if(nargin<3)
                d0 = [0 0];
            end
            
            source_pyramid = cell(obj.pyramid_depth,1);
            current_source = source;
            source_pyramid{1} = current_source;
            for i_p = 2:obj.pyramid_depth
                current_source = imresize(current_source,obj.pyramid_proportion,'box');
                source_pyramid{i_p} = current_source;
            end
            d0_scaled = floor(d0 *(obj.pyramid_proportion^obj.pyramid_depth));
            for i_p = obj.pyramid_depth:-1:1
                [d] = obj.ir_pyramid{i_p}.register_(...
                     source_pyramid{i_p},d0_scaled);
                 d0_scaled = d/obj.pyramid_proportion;
            end
                        
            
        end
    end
    
    properties(SetAccess = immutable, GetAccess = protected)
        ir_constructor
        ir_pyramid
        pyramid_proportion
        pyramid_depth
    end
end