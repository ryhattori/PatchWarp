classdef ImageRegistrator
    properties(SetAccess = immutable, GetAccess=public)
        target
        central_part_to_use
        steps
    end
    properties(SetAccess = immutable, GetAccess=protected)
        central_part_of_target_with_moment
    end
    properties(SetAccess = protected)
        correlationCalculator
    end
    methods(Access = public)
        function obj = ImageRegistrator(target,central_part_to_use) 
            assert(ismatrix(target));
            if(nargin<2||isempty(central_part_to_use))
                central_part_to_use = 0.75;
            end
            
            assert(2>=numel(central_part_to_use));
            
            obj.target = target;
            assert(all(central_part_to_use>0) && all(central_part_to_use<=1));
            obj.central_part_to_use = central_part_to_use;
            
            obj.central_part_of_target_with_moment = ImageWithMoment(...
                obj.get_central(target));
            
            obj.steps = {[0 0],[-1 0],[1 0],[0 -1],[0 1],[-1 -1],[-1 1],[1 -1],[1 1]};
        end
        function central_part = get_central(obj,image,d)
            if(nargin<3 || isempty(d))
                d = [0 0];
            end
            if(numel(d)<2)
                d = [d d];
            end
            edge_part = (1-obj.central_part_to_use)/2;
            central_indices_begin = floor(size(image).*edge_part)+1;
            central_indices_end = ceil(size(image).*(1-edge_part));
            central_part = image((central_indices_begin(1):central_indices_end(1))-d(1),...
                    (central_indices_begin(2):central_indices_end(2))-d(2));
        end
        
        function [d, shifted, c] = register(obj,source,d0)
            if(nargin<3)
                d0 = [0 0];
            end
            if(size(d0,1)<size(source,3))
				d0=repmat(d0(1,:),size(source,3),1);
			end
			if(size(source,3)==1 && size(d0,2)==1)
				d0 = d0';
			end
			
            assert(isequal(size(source,1),size(obj.target,1)));
            assert(isequal(size(source,2),size(obj.target,2)));
            
			I = size(source,3);
			d=zeros(I,2);
			c=zeros(I,1);
			
			for i=1:I
				obj.correlationCalculator = CorrelationCalculator(obj.central_part_of_target_with_moment, source(:,:,i));
				[d(i,:), c(i)] = register_(obj,d0(i,:));
			end
			if(nargout>1)
				shifted = obj.shift(source,d);
			end
        end
    end
    methods(Access = public)
        function [d c] = register_(obj,d0)
            correlations = zeros(numel(obj.steps),1);
            d = d0;
            max_iter = ceil(sum(size(obj.target)-size(obj.central_part_of_target_with_moment.image))/2);
            for j = 1:max_iter
                for i = 1:numel(obj.steps)
                    correlations(i) = obj.correlationCalculator.get(d+obj.steps{i});
                end
                [c,i]=max(correlations);
                if(i==1)
                    break
                else
                    d = d+obj.steps{i};
                end
            end
        end
    end
    methods(Static)
        function shifted = shift(source,d)
            shifted = NaN(size(source),'single');
            shifted((1+max(d(1),0)):(end+min(0,d(1))), (1+max(d(2),0)):(end+min(0,d(2))))...
                = source((1-min(0,d(1))):(end-max(d(1),0)), (1-min(0,d(2))):(end-max(d(2),0)));
        end
    end
    
end