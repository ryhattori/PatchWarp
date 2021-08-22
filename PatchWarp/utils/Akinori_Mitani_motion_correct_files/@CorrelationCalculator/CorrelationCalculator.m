classdef CorrelationCalculator < handle
    properties(SetAccess = immutable, GetAccess = public)
        source
    end
    properties(SetAccess = immutable, GetAccess = protected)
        target_with_moment
    end
    
    properties(SetAccess = immutable, GetAccess = private)
        max_size
        center_offset
        center_index
    end
    properties(Access = private)
        computed_correlation
    end
    methods(Access = public)
        function obj = CorrelationCalculator(target_with_moment,source)
            assert(isa(target_with_moment,'ImageWithMoment'));
            
            obj.target_with_moment = target_with_moment;
            obj.source = source;
            
            obj.max_size = [1 1]+size(source)-size(target_with_moment.image);
            
            assert(all(obj.max_size>0));
            
            obj.computed_correlation = NaN( obj.max_size );
            obj.center_offset = ceil(obj.max_size/2);
            for i = 1:2
                obj.center_index{i} = (1:size(obj.target_with_moment.image,i))+obj.center_offset(i)-1;
            end
        end
        
        function cor = get(obj,d)
            assert(numel(d)==2);
            index = d+obj.center_offset;
            if(any(index<=0)||any(index>obj.max_size))
                cor = -2;
                return;
            end
            if(~isnan(obj.computed_correlation(index(1),index(2))))
                cor = obj.computed_correlation(index(1),index(2));
            else
                shifted_source_with_moment = ImageWithMoment(obj.source(obj.center_index{1}-d(1),obj.center_index{2}-d(2)));
                cor = corr(obj.target_with_moment,shifted_source_with_moment);
                obj.computed_correlation(index(1),index(2)) = cor;
            end
        end
    end
end