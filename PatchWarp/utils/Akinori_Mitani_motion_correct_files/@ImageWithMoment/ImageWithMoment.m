classdef ImageWithMoment < handle

    properties(SetAccess = immutable, GetAccess = public)
        image
        length
    end
    properties(SetAccess = private,GetAccess = public)
        sum
        rss
    end
    methods
        function value = get.sum(obj)
            if(isempty(obj.sum))
                value = sum(sum(obj.image));
                obj.sum = value;
            else
                value = obj.sum;
            end
        end
        function value = get.rss(obj)
            if(isempty(obj.rss))
                value = obj.image(:)'*obj.image(:) - obj.sum ^2 / obj.length;
                obj.rss = value;
            else
                value = obj.rss;
            end
        end
    end
    methods(Access = public)
        function obj = ImageWithMoment(image)
            assert(ismatrix(image));
            obj.image = double(image);
            obj.length = numel(image);
            obj.sum = [];
            obj.rss = [];
        end
        function cor = corr(obj1,obj2)
            cor = covariance(obj1,obj2) ./ sqrt(obj1.rss * obj2.rss);
        end
        function cov = covariance(obj1,obj2)
            if(~isa(obj2,'ImageWithMoment'))
                assert(ismatrix(obj2));
                obj2=ImageWithMoment(obj2);
            end
            assert(isa(obj2,'ImageWithMoment'));
            
            cov = obj1.image(:)' * obj2.image(:) - obj1.sum * obj2.sum / obj1.length;
        end
    end
    
end