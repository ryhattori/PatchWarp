classdef ImageBasis < handle
    properties(SetAccess = immutable, GetAccess = public)
        basis
        basis2d
        sums
        dim
        length
        metric
    end
    properties(Access = public)
        dot_product_to_target
        target_
    end
    properties(Access = public,Dependent)
        target
    end
    methods
        function set.target(obj,value)
            assert(isa(value,'ImageWithMoment'));
            assert(value.length == obj.length);
            obj.target_ = value;
            obj.dot_product_to_target = double(obj.target.image(:))'*double(obj.basis2d) - obj.target.sum * obj.sums/obj.length;
        end
        function value = get.target(obj)
            value = obj.target_;
        end
    end
    methods(Access = public)
        function obj = ImageBasis(basis)
            assert(ndims(basis)<=3);
            obj.basis = basis;
            obj.basis2d = reshape(basis,[size(basis,1)*size(basis,2) size(basis,3)]);
            obj.dim = size(obj.basis2d,2);
            obj.length = size(obj.basis2d,1);
            obj.sums = sum(obj.basis2d);
            obj.metric = double(obj.basis2d') * double(obj.basis2d) - double(obj.sums')*double(obj.sums/(obj.length));
            obj.target_ = [];
        end
        function cor = corr(obj,c1,c2)
            c1 = double(c1(:));
            assert(numel(c1) == obj.dim);
            assert(any(c1));
            c2 = double(c2(:));
            assert(numel(c2) == obj.dim);
            assert(any(c2));
            cor = (c1'*obj.metric*c2)/sqrt((c1'*obj.metric*c1)*(c2'*obj.metric*c2));
            cor = max(min(cor,1),-1);
        end

        function [tc,g] = corr_to_target(obj,c)
            assert(~isempty(obj.target_));
            c = double(c(:));
            assert(numel(c) == obj.dim);
            assert(any(c));
            
            dp = obj.dot_product_to_target*c;
            sqc = sqrt((c'*obj.metric*c));
            sqt = sqrt(obj.target.rss);
            tc = dp/(sqc*sqt);
            if(nargout>1)
                 g = ( obj.dot_product_to_target' /sqc - dp*obj.metric*c/sqc^3)/sqt;
            end
                
        end
    end
    
end