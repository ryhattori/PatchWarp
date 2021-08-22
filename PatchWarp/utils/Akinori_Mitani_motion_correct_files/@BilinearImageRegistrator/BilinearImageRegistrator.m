classdef BilinearImageRegistrator < ImageRegistrator
    methods(Access = public)
        function obj=BilinearImageRegistrator(varargin)
            obj=obj@ImageRegistrator(varargin{:}); 
        end
        function [d, c] = register_(obj,d0,registered)
            if(nargin<3)
                registered=false;
            end
            if(registered)
                d0_int = d0;
            else
                d0_int = register_@ImageRegistrator(obj,round(d0));
            end
            
            cc = zeros(3);
%             for i_dx=1:3
%                 dx=i_dx-2;
%                 for i_dy=1:3
%                     dy=i_dy-2;
%                     cc(i_dx,i_dy)=obj.correlationCalculator.get(d0_int+[dx dy]);
%                 end
%             end
            cc(1,2)=obj.correlationCalculator.get(d0_int+[-1 0]);
            cc(2,1)=obj.correlationCalculator.get(d0_int+[0 -1]);
            cc(2,2)=obj.correlationCalculator.get(d0_int+[0 0]);
            cc(2,3)=obj.correlationCalculator.get(d0_int+[0 1]);
            cc(3,2)=obj.correlationCalculator.get(d0_int+[1 0]);
            [c, i] = max2d_subpixel(cc);
            dx=i(1);
            dy=i(2);
            
            d = d0_int+[dx dy];
        end
    end
            
    methods(Static)
        function shifted = shift(source,d)
			if(size(d,1)<size(source,3))
				d = repmat(d(1,:),size(source,3),1);
			end
			if(size(source,3)==1 && size(d,2)==1)
				d = d';
			end
			if(all(mod(d(:),1)==0))
				shifted = zeros(size(source),class(source));
			else
				shifted = zeros(size(source),'single');
			end
			for i = 1:size(source,3)	
				if(all(mod(d(i,:),1)==0))
					shifted(:,:,i) = shift@ImageRegistrator(source(:,:,i),d(i,:));
				else
					int_shifted = shift@ImageRegistrator(source(:,:,i),fix(d(i,:)));
					shifted(:,:,i) = BilinearImageRegistrator.bilinear_frac_shift(int_shifted,d(i,:)-fix(d(i,:)));
				end
			end
        end
        
        function shifted = bilinear_frac_shift(int_shifted,d)
            assert(all(abs(d)<=1))
            h1 = BilinearImageRegistrator.bilinear_kernel(d(1)-floor(d(1)));
            h2 = BilinearImageRegistrator.bilinear_kernel(d(2)-floor(d(2)));
            frac_shifted = conv2(h1,h2,int_shifted,'full');
            
            shifted = frac_shifted((1:end-1)+(d(1)<0),(1:end-1)+(d(2)<0));
        end
        function h = bilinear_kernel(t)
            assert(numel(t) == 1);
            assert(t>=0 && t<=1);
            h = [1-t t];
        end
    end
end