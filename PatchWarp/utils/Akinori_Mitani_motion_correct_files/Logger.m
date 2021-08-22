classdef Logger < handle
	properties(SetAccess=immutable)
		started_time
		started_cputime
                pid
	end
	properties(SetAccess=private)
		last_time
		last_cputime
	end
	methods
		function obj=Logger()
			obj.started_time = clock();
			obj.started_cputime = cputime();
			obj.last_time = clock();
			obj.last_cputime = cputime();
                        obj.pid = feature('getpid');
		end
		function newline(obj,str,varargin)
			fprintf('[%05d]%s|tot[%5.0fs/%5.0fs] seg[%3.0fs/%3.0fs] %s\n',obj.pid,datestr(clock()), ...
				cputime()-obj.started_cputime,etime(clock(),obj.started_time),...
                cputime()-obj.last_cputime, etime(clock(),obj.last_time),...
				 sprintf(str,varargin{:}) )
			obj.last_time = clock();
			obj.last_cputime = cputime();
		end
	end
end
