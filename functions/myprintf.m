function [] = myprintf(varargin)
%myprintf calles fprintf only if PrintSettings.verbose >0

if PrintSettings.verbose == 1
    fprintf(1,varargin{:});
end


end

