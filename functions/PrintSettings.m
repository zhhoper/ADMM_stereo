classdef PrintSettings
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static = true)
        function val = verbose(newval)
            persistent currentval;
            if nargin >= 1
                currentval = newval;
            end
            val = currentval;
        end
    end
    
    
    
end

