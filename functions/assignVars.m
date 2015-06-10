function assignVars(vars)
% This function takes a struct called 'vars', and assigns their contents
% as variables in the calling envorinment. E.g. If vars contains a field
% called 'idx' that equals 3, then a variable by the name 'idx' with
% contents 3 would appear in the calling environment.

varList = fieldnames(vars);

for i=1:numel(varList)
    assignin('caller',varList{i},vars.(varList{i}));
    
end

end