function obj = myloadobj(file)
%% myloadobj
%
%   obj = myloadobj(file)
%       Loads an object, stored as variables in a .mat file located at
%       'file', into a dcpObject
%
%%

s = load(file);
e = [s.objType '(''' s.sname ''',''' s.datapath ''')'];
obj = eval(e);  %create object
for fn = fieldnames(s)'    %enumerat fields
    try
        obj.(fn{1}) = s.(fn{1});   %and copy
    catch
        warning('Could not copy field %s', fn{1});
    end
end
