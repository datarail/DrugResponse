function new_name = MATLABsafename(name)


% remove the spaces 
new_name = name(name~=' ');

% convert nonalphanumeric characters to underscores
for i=1:length(new_name)
    if ~isalpha_num(new_name(i)), new_name(i) = '_'; end
end
