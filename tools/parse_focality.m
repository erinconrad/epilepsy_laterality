function focality = parse_focality(soz)


if contains(soz,'mesial') && contains(soz,'temporal')
    focality = 'mesial temporal';
elseif contains(soz,'neocortical') && contains(soz,'temporal')
    focality = 'temporal neocortical';
elseif contains(soz,'temporal')
    focality = 'temporal undetermined';
else
    focality = 'not temporal';
end


end