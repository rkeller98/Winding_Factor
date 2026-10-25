function require_analysis_file(file,preparation)
% Give a missing input its full path and the concrete preparation step.
if ~isfile(file)
    error('WindingFactor:MissingInput','Missing input: %s\nRequired preparation: %s',file,preparation);
end
end
