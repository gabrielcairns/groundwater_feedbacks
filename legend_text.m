function leg_text = legend_text(varname,vars,units)

% Round to two significant figures 
vars = round(vars, 2,'significant');

% Make sure they are column vectors
if size(vars,2) ~= 1;     vars = vars';   end

% Convert them to character arrays
if any( isbetween(vars,0,0.001,"open") | (vars>1000))
    % If we're dealing with small or large numbers, use exponential form
    varchars = num2str(vars,'%.1e');
    % But not for zero:
    if any(vars==0)
        varchars(vars==0,:) = '      0';
    end
else
    % If numbers are 'O(1)', just use default settings
    varchars = num2str(vars);
end

% Convert into string array
vars = convertCharsToStrings(cellstr(varchars));

% Change standard form notation to something more tex-friendly
for i = 1:length(vars)
    vars{i} = float_to_tex(vars{i});
end

% Concatenate
leg_text = strcat( '$', varname, '=', vars ,'$', {' '}, units);

end