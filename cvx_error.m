function lines = cvx_error( errmsg, widths, useline, prefix, chop )

% CVX_ERROR   Formats text for inclusion in error messages.
%    This is an internal function used by CVX. It needed to be in the CVX
%    home directory so that it's available during a fresh installation.

if isa( errmsg, 'MException' ),
    errmsg = getReport( errmsg, 'extended', 'hyperlinks', 'off' );
end
lines = {};
errmsg = strtrim( errmsg );
while ~isempty( errmsg ),
    width = widths(1);
    emax = length(errmsg);
    chunk = errmsg(1:min(width+1,emax));
    sndx = find( chunk == 10, 1, 'first' );
    if isempty( sndx ),
        if emax <= width,
            sndx = emax + 1;
        else
            sndx = find( chunk == 32, 1, 'last' );
            if isempty( sndx ),
                sndx = find( errmsg == 32, 1, 'first' );
            end
            if isempty( sndx ),
                sndx = emax + 1;
            end
        end
    end
    if sndx > 1,
        lines{1,end+1} = errmsg(1:sndx-1); %#ok
        if length(widths)>1, widths(1) = []; end
    end
    errmsg(1:min(sndx,emax)) = [];
end
if nargin >= 3 && ( ischar(useline) || useline ),
    line = '-';
    line = line(1,ones(1,max(cellfun(@length,lines))));
    sline = line;
    if ischar( useline ),
        sline(1:length(useline)) = useline;
    end
    lines = [ sline, lines, line ];
end
if nargout == 0 || nargin >= 4,
    if nargin < 4, prefix = ''; end
    lines = sprintf( [ prefix, '%s\n' ], lines{:} );
end
if nargin >= 5 && chop,
    lines(end) = [];
end
if nargout == 0,
    fprintf( lines );
    clear lines
end