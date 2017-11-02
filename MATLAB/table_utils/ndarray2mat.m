function out = ndarray2mat( varargin )
    [nda, outer, colex, simple] = process_args(varargin{:});

    nd = ndims(nda);
    if simple
        nd = nd + 1;
        nda = shiftdim(nda, -1);
    elseif outer
        nda = permute(nda, circshift(1:nd, [0 1]));
    end
    nc = size(nda, 1);

    if ~colex
        nda = permute(nda, circshift(nd:-1:1, [0 1]));
    end

    out = reshape(nda, nc, []).';
end

function [nda, outer, colex, simple] = process_args( varargin )
    p = inputParser;
    p.addRequired('nda');
    p.addParameter('Colex', false);
    p.addParameter('Outer', false);
    p.addParameter('Simple', false);
    p.parse(varargin{:});
    args = p.Results;
    nda = args.nda;
    outer = args.Outer;
    colex = args.Colex;
    simple = args.Simple;
    assert(isequal(simple, false) || isequal(outer, true));
end
