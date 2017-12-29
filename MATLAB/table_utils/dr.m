classdef dr
    %dr  Assorted utilities

    % properties (Constant)
    %     dict = make_hash({{'BASE', 1}});
    % end

    % methods (Static)
    %     function [] = setbase(i)
    %         sprintf('before: %d\n', dr.BASE)
    %         keys(dr.dict);
    %         dr.dict('BASE') = i;
    %         sprintf('after: %d\n', dr.BASE)
    %     end

    %     function b = BASE()
    %         b = dr.dict('BASE');
    %     end
    % end

    properties (Constant)
        BASE = 1;
    end

    methods (Static)
        function [x, xs] = first(x_xs)
        %FIRST split first element off sequence.
            sh = dr.chk_vector(x_xs);
            x = dr.ith_(x_xs, dr.BASE);
            if nargout > 1
                n = prod(sh);
                xs = x_xs(dr.BASE+1:dr.BASE+n-1); %
            end
        end

        function xs = rest(x_xs)
        %REST split first element off sequence and return the rest.
            [~, xs] = dr.first(x_xs);
        end

        function [x, xs] = last(xs_x)
        %LAST split last element off sequence.
            sh = dr.chk_vector(xs_x);
            n = prod(sh);
            x = dr.ith_(xs_x, dr.BASE+n-1);
            if nargout > 1
                xs = xs_x(dr.BASE:dr.BASE+n-2); %
            end
        end

        function xs = least(xs_x)
        %LEAST split last element off sequence and return the rest.
            [~, xs] = dr.last(xs_x);
        end


        function ii = o2ii(offset, shape, colmajor, varargin)
        %O2II convert offset to indices for a given array shape, assuming
        %   specified order
            narginchk(3, 4);
            base = dr.first([varargin {dr.BASE}]);
            ii = dr.o2ii_(offset, shape, colmajor, base);
        end

        function offset = ii2o(idxs, shape, colmajor, varargin)
        %II2O convert indices to offset for a given array shape, assuming
        %   specified order.
            narginchk(3, 4);
            base = dr.first([varargin {dr.BASE}]);
            offset = dr.ii2o_(idxs, shape, colmajor, base);
        end

        function ii = traversal(shape, colmajor, varargin)
        %TRAVERSAL return sequence of indices for given array shape, assuming
        %   specified order of traversal
            narginchk(2, 3);
            base = dr.first([varargin {dr.BASE}]);
            ii = dr.traversal_(shape, colmajor, base);
        end


        function out = mkbox(shape, colmajor, outer, varargin)
        %MKBOX make expanded test nd-array.
            narginchk(3, 4);
            base = dr.first([varargin {dr.BASE}]);
            out = dr.mkbox_(shape, colmajor, outer, base);
        end

        function out = mkslab(shape, colmajor, outer, varargin)
        %MKSLAB make contracted test nd-array.
            narginchk(3, 4);
            base = dr.first([varargin {dr.BASE}]);
            out = dr.mkslab_(shape, colmajor, outer, base);
        end

        function out = box2slab(box, outer)
        %BOX2SLAB convert expanded test nd-array (box) to contracted test
        %   nd-array (slab).
            sh = size(box);
            u = dr.unroll(box, outer);
            s = cell2mat(arraymap(@(i) contract_(u(i, :)), (1:size(u, 1)).'));

            if outer
                out = reshape(s, dr.least(sh));
            else
                out = reshape(s, [1 dr.rest(sh)]);
            end
        end


        function out = slab2box(slab, outer)
        %SLAB2BOX convert contracted test nd-array (slab) to expanded test
        %   nd-array (box).

            if outer
                sh = size(slab);
            else
                sh = size(sqz1(slab));
            end

            u = slab(:).';
            n = numel(u);
            s = cell2mat(arraymap(@(i) expand_(u(i)), (1:n).'));
            d = numel(sh);
            out = reshape(s, [sh d]);

            if ~outer
                out = permute(out, [d + 1 1:d]);
            end
        end


        function out = unroll(ndarr, outer, varargin)
        %UNROLL "unroll" test nd-array.
            narginchk(2, 3);
            simple = nargin > 2 && varargin{1};
            nd = ndims(ndarr);

            if outer
                if simple
                    w = 1;
                else
                    w = size(ndarr, nd);
                end
            else
                w = size(ndarr, 1);
                if simple && w > 1
                    error('inconsistent arguments');
                end

                p = circshift(1:nd, [0 -1]);
                ndarr = permute(ndarr, p);
            end

            out = reshape(ndarr, [], w);
        end



        function i = indexof(x, y)
        %INDEXOF index of first argument in second argument.
            [~, i] = ismember(x, y);
        end

        % function o = ii2o_cm(idxs, shape, varargin)
        % %II2O_CM convert indices to offset for a given array shape, assuming
        % %   column-major order
        %     narginchk(2, 3);
        %     base = dr.first([varargin {dr.BASE}]);
        %     [ii, sh, d] = dr.chk_ii2o_args(idxs, shape, base);
        %     o = dr.BASE + dr.ii2o_cm_(ii, sh, d);
        % end

        % function o = ii2o_rm(idxs, shape, varargin)
        % %II2O_RM convert indices to offset for a given array shape, assuming
        % %   row-major order
        %     narginchk(2, 3);
        %     base = dr.first([varargin {dr.BASE}]);
        %     [ii, sh, d] = dr.chk_ii2o_args(idxs, shape, base);
        %     o = dr.BASE + dr.ii2o_rm_(ii, sh, d);
        % end




        % function out = iis_cm(shape, varargin)
        % %IIS_CM return sequence of indices for given array shape, assuming
        % %   column-major order.
        %     narginchk(1, 2);
        %     base = dr.first([varargin {dr.BASE}]);
        %     out = dr.iis_(shape, true, base);
        % end

        % function out = iis_rm(shape, varargin)
        % %IIS_RM return sequence of indices for given array shape, assuming
        % %   row-major order.
        %     narginchk(1, 2);
        %     base = dr.first([varargin {dr.BASE}]);
        %     out = dr.iis_(shape, false, base);
        % end



        % function ii = o2ii_cm(offset, shape, varargin)
        % %O2II_CM convert offset to indices for a given array shape, assuming
        % %   column-major order
        %     narginchk(2, 3);
        %     base = dr.first([varargin {dr.BASE}]);
        %     [o, sh, d] = dr.chk_o2ii_args(offset, shape, base);
        %     ii = dr.BASE + dr.o2ii_cm_(o, sh, d);
        % end

        % function ii = o2ii_rm(offset, shape, varargin)
        % %O2II_RM convert offset to indices for a given array shape, assuming
        % %   row-major order
        %     narginchk(2, 3);
        %     base = dr.first([varargin {dr.BASE}]);
        %     [o, sh, d] = dr.chk_o2ii_args(offset, shape, base);
        %     ii = dr.BASE + dr.o2ii_rm_(o, sh, d);
        % end


        function nn = vns(tbl, varargin)
            narginchk(1, 2);
            nn = tbl.Properties.VariableNames;
            if nargin == 2
                nn = nn(1, varargin{1});
            end
        end

        function nn = vnames(tbl, varargin)
            narginchk(1, 2);
            nn = tbl.Properties.VariableNames;
            function nm = i2n(i)
                if isnumeric(i)
                    nm = nn{i};
                else
                    nm = i;
                end
            end
            if nargin > 1
                vars = varargin{1};
                if isnumeric(vars)
                    vars = num2cell(vars);
                end
                nn = cellmap(@i2n, vars);
            end
        end

        function ii = vidxs(tbl, vars)
            if isempty(vars)
                ii = [];
                return
            end
            nn = varnames(tbl);

            if isstr_(vars)
                vars = strsplit(vars);
            elseif iscolumn(vars)
                vars = reshape(vars, 1, []);
            else
                assert(isrow(vars), 'second argument is not 1-dimensional');
            end

            function i = n2i(n)
                if isnumeric(n)
                    i = n;
                else
                    [b, i] = ismember(n, nn);
                    if ~b
                        i = n;
                    end
                end
            end
            if iscell(vars)
                ii = cellmap(@n2i, vars);
                unk = ii(arrayfun(@(i) ~isnumeric(ii{i}), 1:numel(ii)));
                dr.unrecognized_variables_(unk);
                ii = cell2mat(ii);
            else
                assert(isnumeric(vars), ...
                       'second argument is neither a cell array nor a numeric array');
                ii = vars;
            end

            nvs = numel(nn);
            aob = ii(ii < 1 | nvs < ii);
            if numel(aob) > 0
              error('DR20:badsubscript', ...
                    'one or more out-of-bound variable indices found: %s', ...
                    strjoin(cellmap(@int2str, num2cell(unique(aob, 'stable')))));
            end

        end

    end



    methods(Static, Access=private)

        function [o, sh, d] = chk_o2ii_args(offset, shape, base)
            [sh, d] = dr.chk_shape(shape);
            o = offset - base;
            assert(0 <= o && o < prod(sh));
        end

        function [ii, sh, d] = chk_ii2o_args(idxs, shape, base)
            [sh, d] = dr.chk_shape(shape);
            [ii, d0] = dr.chk_idxs(idxs - base);
            assert(d0 == d);
            assert(all(0 <= ii & ii < sh));
        end

        function [sh, d] = chk_shape(shape)
            [sh, d] = dr.chk_idxs(shape, true);
            assert(all(sh > 0));
        end

        function [sh] = chk_vector(v)
            assert(isrow(v) || iscolumn(v));
            sh = size(v);
        end

        function [ii, d] = chk_idxs(idxs, varargin)
            narginchk(1, 2);
            [~] = dr.chk_vector(idxs);
            d = numel(idxs);
            assert(d > 0 || nargin > 1 && varargin{1});
            ii = reshape(idxs, 1, []);
            assert(all(isint(ii)));
        end



        function ii = o2ii_(offset, shape, colmajor, base)
        %O2II_ convert linear offset to nd-array indices.
            o2ii__ = dr.mk_o2ii(shape, colmajor, base);
            ii = o2ii__(offset);
        end

        function fn = mk_o2ii(shape, colmajor, base)
        %MK_O2II make offset-to-indices converter.

            if colmajor
                split_ = @dr.first;
                join_  = @(x, y) [x y];
            else
                split_ = @dr.last;
                join_  = @(x, y) [y x];
            end

            function out = o2ii_(offset, shape, d)
                % if d == 1
                %     out = offset;
                if d == 0
                    out = reshape([], [1 0]);
                else
                    [N, Ns] = split_(shape);
                    i = mod(offset, N);
                    out = join_(i, o2ii_((offset - i)/N, Ns, d - 1));
                end
            end

            [sh, d] = dr.chk_shape(shape);
            u = base + prod(sh);
            function ii = wrapper(o)
                assert(base <= o && o < u);
                ii = base + o2ii_(o - base, sh, d);
            end

            fn = @wrapper;
        end


        function offset = ii2o_(idxs, shape, colmajor, base)
            if colmajor
                split_ = @dr.first;
            else
                split_ = @dr.last;
            end

            function out = ii2o__(idxs, shape, d)
                if d == 1
                    out = split_(idxs);
                else
                    [n, ns] = split_(idxs);
                    [N, Ns] = split_(shape);
                    out = n + N * ii2o__(ns, Ns, d - 1);
                end
            end

            [ii, sh, d] = dr.chk_ii2o_args(idxs, shape, base);
            offset = base + ii2o__(ii, sh, d);
        end


        function out = traversal_(shape, colmajor, base)
            o2ii = dr.mk_o2ii(shape, colmajor, base);
            n = prod(shape);
            out = cell2mat(arraymap(o2ii, (base:base + n - 1).'));
        end


        function box = mkbox_(shape, colmajor, outer, base)

            function out = mk_colmajor_box(sh, d, outer, base)
                traverse = dr.traversal(sh, true, base);
                n = prod(sh);
                assert(n == size(traverse, 1));
                ii = base:base + n - 1;
                traverse = arraymap(@(i) num2cell(traverse(i, :)), ii.');
                ss = base:base + d - 1;

                out = zeros(dr.shape_([sh d]), 'uint16');
                for i = ii
                    jj = traverse{i};
                    for s = ss
                        out(jj{:}, s) = jj{s};
                    end
                end

                if ~outer && d > 0
                    out = permute(out, [d+1 1:d]);
                end
            end

            [sh, d] = dr.chk_shape(shape);

            if colmajor
                box = mk_colmajor_box(sh, d, outer, base);
            else

                if outer
                    f = d + 1;
                else
                    f = 1;
                end

                % flip requires R2013b or newer; for earlier versions,
                % one can use flipdim instead of flip, but since we're
                % basing everything on the table datatype, introduced
                % also in R2013b, it's safe to assume that we also
                % have flip.
                box = flip(mk_colmajor_box(fliplr(sh), d, outer, base), f);
            end
        end


        function slab = mkslab_(shape, colmajor, outer, base)
            [sh, d] = dr.chk_shape(shape);
            if d == 0
                slab = 0;
            else
                box = dr.unroll(dr.mkbox(sh, colmajor, false, base), false);
                slab = arrayfun(@(i) contract_(box(i, :)), (1:size(box, 1)).');
                if colmajor
                    slab = reshape(slab, dr.shape_(shape));
                else
                    slab = reshape(slab, dr.shape_(fliplr(shape)));
                end
            end

            if ~outer
                slab = reshape(slab, [1 size(slab)]);
            end
        end

        function sh = shape_(shape)
            d = numel(shape);
            if d == 0 || (d == 1 && shape == 0)
                sh = [1 0];
            else
                sh = [shape 1];
            end
        end


        function out = ith_(x, i)
            if iscell(x)
                out = x{i};
            else
                out = x(i);
            end
        end

        function unrecognized_variables_(unk)
            if numel(unk) > 0
                error('DR20:UnrecognizedTableVariable', ...
                      'unrecognized variable(s): %s', ...
                      strjoin(unique(unk, 'stable'), ', '));
            end
        end

        function chk_vnames_(tbl, vnames)
            unk = setdiff(vnames, dr.vnames(tbl), 'stable');
            dr.unrecognized_variables_(unk);
        end


        % function ii = o2ii_cm_(offset, shape, d)
        %     if d == 1
        %         ii = offset;
        %     else
        %         [N, Ns] = dr.first(shape);
        %         i = mod(offset, N);
        %         ii = [i dr.o2ii_cm_((offset - i)/N, Ns, d - 1)];
        %     end
        % end

        % function ii = o2ii_rm_(offset, shape, d)
        %     if d == 1
        %         ii = offset;
        %     else
        %         [N, Ns] = dr.last(shape);
        %         i = mod(offset, N);
        %         ii = [dr.o2ii_rm_((offset - i)/N, Ns, d - 1) i];
        %     end
        % end


        % function out = range(n, varargin)
        %     narginchk(1, 3);
        %     d = 1;
        %     if nargin > 1
        %         m = n;
        %         n = varargin{1};
        %         if nargin > 2
        %             d = varargin{2};
        %         end
        %     else
        %         m = dr.BASE;
        %     end
        %     out = (m:d:n).';
        % end

        % function o = ii2o_cm_(idxs, shape, d)
        %     if d == 1
        %         o = dr.first(idxs);
        %     else
        %         [n, ns] = dr.first(idxs);
        %         [N, Ns] = dr.first(shape);
        %         o = n + N * dr.ii2o_cm_(ns, Ns, d - 1);
        %     end
        % end

        % function o = ii2o_rm_(idxs, shape, d)
        %     if d == 1
        %         o = dr.last(idxs);
        %     else
        %         [n, ns] = dr.last(idxs);
        %         [N, Ns] = dr.last(shape);
        %         o = n + N * dr.ii2o_rm_(ns, Ns, d - 1);
        %     end
        % end

    end

end
