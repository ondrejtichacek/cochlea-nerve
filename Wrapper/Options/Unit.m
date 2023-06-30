classdef Unit < matlab.mixin.Copyable & matlab.mixin.CustomDisplay
    %UNIT
    properties (Hidden, Access = protected)
        value
    end
    
    properties (Access = public)
        test_set_me
    end
    properties (Hidden, SetAccess = protected)
        test_cant_set_me
    end
    properties (Hidden, GetAccess = protected)
        test_cant_see_me = 'secret'
    end
    
    methods (Access = protected)
        function displayEmptyObject(obj)
            % Custom display of empty Unit
            
            fprintf('     []    (0x0 %s)\n', ...
                class(obj));
        end
        function displayScalarObject(obj)
            % Custom display of scalar-valued Unit
            
            fprintf('     %g %s    (1x1 %s)\n', ...
                obj.(obj.disp_unit), ...
                obj.disp_unit, ...
                class(obj));
        end
        function displayNonScalarObject(obj)
            % Custom display of nonscalar Unit
            v = obj.(obj.disp_unit);
            n = numel(v);
            k = size(v);
            sz = sprintf('%dx', k);
            sz = sz(1:end-1);
            if numel(v) <= 10
                s = sprintf('%g, ', v);
                s = sprintf('[ %s ]', s(1:end-2));
                fprintf('     %s %s    (%s %s)\n', ...
                    s, ...
                    obj.disp_unit, ...
                    sz, ...
                    class(obj));
            else
                fprintf('     [%s %s]\n', ...
                    sz, ...
                    class(obj));
            end
        end
    end
    methods (Static)
        function vars = batch_convert(unit_info, vars)
            vars = cellfun(@(v) Unit.batch_convert_single(unit_info, v), vars, 'UniformOutput', false);
        end
        function var = batch_convert_single(unit_info, var)
            units = fieldnames(unit_info);
            i = find(strcmp(class(var), units));
            if ~isempty(i)
                var = var.(unit_info.(units{i}));
            end
        end
        function v = toSI(v)
            if isa(v, 'Unit')
                v = v.(v.si_unit);
            end
        end
        function k = conversionConstant(from, to)
            % returns pow such that
            % T = k * F
            % where F [from unit] and T [to unit]
            
            k = 10^Unit.conversionPower(from, to);
        end
        function pow = conversionPower(from, to)
            % returns pow such that
            % T = F * 10^pow
            % where F [from unit] and T [to unit]
            
            % assume unit prefix is one letter
            
            if strcmp(from, to)
                pow = 0;
                return
            end
            
            b1 = from(2:end);
            b2 = to(2:end);
            
            n1 = numel(from);
            n2 = numel(to);
            
            if n1 < n2
                b1 = from;
            elseif n2 < n1
                b2 = to;
            end
                        
            assert(strcmp(b1, b2), 'Base unit must be identical');
            
            pow = Unit.handlePowerPrefixes(from, b1) ...
                - Unit.handlePowerPrefixes(to, b2);
        end
        function pow = handlePowerPrefixes(unit, base)
            switch unit
                case ['P', base], pow =  15;
                case ['T', base], pow =  12;
                case ['G', base], pow =   9;
                case ['M', base], pow =   6;
                case ['k', base], pow =   3;
                case base,        pow =   0;
                case ['d', base], pow =  -1;
                case ['c', base], pow =  -2;
                case ['m', base], pow =  -3;
                case ['u', base], pow =  -6;
                case ['n', base], pow =  -9;
                case ['p', base], pow = -12;
                case ['f', base], pow = -15;
                otherwise, error('Unit %s not recognised.', unit);
            end
        end
        function [a, b, unit, constructor] = binaryOperatorCommon(A, B)
            assert(isa(A,class(B)))
            
            unit = A.internal_unit;
            
            if numel(A) == 1
                a = A.(unit);
            else
                a = reshape([A.(unit)], size(A));
            end
            if numel(B) == 1
                b = B.(unit);
            else
                b = reshape([B.(unit)], size(B));
            end
            
            constructor = A.constructor();
            
            if strcmp(func2str(constructor),'Quantity')
                constructor = @(v, u) constructor(v, u, A.disp_unit);
            end
        end
        function [a, unit, constructor] = unaryOperatorCommon(A)
            
            unit = A.internal_unit;
            
            if numel(A) == 1
                a = A.(unit);
            else
                a = reshape([A.(unit)], size(A));
            end
            
            constructor = A.constructor();
            if strcmp(func2str(constructor),'Quantity')
                constructor = @(v, u) constructor(v, u, A.disp_unit);
            end
        end
        function [a, b] = binaryComparisonCommon(A, B)
            if isnumeric(A) && all(A == 0)
                a = A;
                [b, ~, ~] = Unit.unaryOperatorCommon(B);
            elseif isnumeric(B) && all(B == 0)
                b = B;
                [a, ~, ~] = Unit.unaryOperatorCommon(A);
            else
                [a, b, ~, ~] = Unit.binaryOperatorCommon(A,B);
            end
        end
        function u = simplifyUnit(unit)
            switch unit
                    case '1/s'
                        u = 'Hz';
                    case '1/Hz'
                        u = 's';
                    otherwise
                        u = unit;
            end
        end
        function con = getConstructor(unit)
            switch unit
                    case 's'
                        con = @Time;
                    case 'Hz'
                        con = @Frequency;
                    otherwise
                        con = @Quantity;
            end
        end
        function R = newQuantity(r, unit, displ_unit)
            if nargin < 3
                displ_unit = unit;
            end
            simpl_unit = Unit.simplifyUnit(displ_unit);
            constructor = Unit.getConstructor(simpl_unit);
            if strcmp(func2str(constructor), 'Quantity')
                R = constructor(r, unit, displ_unit);
            else
                R = constructor(r, simpl_unit);
            end
        end
        function val = prefixMultiply(val, pow)
            if pow ~= 0
                val = val .* 10^pow;
            end
        end
        function varargout = min_max(FUN, A, varargin)
            nout = max([1, nargout]);
            varargout = cell(1, nout);
            
            if numel(varargin) >= 1
                if isempty(varargin{1})
                    fun = 1;
                else
                    fun = 2;
                end
            else
                fun = 1;
            end
            
            if fun == 1
                [a, u, constructor] = Unit.unaryOperatorCommon(A);
                [varargout{:}] = FUN(a, varargin{:});
            elseif fun == 2
                B = varargin{1};
                [a, b, u, constructor] = Unit.binaryOperatorCommon(A,B);
                [varargout{:}] = FUN(a,b,varargin{2:end});
            end
            
            varargout{1} = constructor(varargout{1}, u);
        end
    end
    
    methods
        function hash = uint8(obj)
            
            defining_data = struct( ...
                'class', class(obj), ...
                'internal_unit', obj.internal_unit, ...
                'value', obj.value);
            
            OPT.Method = 'SHA-1';
            OPT.Format = 'uint8';
            hash = DataHash(defining_data, OPT);
            
        end
        
        % -----------------------------------------------------------------
        % CONSTRUCTOR
        function obj = Unit(value, unit, base_unit, internal_unit_power)
            
            assert(isnumeric(value), sprintf('Input must be of a numeric class, not %s.', class(value)))
            
            if isempty(unit)
                % not specifying a unit is allowed if the value is 0 or Inf
                if any(value(:) ~= 0 & ~isinf(value(:)) & ~isnan(value(:)))
                    error('Unit not set')
                end
            else
                
                pow = Unit.handlePowerPrefixes(unit, base_unit);                
                pow = pow - internal_unit_power;
                
                if pow ~= 0
                    value = value * 10^pow;
                end
                
            end            
            obj.value = value;
        end
        % -----------------------------------------------------------------
        function n = numArgumentsFromSubscript(obj, s, indexingContext)
        %length(s(1).subs{:});
            switch indexingContext
                case matlab.mixin.util.IndexingContext.Statement
                    n = 1; % nargout for indexed reference used as statement
                case matlab.mixin.util.IndexingContext.Expression
                    n = 1; % nargout for indexed reference used as function argument
                case matlab.mixin.util.IndexingContext.Assignment
                    n = 1; % nargin for indexed assignment
            end
        end
        % -----------------------------------------------------------------
        % OVERLOADED OPERATORS
        % A + B
        function R = plus(A, B)
            [a, b, u, constructor] = Unit.binaryOperatorCommon(A,B);
            R = constructor(a + b, u);
        end
        % A - B
        function R = minus(A, B)
            [a, b, u, constructor] = Unit.binaryOperatorCommon(A,B);
            R = constructor(a - b, u);
        end
        % +A
        function R = uplus(A)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            R = constructor(+a, u);
        end
        % -A
        function R = uminus(A)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            R = constructor(-a, u);
        end
        % A .* B
        function R = times(A, B)
            
            % (number A) .* B
            if isnumeric(A) || islogical(A)
                [b, u, constructor] = Unit.unaryOperatorCommon(B);
                
                R = constructor(A.*b, u);
                
            % A .* (number B)
            elseif isnumeric(B) || islogical(B)
                R = B .* A;
                
            % A .* B
            else
                unit_a = A.si_unit;
                unit_b = B.si_unit;
                
                a = reshape([A.(unit_a)], size(A));
                b = reshape([B.(unit_b)], size(B));
                
                r = a .* b;
                
                R = Unit.newQuantity(r, sprintf('%s*%s', unit_a, unit_b));
            end
        end
        % A * B
        function R = mtimes(A, B)
            
            % (number A) * B
            if isnumeric(A) || islogical(A)
                [b, u, constructor] = Unit.unaryOperatorCommon(B);
                
                R = constructor(A*b, u);
                
            % A * (number B)
            elseif isnumeric(B) || islogical(B)
                R = B * A;
                
            % A * B
            else
                unit_a = A.si_unit;
                unit_b = B.si_unit;
                
                a = reshape([A.(unit_a)], size(A));
                b = reshape([B.(unit_b)], size(B));
                
                r = a * b;
                
                R = Unit.newQuantity(r, sprintf('%s*%s', unit_a, unit_b));
            end
        end
        % A ./ B
        function R = rdivide(A, B)
            
            % A ./ (number B)
            if isnumeric(B) || islogical(B)
                [a, u, constructor] = Unit.unaryOperatorCommon(A);
                
                R = constructor(a./B, u);
                
            else
                % (number A) ./ B
                if isnumeric(A) || islogical(A)
                    unit_a = '1';
                    a = A;
                    
                % A ./ B
                else
                    unit_a = A.si_unit;
                    a = reshape([A.(unit_a)], size(A));
                end
                
                unit_b = B.si_unit;
                b = reshape([B.(unit_b)], size(B));
                
                r = a ./ b;
                
                if isa(A,class(B))
                    R = r;
                else
                    R = Unit.newQuantity(r, sprintf('%s/%s', unit_a, unit_b));
                end
            end
        end
        % A / B
        function R = mrdivide(A, B)
            
            % A / (number B)
            if isnumeric(B) || islogical(B)
                [a, u, constructor] = Unit.unaryOperatorCommon(A);
                
                R = constructor(a/B, u);
                
            else
                % (number A) / B
                if isnumeric(A) || islogical(A)
                    unit_a = '1';
                    a = A;
                    
                % A / B
                else
                    unit_a = A.si_unit;
                    a = reshape([A.(unit_a)], size(A));
                end
                
                unit_b = B.si_unit;
                b = reshape([B.(unit_b)], size(B));
                
                r = a / b;
                
                if isa(A,class(B))
                    R = r;
                else
                    R = Unit.newQuantity(r, sprintf('%s/%s', unit_a, unit_b));
                end
            end
        end
        % A.^B
        function R = power(A, b)
            
            assert(isnumeric(b))

            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            
            if all(b == 1)
                R = constructor(a.^b, u);
            else
                R = Unit.newQuantity(a.^b, sprintf('pow__%s__%g', u, b), sprintf('%s^%g', u, b));
            end
        end
        % A^B
        function R = mpower(A, b)
            
            assert(isnumeric(b))

            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            
            if all(b == 1)
                R = constructor(a^b, u);
            else
                R = Unit.newQuantity(a^b, sprintf('%s^%g', u, b));
            end
        end
        % A < B
        function r = lt(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a < b;
        end
        % A > B
        function r = gt(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a > b;
        end
        % A <= B
        function r = le(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a <= b;
        end
        % A >= B
        function r = ge(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a >= b;
        end
        % A ~= B
        function r = ne(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a ~= b;
        end
        % A == B
        function r = eq(A, B)
            [a, b] = Unit.binaryComparisonCommon(A, B);
            r = a == b;
        end
        % A'
        function R = ctranspose(A)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            R = constructor(ctranspose(a), u);
        end
        % A.'
        function R = transpose(A)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            R = constructor(transpose(a), u);
        end
        % [A, B]
        function C = horzcat(varargin)
            A = varargin{1};
            unit = A.internal_unit;
            constructor = A.constructor();
            
            values = cellfun(@(x) x.(unit), varargin, 'UniformOutput', false);
            C = constructor(horzcat(values{:}), unit);
        end
        % [A; B]
        function C = vertcat(varargin)
            A = varargin{1};
            unit = A.internal_unit;
            constructor = A.constructor();
            
            values = cellfun(@(x) x.(unit), varargin, 'UniformOutput', false);
            C = constructor(vertcat(values{:}), unit);
        end
        % A(end)
        function ind = end(obj, k, n)
            % k is the index of the expression containing end
            % n is the total number of indices in the expression
            
            % S = size(obj);
            % if numel(S) ~= n
            %     error('Use all indices when indexing to an instance of %s', class(obj))
            % end
            
            ind = builtin('end', obj.value, k, n);
        end
        % -----------------------------------------------------------------
        % Overloading basic builtin functions
        %
        % size(A)
        function varargout = size(obj, varargin)
            nout = max([1, nargout]);
            varargout = cell(1, nout);
            [varargout{:}] = size([obj.value], varargin{:});
        end
        % numel(A)
        function val = numel(obj)
            val = numel([obj.value]);
        end
        % length(A)
        function val = length(obj)
            val = length([obj.value]);
        end
        % isinf(A)
        function ind = isinf(obj)
            ind = isinf(obj.value);
        end
        % isnan(A)
        function ind = isnan(obj)
            ind = isnan(obj.value);
        end
        % diff(A)
        function B = diff(A, varargin)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            B = constructor(diff(a, varargin{:}), u);
        end
        % abs(A)
        function B = abs(A)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            B = constructor(abs(a), u);
        end
        % min(A)
        function varargout = min(A, varargin)
            nout = max([1, nargout]);
            varargout = cell(1, nout);
            
            [varargout{:}] = Unit.min_max(@min, A, varargin{:});
        end
        % max(A)
        function varargout = max(A, varargin)
            nout = max([1, nargout]);
            varargout = cell(1, nout);
            
            [varargout{:}] = Unit.min_max(@max, A, varargin{:});
        end
        % sum(A)
        function B = sum(A, varargin)
            [a, u, constructor] = Unit.unaryOperatorCommon(A);
            B = constructor(sum(a, varargin{:}), u);
        end
        % double(A)
        function a = double(A, flag)
            if nargin == 1 || ~ strcpm(flag, 'nowarn')
                warning('Implicitely converting %s to double. Turn off by providing ''nowarn'' flag as a second argument', class(A))
            end
            a = [A.value];
        end
        function a = char(A)
            u = A.disp_unit;
            a = sprintf('%g %s', [A.(u)], u);
        end
        function C = linspace(A, B, varargin)
            % for larger arrays, this is much faster
            [a, b, u, constructor] = Unit.binaryOperatorCommon(A,B);
            C = constructor(linspace(a,b,varargin{:}), u);
        end
        % -----------------------------------------------------------------
        % Overloading SUBSREF
        % http://www.mathworks.com/help/matlab/matlab_oop/code-patterns-for-subsref-and-subsasgn-methods.html
        %
        function varargout = subsref(obj, s)
            
            nout = max([1, nargout]);
            varargout = cell(1, nout);
            
            switch s(1).type
                case '.'
                    if length(s) == 1
                        % Implement obj.PropertyName
                        try
                            [varargout{:}] = obj.protectPropertyGetAccess(s);
                            % varargout = {obj.protectPropertyGetAccess(s)};
                        catch ME
                            switch ME.identifier
                                case 'MATLAB:maxlhs'
                                    obj.protectPropertyGetAccess(s);
                                    varargout = cell(0, 0);
                                otherwise
                                    rethrow(ME)
                            end
                        end
                        
                    elseif length(s) == 2 && strcmp(s(2).type,'()')
                        % Implement obj.PropertyName(indices)
                        try
                            [varargout{:}] = obj.protectPropertyGetAccess(s);
                            % varargout = {obj.protectPropertyGetAccess(s)};
                        catch ME
                            switch ME.identifier
                                case 'MATLAB:maxlhs'
                                    obj.protectPropertyGetAccess(s);
                                    varargout = cell(0, 0);
                                otherwise
                                    rethrow(ME)
                            end
                        end
                        
                    else
                        [varargout{:}] = builtin('subsref', obj, s);
                    end
                    
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        [a, u, constructor] = Unit.unaryOperatorCommon(obj);
                        [varargout{:}] = constructor(subsref(a, s), u);
                        
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj(ind).PropertyName
                        
                        tmp = obj.protectPropertyGetAccess(s(2));
                        tmp = builtin('subsref', tmp, s(1));
                        
                        % ??? This might need redoing
                        if nargout > 1
                            if numel(tmp) == nargout
                                if isnumeric(tmp)
                                    varargout = num2cell(tmp);
                                else
                                    error('Not yet implemented')
                                end
                            else
                                error('Don''t know what to do - sizes don''t match')
                            end
                        else
                            varargout = {tmp};
                        end
                            
                    elseif length(s) == 3 && strcmp(s(2).type, '.') && strcmp(s(3).type, '()')
                        % Implement obj(indices).PropertyName(indices)
                        error('Not implemented yet');
                    else
                        % Use built-in for any other expression
                        varargout = {builtin('subsref', obj, s)};
                    end
                case '{}'
                    if length(s) == 1
                        % Implement obj{indices}
                        error('Not implemented yet');
                    elseif length(s) == 2 && strcmp(s(2).type, '.')
                        % Implement obj{indices}.PropertyName
                        error('Not implemented yet');
                    else
                        % Use built-in for any other expression
                        varargout = {builtin('subsref', obj, s)};
                    end
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        % -----------------------------------------------------------------
        % Overloading SUBSASGN
        % http://www.mathworks.com/help/matlab/matlab_oop/code-patterns-for-subsref-and-subsasgn-methods.html
        %
        function obj = subsasgn(obj, s, varargin)
            switch s(1).type
                case '.'
                    if length(s) == 1
                        % Implement obj.PropertyName = varargin{:};
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    elseif length(s) == 2 && strcmp(s(2).type, '()')
                        % Implement obj.PropertyName(indices) = varargin{:};
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    else
                        % Call built-in for any other case
                        obj = builtin('subsasgn', obj, s, varargin);
                        
                    end
                    
                case '()'
                    if length(s) == 1
                        % Implement obj(indices) = varargin{:};                        
                        obj.value = builtin('subsasgn', obj.value, s, varargin{:}.value);
                        
                    elseif length(s) == 2 && strcmp(s(2).type, '.')
                        % Implement obj(indices).PropertyName = varargin{:};
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    elseif length(s) == 3 && strcmp(s(2).type, '.') && strcmp(s(3).type, '()')
                        % Implement obj(indices).PropertyName(indices) = varargin{:};
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    else
                        % Use built-in for any other expression
                        obj = builtin('subsasgn', obj, s, varargin);
                        
                    end
                    
                case '{}'
                    if length(s) == 1
                        % Implement obj{indices} = varargin{:}
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    elseif length(s) == 2 && strcmp(s(2).type, '.')
                        % Implement obj{indices}.PropertyName = varargin{:}
                        obj = protectPropertySetAccess(obj, s, varargin{:});
                        
                    else
                        % Use built-in for any other expression
                        obj = builtin('subsasgn', obj, s, varargin);
                        
                    end
                    
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        % -----------------------------------------------------------------
        function varargout = protectPropertyGetAccess(obj, s)
            varargout = cell(1, nargout);
            
            mp = findprop(obj, s(1).subs);
            if isempty(mp)
                propery_not_found = true;
            else
                switch mp.GetAccess
                    case {'private', 'protected'}
                        error('Property `%s` is %s (GetAccess).', s.subs, mp.GetAccess)
                    case 'public'
                        [varargout{:}] = obj.(s.subs);
                        % varargout = {obj.(s.subs)};
                    otherwise
                        error('Unknown GetAccess value `%s`.', mp.GetAccess)
                end
                return
            end
            
            mp = methods(obj);
            if any(strcmp(mp, s(1).subs))
                if nargout == 0
                   builtin('subsref',obj,s);
                else
                   [varargout{:}] = builtin('subsref',obj,s);
                end
                return
            else
                method_not_found = true;
            end
            
            if propery_not_found && method_not_found
               error('No public property / method `%s` exists for class `%s`.', s.subs, class(obj))
            end
        end
        % -----------------------------------------------------------------
        function obj = protectPropertySetAccess(obj, s, varargin)
            
             mp = findprop(obj, s.subs);
             if isempty(mp)
                 error('No public property `%s` exists for class `%s`.', s.subs, class(obj))
             else
                 switch mp.SetAccess
                     case {'private', 'protected'}
                         error('Property `%s` is %s (SetAccess).', s.subs, mp.SetAccess)
                     case 'public'
                         obj.(s.subs) = varargin{:};
                     otherwise
                         error('Unknown SetAccess value `%s` (SetAccess).', mp.SetAccess)
                 end
             end
        end
        % =================================================================
        % TEST FUNCTIONS
        %
        function test_method_no_arg_no_return(obj)
            fprintf('AHOJ %s\n', obj.internal_unit)
        end
        function t = test_method_no_arg_one_return(obj)
            t = rand(4);
        end
        function [t, r] = test_method_no_arg_two_return(obj)
            t = rand(4);
            r = {'asd','rew'};
        end
        function test_method_one_arg_no_return(obj, arg)
            fprintf('AHOJ %s %d\n', obj.internal_unit, arg)
        end
        function t = test_method_one_arg_one_return(obj, arg)
            t = rand(arg);
        end
        function [t, r] = test_method_one_arg_two_return(obj, arg)
            t = rand(arg);
            r = {'asd','rew'};
        end
        function test_method_two_arg_no_return(obj, arg1, arg2)
            fprintf('AHOJ %s %d\n', obj.internal_unit, arg1, arg2)
        end
        function t = test_method_two_arg_one_return(obj, arg1, arg2)
            t = rand(arg1, arg2);
        end
        function [t, r] = test_method_two_arg_two_return(obj, arg1, arg2)
            t = rand(arg1, arg2);
            r = {'asd','rew'};
        end
        % =================================================================
    end
    
    
end

