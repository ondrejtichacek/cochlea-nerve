classdef CircuitElement < Opt
    %CIRCUITELEMENT 
    
    properties
        name char
        type char
        unit char
        si_unit char
        node1
        node2
        dep
        dep_fun
        dep_fun_der
        dep_fun_args
        dep_fun_der_respect
        pre struct
        size_dep = 'xpos'
        size_fun = @(x) ones(size(x))
        ref_value = 1
        long_value
        special char = 'none'
        debug
    end
    properties (Hidden)
        definingProperties = { ...
        'name', ...
        'type', ...
        'unit', ...
        'si_unit', ...
        'node1', ...
        'node2', ...
        'dep', ...
        'dep_fun', ...
        'dep_fun_der', ...
        'dep_fun_args', ...
        'dep_fun_der_respect', ...
        'pre', ...
        'size_dep', ...
        'size_fun', ...
        'ref_value', ...
        'long_value', ...
        'special', ...
        }
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end
    methods
        %function val = get.definingProperties(obj)
        %    val = properties(obj);
        %end
        function val = depends_on(obj, varargin)
            val = zeros([numel(obj),numel(varargin)], 'logical');
            for i = 1:numel(obj)
                for j = 1:numel(varargin)
                    val(i,:) = any(strcmp(varargin{j}, obj(i).dep));
                end
            end
        end
%         function val = value(obj, varargin)
%             val = obj.ref_value;
%             
%             assert(numel(varargin) == numel(obj.dep), 'Mismatch in number of inputs')
%             
%             for i = 1:numel(obj.dep)
%                 val = val .* obj.dep_fun{i}(varargin{i});
%             end
%         end
        function val = fixSize(obj,varargin)
            provided_args = varargin(1:2:end-1);
            provided_values = varargin(2:2:end);
            
            ix = find(strcmp(obj.size_dep, provided_args), 1);
            if ix
                xpos = provided_values{ix};
                val = obj.size_fun(xpos);
            else
                val = 1;
            end
        end
        function val = value(obj, varargin)
            % Returns the element value after specifying values of all
            % dependencies. If provided with variables that are not in the
            % dependency list, it ignores these values.
            val = value_fun_internal(obj, false, true, {}, [], varargin{:});
            
            % This is a temporary fix that makes sure that calling
            % value('xpos', x, ...) has the right shape. In other words,
            % there is always a dependance on x: ones(size(x))
            val = val .* obj.fixSize(varargin{:});
        end
        function fun = value_fun(obj, varargin)
            fun = value_fun_internal(obj, false, false, {}, [], varargin{:});
        end
        function [fun, args] = value_fun_pre(obj)
            val = obj.ref_value;
            args = obj.pre.dep_fun.args;
            F = obj.pre.dep_fun.fun;
            fun = @(varargin) val .* F(varargin{:});
        end
        function [fun, args] = der_fun_pre(obj)
            val = obj.ref_value;
            args = obj.pre.dep_fun_der.args;
            F = obj.pre.dep_fun_der.fun;
            fun = @(varargin) val .* F(varargin{:});
        end
        function [fun, args] = der_fun(obj, derivative_args, varargin)
            if ~isempty(obj.dep_fun_der_respect)
                assert(iscell(derivative_args));
                
                [~,fun_der] = cellfun(@(x) ismember(x, obj.dep_fun_der_respect), derivative_args, 'UniformOutput', true);

                [~, val, fun_type, fun_evaluated] = value_fun_internal(obj, false, false, derivative_args, fun_der, varargin{:});
                
                fun_evaluated(fun_der) = true;
                fun_type(fun_der) = [];
                
                fun_not_evaluated = find(fun_evaluated == false);
                
                fun_type = [fun_type, repmat({'dep_fun_der'}, size(fun_der))];

                [fun, args] = create_anonymous_fcn(obj, val, fun_type, [fun_not_evaluated, fun_der]);
            else
                fun = NaN;
                args = {};
            end
        end
        function fun = der_fun_old(obj, derivative_args, varargin)
            assert(iscell(derivative_args));
            [~, val, fun_type, args] = value_fun_internal_old(obj, false, false, derivative_args, varargin{:});
            fun_type = [fun_type, repmat({'dep_fun_der'}, size(derivative_args))];
            fun = create_anonymous_fcn_old(obj, val, fun_type, [args, derivative_args]);
        end
        function [fun, val, fun_type, args] = value_fun_internal_old(obj, flag_err_arg_not_foud, flag_err_arg_not_provided, omit_args, varargin)
            val = obj.ref_value;
            
            provided_args = varargin(1:2:end-1);
            provided_values = varargin(2:2:end);
            
            assert(numel(provided_args) == numel(provided_values), 'The number of arguments does not add up to the name-value pairs')
            
            for i = 1:numel(provided_args)                
                j = find(strcmp(provided_args{i}, obj.dep),1);
                if j
                    val = val .* obj.dep_fun{j}(provided_values{i});
                else
                    if flag_err_arg_not_foud
                        error('Element value is not dependent on argument %s.', provided_args{i})
                    end
                end
            end
            
            args = setdiff(obj.dep, [provided_args, omit_args], 'stable'); % keep the order defined by obj.dep
            
            n = numel(args);
            
            if n > 0 && flag_err_arg_not_provided
                error('Insuficient number of input arguments provided. Check dependency of the value.')
            end
            
            fun_type = repmat({'dep_fun'}, size(args));
            
            fun = create_anonymous_fcn(obj, val, fun_type, args);
            
        end
        function [fun, val, fun_type, fun_evaluated] = value_fun_internal(obj, flag_err_arg_not_foud, flag_err_arg_not_provided, omit_args, omit_fun, varargin)
            val = obj.ref_value;
            
            provided_args = varargin(1:2:end-1);
            provided_values = varargin(2:2:end);
            
            assert(numel(provided_args) == numel(provided_values), 'The number of arguments does not add up to the name-value pairs')
            
            for i = 1:numel(omit_args)                
                j = find(strcmp(omit_args{i}, provided_args),1);
                if j
                    provided_args(j) = [];
                    provided_values(j) = [];
                end
            end
            
            for i = 1:numel(provided_args)                
                j = find(strcmp(provided_args{i}, obj.dep),1);
                if j
                else
                    if flag_err_arg_not_foud
                        error('Element value is not dependent on argument %s.', provided_args{i})
                    end
                end
            end
            
            fun_evaluated = zeros(size(obj.dep_fun), 'logical');
            
            for i = 1:numel(obj.dep_fun)
                if ~any(i == omit_fun)
                    if nargin(obj.dep_fun{i}) == 1
                        arg = obj.dep{i};

                        j = find(strcmp(arg, provided_args),1);
                        if j
                            val = val .* obj.dep_fun{i}(provided_values{j});
                            fun_evaluated(i) = true;
                        else
                            if flag_err_arg_not_provided
                                error('Dependent argument %s not provided.', arg)
                            end
                        end

                    else
                        arg = obj.dep_fun_args{i};
                        [~,m] = cellfun(@(x) ismember(x, provided_args), arg, 'UniformOutput', true);
                        if any(m == 0)
                            if flag_err_arg_not_provided
                                error('Dependent argument(s) not provided.')
                            end
                        else
                            val = val .* obj.dep_fun{i}(provided_values{m});
                            fun_evaluated(i) = true;
                        end
                    end
                end
            end
            
            fun_not_evaluated = find(fun_evaluated == false);
            
            fun_type = repmat({'dep_fun'}, size(fun_not_evaluated));
            
            if any(fun_evaluated == false)
                fun = create_anonymous_fcn(obj, val, fun_type, fun_not_evaluated);
            else
                fun = val;
            end
            
        end
        function [fun, args] = create_anonymous_fcn(obj, val, fun_type, fun_not_evaluated)
            
            args = {};
            
            str1 = '@(';
            str2 = 'val.*';
            
            for j = 1:numel(fun_not_evaluated)
                i = fun_not_evaluated(j);
                if nargin(obj.(fun_type{j}){i}) == 1
                    arg = {obj.dep{i}};
                else
                    arg = obj.dep_fun_args{i};
                end

                tmp = [arg; num2cell(repmat(j,size(arg)))];
                        
                sarg = sprintf('%s%d,', tmp{:});
                str1 = [str1, sarg];
                str2 = [str2, sprintf('obj.%s{%d}(%s)',fun_type{j},i,sarg(1:end-1)), '.*'];
                
                args = [args, arg];
                
            end
            
            str = ['fun = ',str1(1:end-1), ')', str2(1:end-2),';'];
            eval(str);
            
        end
        function fun = create_anonymous_fcn_old(obj, val, fun_type, args)
            
            n = numel(args);
            
            if n == 0
                fun = val;
            elseif n == 1
                j = find(strcmp(args{1}, obj.dep),1);
                fun = @(x1) val .* obj.(fun_type{1}){j}(x1);
            elseif n == 2
                j1 = find(strcmp(args{1}, obj.dep),1);
                j2 = find(strcmp(args{2}, obj.dep),1);
                fun = @(x1,x2) val .* obj.(fun_type{1}){j1}(x1) .* obj.(fun_type{2}){j2}(x2);
            else
                % this is bad practice, but I don't see other way
                str1 = '@(';
                str2 = 'val.*';
                for i = 1:n
                    j = find(strcmp(args{i}, obj.dep),1);
                    str1 = [str1, args{i}, ','];
                    str2 = [str2, sprintf('obj.%s{%d}(%s)',fun_type{i},j,args{i}), '.*'];
                end
                str = ['fun = ',str1(1:end-1), ')', str2(1:end-2),';'];
                eval(str);
            end
        end
        function obj = CircuitElement(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
%             obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
    end
    
end

