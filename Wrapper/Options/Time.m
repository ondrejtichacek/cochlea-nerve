classdef Time < Unit
    %TIME

    properties (Hidden)
        disp_unit = 'ms';
    end
    
    properties (Constant, Hidden)
        si_unit = 's';
        internal_unit = 'us';
        internal_unit_power = -6;
    end
    
    properties (Dependent)
        s
        ms
        us
        ns
        ps
        fs
    end
    
%     properties (Hidden)
%         time
%     end

    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Time;
        end
        
        function frequency_unit = inverse_unit(time_unit)
            
            T = {'s',  'ms',  'us',  'ns',  'ps',  'fs'};
            F = {'Hz', 'kHz', 'MHz', 'GHz', 'THz', 'PHz'};
            
            i = find(strcmp(time_unit, T), 1);
            if i            
                frequency_unit = F{i};
            else
                error('Inverse unit not found for the unit %s', time_unit);
            end
            
        end
    end

    methods
        
        function obj = Time(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, 's', Time.internal_unit_power);
            
        end
        
        function R = times(A, B)
            
            if isa(A, 'Time') && isa(B, 'Frequency')
                a = A.(Time.si_unit);
                b = B.(Frequency.si_unit);
                
                R = a .* b;
            else
                R = times@Unit(A,B);
            end
            
        end
        function R = mtimes(A, B)
            
            if isa(A, 'Time') && isa(B, 'Frequency')
                a = A.(Time.si_unit);
                b = B.(Frequency.si_unit);
                
                R = a * b;
            else
                R = mtimes@Unit(A,B);
            end
            
        end
        
        
        function t = get.s(obj)
            pow = obj.internal_unit_power + 0;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
        function t = get.ms(obj)
            pow = obj.internal_unit_power + 3;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
        function t = get.us(obj)
            pow = obj.internal_unit_power + 6;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
        function t = get.ns(obj)
            pow = obj.internal_unit_power + 9;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
        function t = get.ps(obj)
            pow = obj.internal_unit_power + 12;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
        function t = get.fs(obj)
            pow = obj.internal_unit_power + 15;
            t = Unit.prefixMultiply(obj.value, pow); %t = obj.value * 10^pow;
        end
    end
    
end

