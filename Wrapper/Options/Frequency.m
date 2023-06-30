classdef Frequency < Unit
    %FREQUENCY
    
    properties (Constant, Hidden)
        si_unit = 'Hz';
        internal_unit = 'Hz';
        internal_unit_power = 0;
        disp_unit = 'kHz';
    end
    
    properties (Dependent)
        Hz
        kHz
        MHz
        GHz
    end
    
%     properties (Hidden)
%         frequency
%     end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Frequency;
        end
        
        function time_unit = inverse_unit(frequency_unit)
            
            T = {'s',  'ms',  'us',  'ns',  'ps',  'fs'};
            F = {'Hz', 'kHz', 'MHz', 'GHz', 'THz', 'PHz'};
            
            i = find(strcmp(frequency_unit, F), 1);
            if i            
                time_unit = T{i};
            else
                error('Inverse unit not found for the unit %s', frequency_unit);
            end
            
        end
    end
    
    
    methods
        function obj = Frequency(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, 'Hz', Frequency.internal_unit_power);
            
        end
        
        
        function R = times(A, B)
            
            if isa(A, 'Frequency') && isa(B, 'Time')
                a = A.(Frequency.si_unit);
                b = B.(Time.si_unit);
                
                R = a .* b;
            else
                R = times@Unit(A,B);
            end
            
        end
        function R = mtimes(A, B)
            
            if isa(A, 'Frequency') && isa(B, 'Time')
                a = A.(Frequency.si_unit);
                b = B.(Time.si_unit);
                
                R = a .* b;
            else
                R = mtimes@Unit(A,B);
            end
            
        end
        
        function f = get.Hz(obj)
            pow = obj.internal_unit_power + 0;
            f = obj.value * 10^pow;
        end
        function f = get.kHz(obj)
            pow = obj.internal_unit_power - 3;
            f = obj.value * 10^pow;
        end
        function f = get.MHz(obj)
            pow = obj.internal_unit_power - 6;
            f = obj.value * 10^pow;
        end
        function f = get.GHz(obj)
            pow = obj.internal_unit_power - 9;
            f = obj.value * 10^pow;
        end
    end
    
end

