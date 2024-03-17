classdef Current < Unit
    %CURRENT
    
    properties (Constant, Hidden)
        si_unit = 'A';
        internal_unit = 'pA';
        disp_unit = 'pA';
        
        internal_unit_power = -12;
        
        quantity = 'current';
    end
    properties (Dependent)
        PA
        TA
        GA
        MA
        kA
        A
        mA
        uA
        nA
        pA
        fA
    end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Current;
        end
        
    end
    
    
    methods
        function obj = Current(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, Current.si_unit, Current.internal_unit_power);
            
        end
        function val = get.PA(obj)
            pow = obj.internal_unit_power - 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.TA(obj)
            pow = obj.internal_unit_power - 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.GA(obj)
            pow = obj.internal_unit_power - 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.MA(obj)
            pow = obj.internal_unit_power - 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.kA(obj)
            pow = obj.internal_unit_power - 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.A(obj)
            pow = obj.internal_unit_power + 0;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.mA(obj)
            pow = obj.internal_unit_power + 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.uA(obj)
            pow = obj.internal_unit_power + 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.nA(obj)
            pow = obj.internal_unit_power + 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.pA(obj)
            pow = obj.internal_unit_power + 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.fA(obj)
            pow = obj.internal_unit_power + 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
    end
    
end

