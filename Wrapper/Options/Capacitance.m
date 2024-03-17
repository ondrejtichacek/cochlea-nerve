classdef Capacitance < Unit
    %CAPACITANCE
    
    properties (Constant, Hidden)
        si_unit = 'F';
        internal_unit = 'pF';
        disp_unit = 'pF';
        
        internal_unit_power = -12;
        
        quantity = 'capacitance';
    end
    properties (Dependent)
        PF
        TF
        GF
        MF
        kF
        F
        mF
        uF
        nF
        pF
        fF
    end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Capacitance;
        end
        
    end
    
    
    methods
        function obj = Capacitance(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, Capacitance.si_unit, Capacitance.internal_unit_power);
            
        end
        function val = get.PF(obj)
            pow = obj.internal_unit_power - 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.TF(obj)
            pow = obj.internal_unit_power - 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.GF(obj)
            pow = obj.internal_unit_power - 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.MF(obj)
            pow = obj.internal_unit_power - 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.kF(obj)
            pow = obj.internal_unit_power - 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.F(obj)
            pow = obj.internal_unit_power + 0;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.mF(obj)
            pow = obj.internal_unit_power + 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.uF(obj)
            pow = obj.internal_unit_power + 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.nF(obj)
            pow = obj.internal_unit_power + 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.pF(obj)
            pow = obj.internal_unit_power + 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.fF(obj)
            pow = obj.internal_unit_power + 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
    end
    
end

