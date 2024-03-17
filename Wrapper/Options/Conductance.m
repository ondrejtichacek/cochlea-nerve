classdef Conductance < Unit
    %CONDUCATANCE
    
    properties (Constant, Hidden)
        si_unit = 'S';
        internal_unit = 'nS';
        disp_unit = 'nS';
        
        internal_unit_power = -9;
        
        quantity = 'conductance';
    end
    properties (Dependent)
        PS
        TS
        GS
        MS
        kS
        S
        mS
        uS
        nS
        pS
        fS
    end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Conductance;
        end
        
    end
    
    
    methods
        function obj = Conductance(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, Conductance.si_unit, Conductance.internal_unit_power);
            
        end
        function val = get.PS(obj)
            pow = obj.internal_unit_power - 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.TS(obj)
            pow = obj.internal_unit_power - 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.GS(obj)
            pow = obj.internal_unit_power - 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.MS(obj)
            pow = obj.internal_unit_power - 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.kS(obj)
            pow = obj.internal_unit_power - 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.S(obj)
            pow = obj.internal_unit_power + 0;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.mS(obj)
            pow = obj.internal_unit_power + 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.uS(obj)
            pow = obj.internal_unit_power + 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.nS(obj)
            pow = obj.internal_unit_power + 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.pS(obj)
            pow = obj.internal_unit_power + 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.fS(obj)
            pow = obj.internal_unit_power + 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
    end
    
end

