classdef Voltage < Unit
    %VOLTAGE
    
    properties (Constant, Hidden)
        si_unit = 'V';
        internal_unit = 'nV';
        disp_unit = 'mV';
        
        internal_unit_power = -9;
        
        quantity = 'voltage';
    end
    properties (Dependent)
        PV
        TV
        GV
        MV
        kV
        V
        mV
        uV
        nV
        pV
        fV
    end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(obj)
            proplist = {obj.disp_unit};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
    
    methods (Static)
        
        function con = constructor()
            con = @Voltage;
        end
        
    end
    
    
    methods
        function obj = Voltage(value, unit)
            
            if nargin == 1
                unit = '';
            end
            
            obj = obj@Unit(value, unit, Voltage.si_unit, Voltage.internal_unit_power);
            
        end
        function val = get.PV(obj)
            pow = obj.internal_unit_power - 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.TV(obj)
            pow = obj.internal_unit_power - 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.GV(obj)
            pow = obj.internal_unit_power - 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.MV(obj)
            pow = obj.internal_unit_power - 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.kV(obj)
            pow = obj.internal_unit_power - 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.V(obj)
            pow = obj.internal_unit_power + 0;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.mV(obj)
            pow = obj.internal_unit_power + 3;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.uV(obj)
            pow = obj.internal_unit_power + 6;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.nV(obj)
            pow = obj.internal_unit_power + 9;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.pV(obj)
            pow = obj.internal_unit_power + 12;
            val = Unit.prefixMultiply(obj.value, pow);
        end
        function val = get.fV(obj)
            pow = obj.internal_unit_power + 15;
            val = Unit.prefixMultiply(obj.value, pow);
        end
    end
    
end

