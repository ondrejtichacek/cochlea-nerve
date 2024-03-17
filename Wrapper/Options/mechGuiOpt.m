classdef mechGuiOpt < Opt
    %MECHGUIOPT 
    
    properties
        H_FIG
        y0
        abs_tol
        ax_sig
        ax
        ax_twa
        line_sig
        line_0
        lines
        text_label
        patches
        line_twa
        last_t
        last_draw = -Inf;
        size_info
        Numstacks
        uicontrol
        loader_t_index
        loader_t_click
        draw_period = Time(0)
        aborted = false
        paused = false
        draw_patches (1,1) logical
        preloaded = false
        preloaded_vars = struct()
    end
    properties (Hidden)
        definingProperties = {};
    
        requiredParameters = {};
    end
    
    methods
        
        function obj = mechGuiOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
        end
        
    end
    
end

