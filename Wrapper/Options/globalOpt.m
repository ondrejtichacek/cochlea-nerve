classdef globalOpt < Opt
    properties
        username char
        computer char
        
        IDprefix
        cochleadir char
        cacheDir char
        resDir char
        tmpDir char
        bpDir char
        outStream
        
        scratchdir char
        matlabdir char
        
        allowedFolders cell = {}
        
        MNADir char
        
        num_cores int32
    end
    properties (Dependent, Hidden)
        definingProperties = {};
    end
    methods
        function opt = globalOpt( USERNAME, COMPUTERNAME )
        %GLOBALOPT Returns global options for cochlear model based on user and coputer name
        %   Options are based on in code settings
        %   IN:
        %       pocitac     = (string, optional) name of the computer
        %       uzivatel    = (string, optional) name of the user
        %                       if either argument is omitted, computer name is generated automatically by userAndComputerName.m
        %
        %   OUT:
        %       opt.computer    = input parameter
        %       opt.cochleadir  = directory containing CochleaNerve toolbox
        %       opt.MNADir      = directory containing MNA functions
        %       opt.resDir      = directory to save results
        %       opt.tmpDir      = directory to save temporary files
        %       opt.scratchdir  = directory to scratch drive (custom, TODO)
        %
        %       opt.outStream   = '_path_to_file_' file to save standard output
        %                       = 1 print on screen
        %
        %   Examples:
        %   opt = GLOBALOPT() returns settings for current computer
        %   opt = GLOBALOPT('Ondra') returns settings for current computer and user Ondra
        %   opt = GLOBALOPT([],'metacentrum') returns settings for metacentrum grid and current user
        %   opt = GLOBALOPT('tichacek','magnesium') returns settings for magnesium cluster and user tichacek
        %
        %   Comments:
        %       * opt.tmpDir needs to be at fast disc (such as scratch disc)
        %       * opt.resDir does not need to be at fast disc
        %       * the two comments above are irrelevant on desktop pc
        %


        %% Process inputs

        [opt.username, opt.computer] = userAndComputerName();

        global computername
        if nargin < 2 && (exist('computername','var') && ~isempty(computername))
            opt.computer = computername;
        end

        if nargin > 0 && ~isempty(USERNAME)
                opt.username = USERNAME;
        end
        if nargin > 1 && ~isempty(COMPUTERNAME)
                opt.computer = COMPUTERNAME;
        end

        %% Default options

        try
            %opt.matlabdir = userpath;
            %opt.matlabdir = opt.matlabdir(1:end-1); %Remove (semi)colon at end of userpath
        catch ME
            switch ME.identifier
                case 'Octave:undefined-function'
                    opt.matlabdir = '~/Box/MATLAB';
                otherwise
                    rethrow(ME)
            end
        end

        identifier = [opt.username, '@', opt.computer];
        
        switch identifier

            case 'ondrej@aviendha'
                opt.resDir = '~/CNresults/res';
                opt.tmpDir = '~/CNresults/tmp';
                opt.scratchdir = '';
                
            otherwise
                error('The system %s is not configured. Add it to globalOpt to configure.', [opt.username, '@', opt.computer]);

        % *** TEMPLATES ***

        %     case 'user@WINDOWS-PC'
        %         opt.resDir = 'C:\Users\USERNAME\...\CNresults\res';
        %         opt.tmpDir = 'C:\Users\USERNAME\...\CNresults\tmp';
        %         opt.scratchdir = '';

        %     case 'user@UNIX-MACHINE'
        %         opt.resDir = '~/CNresults/res';
        %         opt.tmpDir = '~/CNresults/tmp';
        %         opt.scratchdir = '/scratch/';

        end

        %opt.MNADir = fullfile(opt.cochleadir, 'CochlearModelTime', 'time_model');

        opt.cacheDir = opt.resDir;

        end

        function extra_dirs = get_extra_dirs(obj)
            extra_dirs = obj.allowedFolders;
            for i = 1:numel(extra_dirs)
                extra_dirs{i} = fullfile(extra_dirs{i}, 'res');
            end
            
            extra_dirs = [{obj.resDir}, extra_dirs, {obj.resDir}];
        end
    end
end