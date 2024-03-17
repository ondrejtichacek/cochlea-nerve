function [ username, computername, system_spec ] = userAndComputerName(  )
%USERANDCOMPUTERNAME Returns user and computer name independently on system

if isunix
    
    [~,username] = system('id -u -n'); % username    
    [~,computername] = system('uname -n'); % computername
    
    username = strtrim(username);
    computername = strtrim(computername);
    system_spec = 'unix';
    
elseif ispc

%     [s,c] = dos('ECHO %COMPUTERNAME%');
%     computername = strtrim(c);

    username = getenv('UserName');
    computername = getenv('ComputerName');
    system_spec = 'pc';
    
elseif ismac    
    
    system_spec = 'mac';
    username = getenv('USERNAME');
    computername = getenv('HOSTNAME');

end

end

