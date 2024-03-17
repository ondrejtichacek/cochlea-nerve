function [randString] = randstr(len, characters)
arguments
    len (1,1) int32 = 1 % length of random string to generate
    characters (1,:) char = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
end

%find number of random characters to choose from
numRands = length(characters); 

%generate random string
randString = characters( ceil(rand(1,len)*numRands) );

end