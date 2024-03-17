function OnlineStatsTest(n)

if nargin == 0
    n = 1000;
end

%%

disp('Testing OnlineMinMax')

%%

data = rand(n,4);

[MIN, MAX] = deal([]);
for i = 1:size(data,1)
    [ MIN, MAX ] = OnlineMinMax( data(i,:), MIN, MAX );
end

compare_results(MIN, min(data));
compare_results(MAX, max(data));

%%

data = rand(1000,4,8);

[MIN, MAX] = deal([]);
for i = 1:size(data,1)
    [ MIN, MAX ] = OnlineMinMax( data(i,:,:,:), MIN, MAX );
end

compare_results(MIN, min(data));
compare_results(MAX, max(data));

%%

data = rand(1000,4,8,16);

[MIN, MAX] = deal([]);
for i = 1:size(data,3)
    [ MIN, MAX ] = OnlineMinMax( data(:,:,i,:), MIN, MAX );
end

compare_results(MIN, min(data,[],3));
compare_results(MAX, max(data,[],3));

%%

disp('Testing OnlineVariance')

%%

data = rand(n,4);

[N, MEAN, M2] = deal([]);
for i = 1:size(data,1)
    [ N, MEAN, M2, VAR, SD ] = OnlineVariance( data(i,:), N, MEAN, M2 );    
end

compare_results(MEAN, mean(data))
compare_results(VAR, var(data))
compare_results(SD, std(data))

%%

data = rand(n,4,8);

[N, MEAN, M2] = deal([]);
for i = 1:size(data,1)
    [ N, MEAN, M2, VAR, SD ] = OnlineVariance( data(i,:,:,:), N, MEAN, M2 );    
end

compare_results(MEAN, mean(data))
compare_results(VAR, var(data))
compare_results(SD, std(data))

%%

data = rand(n,4,8,16);

[N, MEAN, M2] = deal([]);
for i = 1:size(data,3)
    [ N, MEAN, M2, VAR, SD ] = OnlineVariance( data(:,:,i,:), N, MEAN, M2 );    
end

compare_results(MEAN, mean(data,3))
compare_results(VAR, var(data,0,3))
compare_results(SD, std(data,0,3))

%%
    function compare_results(r1, r2)
        
        d = abs(r1 - r2);
        if all(d(:) == 0)
            disp('[OK]   ... test passed');
        elseif all(d(:) < 100*eps)
            disp('[WARN] ... test passed up to 100 eps')
        else
            error('test failed');
        end
    end

end
