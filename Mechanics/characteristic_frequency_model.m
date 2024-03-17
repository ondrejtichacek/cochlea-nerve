function [ cx, cf ] = characteristic_frequency_model(version, args)
%CHARACTERISTIC_FREQUENCY_MODEL
arguments
    version (1,:) char = 'Greenwood';
    args.plotflag (1,1) logical = false
end

switch version
    case 'Vetesnik v1'
        base = 2.4531;
        alpha = -6.36;
        Fo = 21100;
        
        cf = @(x) Fo * base.^(alpha*(x));
        cx = @(c) log(c/Fo)/log(base^alpha);

    case 'Vetesnik v2'
        base = 2.2022;
        alpha = -7.3094;
        Fo = 20899;
        
        cf = @(x) Fo * base.^(alpha*(x));
        cx = @(c) log(c/Fo)/log(base^alpha);

    case {'Greenwood', 'Greenwood 1990'}
        % Greenwood 1990
        % DOI:10.1121/1.399052
        A = 165.4;
        alpha = 2.1;
        K = 0.88;

        cf = @(x) A * (10.^(alpha*(1-x)) - K);
        cx = @(f) 1 - log10(f/A + K) / alpha;

    case {'Li', 'Li 2021'}
        % Li et al. 2021 
        % DOI:10.1038/s41598-021-83225-w
        A = 991.395;
        alpha = 1.283;
        K = 0.98;

        cf = @(x) A * (10.^(alpha*(1-x)) - K);
        cx = @(f) 1 - log10(f/A + K) / alpha;

    otherwise
        error('Something is wrong')
end

if args.plotflag == true
    plot_all();
end

end

function plot_all()
x = linspace(0,1,200);

versions = {'Vetesnik v1', 'Vetesnik v2', 'Greenwood', 'Li'};

figure
hold on
for i = 1:numel(versions)
    
    [~, cf] = characteristic_frequency_model(versions{i}, 'plotflag', false);
    plot(x, cf(x), 'DisplayName', versions{i});
end
legend()

end

