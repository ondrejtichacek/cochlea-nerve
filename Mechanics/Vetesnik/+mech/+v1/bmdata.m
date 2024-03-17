function [ GG, x, stiff, dampn, damp0, ShSp, BMinput, sh, G, M ] = bmdata(x, args)
arguments
    x = []
    args.stiff_exp (1,:) double = -3.6;
    args.stiff_0 (1,:) double = 2e+3;
    args.verflag (1,:) char = 'v1'
    args.plotflag (1,1) logical = false
end
%
% based on the function BMADATARG
%
% Generates the distributed parameters for the basilar membrane (BM) motion equation.
%
% X = coordinates of a set of points along the BM in BETA units. BETA is the BM
% length, therefore X range from 0 to 1. The distributed parameters are computed
% for any set of points. X may be an equally spaced set of  values covering the
% interval 0-1, or a set of values non-uniformly spaced along the BM (see function VARGRID.M)

if isempty(x)
    x = (1:400)/400;
    args.verflag = 'v1';
    args.plotflag = 1;
end

x = x(:);
N = length(x);
dx = [x(1); diff(x)];

% The routine GREENF.M generates the Green's function matrix G, the stapes-BM
% coupling coefficient vector Gs, and the shearing viscosity matrix Sh for the
% human cochlea.
[Gs, G, Sh, sh] = mech.greenf(x, args.plotflag);

ShSp = sparse(Sh);

% MASS.M computes the organ of Corti local mass as a function of position.
M = diag(mech.mass(x, args.plotflag));

% Ginv = inv(G+M);
% Ginv is  the inverse of the BM integro-differential motion equation kernel G+M.
% The unit length dimension of the model is BETA = BM length; consequently the
% units of the quantities are Ginv [BETA/Kg]; Gs [Kg/BETA]; Sh [Kg/(BETA*sec)]; x [BETA]

GG = G + M;

% BMinput = Ginv*Gs'; % Adimensional quantity.
BMinput = Gs';

% BMinput represents the ratio between the BM displacement and stapes 
% displacement as the direct effect of stapes motion on the BM and the 
% BM-to-BM reaction at zero BM stiffness and viscosity (BMinput decays 
% rapidly from base to apex starting from about 250)

%% STIFFNESS
% disp('    Computing stiffness ...')

switch args.verflag
    case 'v1'
        % original (probably by Nobili)

        expk = -3.6*log(10); 	% Stiffness exponent  (Base E)
        ko = 2e+3;		  % Kg/(BETA*sec^2) , BETA = BM length taken as lenght unit
        
        % Averaging stiffness over dx
        k = ko*exp(expk*x);
        k = k(:);
        stiff = [k(1)-ko; diff(k)]./(expk*dx);
    
    case 'v4'

        stiff = 0.45e-6*(2*pi*165.4*(10.^(2.1*(1-x)) - 0.88)).^2;

    case 'v2'
        alpha = args.stiff_exp;
        
        expk = alpha*log(10); 	% Stiffness exponent  (Base E)
        
        ko = args.stiff_0;		  % Kg/(BETA*sec^2) , BETA = BM length taken as lenght unit
        
        k = ko*exp(expk*x);
        k = k(:);
        stiff = [k(1)-ko; diff(k)]./(expk*dx);
    
    case 'v3'
        alpha = args.stiff_exp;
        
        expk = alpha*log(10); 	% Stiffness exponent  (Base E)
        
        ko = args.stiff_0;		  % Kg/(BETA*sec^2) , BETA = BM length taken as lenght unit
        
        stiff = zeros(size(x));
        for i = 1:numel(alpha)
            k = ko(i)*exp(expk(i)*x);
            k = k(:);
            st = [k(1)-ko(i); diff(k)]./(expk(i)*dx);

            stiff = stiff + st;
        end

        assert(all(stiff > 0));
        
    case 'v1.1'
        % modified by OT

        % cf = 165.4 * (10.^(2.1*(1-x)) - 0.88); % Greenwood
        % cf = 991.395 * (10.^(1.283*(1-x)) - 0.98); % Li et al. 2021 10.1038/s41598-021-83225-w

        alpha = 2.1; % Greenwood
        % alpha = 1.283; % Li et al. 2021

        % offset = 0.0;
        offset = 1.1;
        % offset = 1.5;

        expk = -(alpha + offset)*log(10); 	% Stiffness exponent  (Base E)
        ko = 2.0e+3;		  % Kg/(BETA*sec^2) , BETA = BM length taken as lenght unit
        
        k0 = 0.88/165.4 * ko / (10^2.1/165.4 - 0.88/165.4);
        k0 = k0*0.2;
        k0 = k0*0;

        k = ko*exp(expk*x) - k0;
        k = k(:);
        stiff = k;

%         figure
%         plot(x, k)

%         figure
%         hold on
%         cf = 165.4 * (10.^(2.1*(1-x)) - 0.88);
%         plot(x, cf/cf(1))
%         % cf = 991.395 * (10.^(1.283*(1-x)) - 0.98);
%         % plot(x, cf/cf(1))
%         plot(x, k/k(1));

    otherwise
        error('something is wrong')
end




if args.plotflag
    figure
    plot(x, stiff)
    hold on
    plot(x, k)
    ylabel('stiffness')
end

%% DAMPING
% damping constant is imagined to depend mainly the
% intrinsic viscosity of cochlear motor (OHC+DEITERS CELLS)
% h = eta*L/H ; eta = viscosity coefficient
% [3.35e-5 kg/(BETA*sec) for water, see GREENF.M]
% L/H = shearing ratio between length L of viscous segment and
% thickness H of shearing medium (about 10)
%
%		Author: Renato Nobili - Padova University, Italy (October 2000)

% disp('    Computing damping ...')

ho = 10e-3; 	% Kg/(BETA sec) [about 30 times water viscosity]
exph = -log(4);	% decreases 4 times
dampn = ho*exp(exph*x);

hbo = 1.1*ho;
heo = 1.3*dampn(N);
expho = log(heo/hbo);
damp0 = hbo*exp(expho*x);

% Vencovsky 2019
% damp0 = ho*sqrt(stiff/stiff(1));
% dampn = damp0;

% if plotflag
%     figure
%     hold on
%     plot(x, damp0 - dampn)
% end

%% Extra damping
% increase apical damping to prevent long-lasting oscillations

% version by Vetesnik
Ce = 2*dampn(N);
dampc_v1 = Ce*exp((x-1)/0.075);
dampn_v1 = dampn + dampc_v1;

% version by OT
Ce = 4*dampn(N);
dampc_v2 = Ce*exp(40*(x-1));
dampn_v2 = dampn + dampc_v2;

switch args.verflag
    case 'v1'
        dampn = dampn_v1;
    case {'v2', 'v3', 'v4'}
        dampn = dampn_v2;
    otherwise
        error('version must be one of v1, v2, not %s', verflag);
end
   
dampn = dampn(:);
damp0 = damp0(:);

%%

if args.plotflag
    figure
    hold on
    plot(x, damp0)
    plot(x, dampn_v2)
    plot(x, dampc_v2)
    plot(x, dampn_v1)
    plot(x, dampc_v1)
    ylabel('damping')
    legend({'damp0', 'dampn_v2', 'dampc_v2', 'dampn_v1', 'dampc_v1'})
end
