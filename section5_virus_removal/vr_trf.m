function xout = vr_trf(xin,mode,qs)
% function xout = vr_trf(xin,mode,qs)
%
% Transform between abundance and carbon variables
%
% Inputs:
% . xin : variables needing to be transformed
% . mode : 
%       . mode==1: variables from abundance to carbon
%       . mode==2: variables from carbon to abundance
%       . mode==3: parameters from abundance to carbon
% . qs : carbon quotas from the parameter sets serving as conversion factors
% Outputs : 
% . xout : variables after transformation

% From abundance to carbon
if mode==1
    % transform variables from abundance to carbon
    vec=[qs.qP qs.qB qs.qVP qs.qVB qs.qG 1];
    if size(xin,1) == 6
        xout=vec'.*xin;
    elseif size(xin,2) == 6
        xout=vec.*xin;
    else
        disp('Input has wrong size')
        xout=nan;
        return
    end

% From carbon to abundance
elseif mode==2
    % transform variables from carbon to abundance
   vec=1./[qs.qP qs.qB qs.qVP qs.qVB qs.qG 1];
   if size(xin,1) == 6
       xout=vec'.*xin;
   elseif size(xin,2) == 6
       xout=vec.*xin;
   else
       disp('Input has wrong size')
       xout=nan;
       return
   end

% Parameters from abundance to carbon
elseif mode==3
    % transform parameters from abundance to carbon
    xout=struct('muPc',xin.muP, ...
                'KPc',qs.qP*xin.KP, ...
                'm1Pc',xin.m1P, ...
                'r1Pc',xin.r1P, ...
                'psPc',xin.psP/qs.qG, ...
                'muBc',xin.muB, ...
                'KBc',xin.KB, ...
                'eBc',xin.eB, ...
                'm1Bc',xin.m1B, ...
                'r1Bc',xin.r1B, ...
                'psBc',xin.psB/qs.qG, ...
                'm2Gc',xin.m2G/qs.qG,...
                'r1Gc',xin.r1G, ...
                'pic',xin.pi, ...
                'pec',xin.pe, ...
                'phPc',xin.phP/qs.qVP, ...
                'bePc',qs.qVP*xin.beP/qs.qP, ...
                'mVPc',xin.mVP, ...
                'phBc',xin.phB/qs.qVB, ...
                'beBc',qs.qVB*xin.beB/qs.qB, ...
                'mVBc',xin.mVB);
end

end