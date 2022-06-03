function tM=t(M)
%tM=t(M)
%__________________________________________________________________________
% Adapt pagetranspose depending on the Matlab version

vv = ver('MATLAB');
vvNum = str2double(vv.Release(3:end-2));
vvId = vv.Release(end-1:end-1);

vvOld=false;

if vvNum<2020
    vvOld=True;
elseif vvNum==2020 && vvId < 'b'
    vvOld=True;
end

if vvOld
    tM=permute(M,[2 1 3 4 5 6 7 8 9]); % sinon
else
    tM=pagetranspose(M); % si version Matlab >= R2020B
end
end