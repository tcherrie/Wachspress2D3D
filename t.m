function tM=t(M)
%tM=t(M)
%__________________________________________________________________________
% Raccourci pour pagetranspose qui calcule la transpos�e.
% Si la version de Matlab est inf�rieure � 2020B, remplacer le corps de
% cette fonction par " tM=permute(M,[2 1 3 4 5 6 7 8 9]); " 
%__________________________________________________________________________


tM=pagetranspose(M); % si version Matlab >= R2020B
%tM=permute(M,[2 1 3 4 5 6 7 8 9]); % sinon

end