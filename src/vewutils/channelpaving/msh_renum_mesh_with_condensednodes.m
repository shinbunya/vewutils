function [ma, map_renum] = msh_renum_mesh_with_condensednodes(ma)
[ma,map_renum] = ma.renum();

% Renum condensed_nodes values
for k=1:ma.f13.nAttr
    if strcmp(ma.f13.userval.Atr(k).AttrName,'condensed_nodes')
        val = ma.f13.userval.Atr(k).Val(2:end,:);
        idx = find(val~=0);
        val(idx) = map_renum(val(idx));
        ma.f13.userval.Atr(k).Val(2:end,:) = val;
    end
end
end