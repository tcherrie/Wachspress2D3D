function couleurs=legend2color(legende,inv)
if nargin<=1 || isempty(inv)
    inv=0;
end

dict_colors.fer=[0 0 0]-0.5;
dict_colors.iron=[0 0 0]-0.5;
dict_colors.fesi=[0 0 0]-0.5;
dict_colors.fesi1=[0 0 0]-0.5;
dict_colors.fesi2=[0 0 0]-0.5;
dict_colors.feco=[0 0 0]-0.5;
dict_colors.feni=[0 0 0]-0.5;
dict_colors.smc=[0 0 0]-0.5;

dict_colors.air=[1 1 1]-0.5;

dict_colors.ap=[1 0 0]-0.5;
dict_colors.am=[1 0 1]-0.5;
dict_colors.bp=[0 1 0]-0.5;
dict_colors.bm=[1 1 0]-0.5;
dict_colors.cp=[0 0 1]-0.5;
dict_colors.cm=[0 1 1]-0.5;


dict_colors.("v1")=[0 0 0]-0.5;
dict_colors.("v2")=[1 0 0]-0.5;
dict_colors.("v3")=[1 0 1]-0.5;
dict_colors.("v4")=[0 1 0]-0.5;
dict_colors.("v5")=[1 1 0]-0.5;
dict_colors.("v6")=[0 0 1]-0.5;
dict_colors.("v7")=[0 1 1]-0.5;
dict_colors.("v8")=[1 1 1]-0.5;
dict_colors.("v9")=[0.8 0.2 0.2]-0.5;
dict_colors.("v10")=[0.8 0.2 0.8]-0.5;
dict_colors.("v11")=[0.2 0.8 0.2]-0.5;
dict_colors.("v12")=[0.8 0.8 0.2]-0.5;
dict_colors.("v13")=[0.2 0.2 0.8]-0.5;
dict_colors.("v14")=[0.2 0.8 0.8]-0.5;
dict_colors.("v15")=[0.6 0.1 0.1]-0.5;
dict_colors.("v16")=[0.1 0.1 0.6]-0.5;
dict_colors.("v17")=[0.1 0.6 0.1]-0.5;
dict_colors.("v18")=[0.6 0.6 0.1]-0.5;
dict_colors.("v19")=[0.1 0.1 0.6]-0.5;
dict_colors.("v20")=[0.1 0.1 0.6]-0.5;

dict_colors.dp=[0.8 0.2 0.2]-0.5;
dict_colors.dm=[0.8 0.2 0.8]-0.5;
dict_colors.ep=[0.2 0.8 0.2]-0.5;
dict_colors.em=[0.8 0.8 0.2]-0.5;
dict_colors.fp=[0.2 0.2 0.8]-0.5;
dict_colors.fm=[0.2 0.8 0.8]-0.5;
dict_colors.gp=[0.6 0.1 0.1]-0.5;
dict_colors.gm=[0.1 0.1 0.6]-0.5;
dict_colors.hp=[0.1 0.6 0.1]-0.5;
dict_colors.hm=[0.6 0.6 0.1]-0.5;
dict_colors.ip=[0.1 0.1 0.6]-0.5;
dict_colors.im=[0.1 0.1 0.6]-0.5;

dict_colors.mag1=[1 0 0]-0.5;
dict_colors.mag2=[1 0 1]-0.5;
dict_colors.mag3=[0 1 0]-0.5;
dict_colors.mag4=[1 1 0]-0.5;
dict_colors.mag5=[0 0 1]-0.5;
dict_colors.mag6=[0 1 1]-0.5;
dict_colors.mag7=[1 0 0]-0.5;
dict_colors.mag8=[1 0 1]-0.5;
dict_colors.mag9=[0 1 0]-0.5;
dict_colors.mag10=[1 1 0]-0.5;
dict_colors.mag11=[0 0 1]-0.5;
dict_colors.mag12=[0 1 1]-0.5;

dict_colors.magnet1=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet2=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet3=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet4=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet5=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet6=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet7=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet8=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet9=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet10=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet11=[0.3 0.3  0.3 ]-0.5;
dict_colors.magnet12=[0.3 0.3  0.3 ]-0.5;


dict_colors.mag_h=[1 0 0]-0.5;
dict_colors.mag_b=[1 0 1]-0.5;
dict_colors.mag_l=[0 1 0]-0.5;
dict_colors.mag_r=[1 1 0]-0.5;

dict_colors.mag_n=[1 0 0]-0.5;
dict_colors.mag_s=[1 0 1]-0.5;
dict_colors.mag_e=[0 1 0]-0.5;
dict_colors.mag_w=[1 1 0]-0.5;

dict_colors.mag_up=[1 0 0]-0.5;
dict_colors.mag_down=[1 0 1]-0.5;
dict_colors.mag_right=[0 1 0]-0.5;
dict_colors.mag_left=[1 1 0]-0.5;


dict_colors.fer1=[0 0 0]-0.5;
dict_colors.fer2=[0 0 1]-0.5;
dict_colors.fer3=[0 1 0]-0.5;
dict_colors.fer4=[0 1 1]-0.5;
dict_colors.fer5=[1 0 0]-0.5;
dict_colors.fer6=[1 0 1]-0.5;
dict_colors.fer7=[1 1 0]-0.5;

dict_colors.conducteur=[1 0.5 0]-0.5;
dict_colors.aimant=[1 0 0]-0.5;
dict_colors.conducteurs=[1 0.5 0]-0.5;
dict_colors.aimants=[1 0 0]-0.5;
dict_colors.mag=[1 0 0]-0.5;
dict_colors.fers=[1 0 0]-0.5;

dict_colors.inductor=[0.8 0 0.2]-0.5;
dict_colors.inducteur=[0.8 0 0.2]-0.5;
dict_colors.indp=[0.8 0 0.2]-0.5;
dict_colors.indm=[0.2 0 0.8]-0.5;


if nargin==0 || isempty(legende)
    couleurs=convertCharsToStrings(fieldnames(dict_colors));
else
legende(lower(legende)=="a+")="ap";
legende(lower(legende)=="a-")="am";
legende(lower(legende)=="b+")="bp";
legende(lower(legende)=="b-")="bm";
legende(lower(legende)=="c+")="cp";
legende(lower(legende)=="c-")="cm";
legende(lower(legende)=="d+")="dp";
legende(lower(legende)=="d-")="dm";
legende(lower(legende)=="e+")="ep";
legende(lower(legende)=="e-")="em";
legende(lower(legende)=="f+")="fp";
legende(lower(legende)=="f-")="fm";
legende(lower(legende)=="g+")="gp";
legende(lower(legende)=="g-")="gm";
legende(lower(legende)=="h+")="hp";
legende(lower(legende)=="h-")="hm";
legende(lower(legende)=="i+")="ip";
legende(lower(legende)=="i-")="im";

legende(lower(legende)=="ind-")="indm";
legende(lower(legende)=="ind+")="indp";

legend2=legende;
if inv
    legend2(lower(legende)=="ap")="am";
    legend2(lower(legende)=="am")="ap";
    legend2(lower(legende)=="bm")="bp";
    legend2(lower(legende)=="bp")="bm";
    legend2(lower(legende)=="cm")="cp";
    legend2(lower(legende)=="cp")="cm";
    legend2(lower(legende)=="dp")="dm";
    legend2(lower(legende)=="dm")="dp";
    legend2(lower(legende)=="em")="ep";
    legend2(lower(legende)=="ep")="em";
    legend2(lower(legende)=="fm")="fp";
    legend2(lower(legende)=="fp")="fm";
    legend2(lower(legende)=="gp")="gm";
    legend2(lower(legende)=="gm")="gp";
    legend2(lower(legende)=="hm")="hp";
    legend2(lower(legende)=="hp")="hm";
    legend2(lower(legende)=="im")="ip";
    legend2(lower(legende)=="ip")="im";
    
    legend2(lower(legende)=="mag_up")="mag_down";
    legend2(lower(legende)=="mag_down")="mag_up";
    legend2(lower(legende)=="mag_n")="mag_s";
    legend2(lower(legende)=="mag_s")="mag_n";
    legend2(lower(legende)=="mag_e")="mag_w";
    legend2(lower(legende)=="mag_w")="mag_e";
    legend2(lower(legende)=="mag_right")="mag_left";
    legend2(lower(legende)=="mag_left")="mag_right";
    legend2(lower(legende)=="mag_h")="mag_b";
    legend2(lower(legende)=="mag_b")="mag_h";
    legend2(lower(legende)=="mag_l")="mag_r";
    legend2(lower(legende)=="mag_r")="mag_l";
    
    legend2(lower(legende)=="indp")="indm";
    legend2(lower(legende)=="indm")="indp";
end


for i=1:length(legend2)
    try
    couleurs(i,:)=dict_colors.(lower(legend2(i)));
    catch
        couleurs(i,:)=[0,0,0];
    end
end

end


end