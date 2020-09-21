function [endm critere_opt] = nfindr(y_proj,R)

% size(y_red) = Nb_pix x Nb_bandes
[Nb_pix Nb_bandes] = size(y_proj);

Nb_endm = Nb_bandes+1;

% quick hull if size allows it
if Nb_endm<7
    tmp = convhulln(y_proj);
    ind_envlp = unique(tmp(:));
else
    ind_envlp = 1:Nb_pix;
end
envlp = y_proj(ind_envlp,:);
Nb_pix_envlp = length(ind_envlp);

% INITIALIZATION
% choice of Nb_endm pixels
ind_perm = randperm(Nb_pix_envlp);
combi = sort(ind_perm(1:Nb_endm));

% cacul du critère
% keyboard
candidat_opt = y_proj(combi,1:Nb_bandes)';
critere_opt = -abs(det([candidat_opt;ones(1,Nb_endm)]));

% Algorithm
h = waitbar(0,['N-FINDR...'],'CreateCancelBtn','closereq', 'name', 'N-FINDR');
for ind_endm=1:Nb_endm
    for ind_pix=1:Nb_pix_envlp
        waitbar((Nb_pix_envlp*(ind_endm-1)+ind_pix)/(Nb_pix_envlp*Nb_endm),h)

        combi_cand = combi;
        combi_cand(ind_endm) = ind_pix;
        if length(unique(combi_cand)) == Nb_endm
            candidat_cand = envlp(combi_cand,1:Nb_bandes)';
            critere_cand = -abs(det([candidat_cand;ones(1,Nb_endm)]));

            %test
            if critere_cand<critere_opt
                critere_opt = critere_cand;
                candidat_opt = candidat_cand;
                combi = combi_cand;
            end
        end
    end
end
delete(h)
endm = candidat_opt;


