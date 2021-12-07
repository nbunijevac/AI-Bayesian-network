clc
clear all 
close all

% Pi_J =  [i-|j-; i-|j+]

Pb_A = [0.7; 0.2];
Pc_A = [0.5; 0.7];
Pd_A = [0.3; 0.8];
Pg_A = [0.4; 0.6];
Pf_D = [0.8; 0.3];
Pe_AB = [0.6 0.5; 0.3 0.8];
Pa = 0.3;

N = 50000;
Nr = 100;

%% 1 zadatak b) uzorkovanje sa odbacivanjem 

% P(e-|f+) ~ #(a*,b*,d*,e-,f+)/#(a*,b*,d*,e*,f+)
% x = 1/2 <=> x = x-/x+

p_emfp_b = zeros(1,Nr);
e = 1; f = 2;
    
for k = 1:Nr
    % ukupni (imenilac iz gornjeg razlomka)
    cnt_fp = 0;
    % odgovarajuci (brojilac iz gornjeg razlomka)
    cnt_emfp = 0;

    for i = 1:N
        a = 1; b = 1; d = 1;
        % uzorkovanje A
        if (rand()>Pa)
            a = 2;
        end
        % uzorkovanje B
        if (rand()>Pb_A(a))
            b = 2;
        end
        % uzorkovanje d
        if (rand()>Pd_A(a))
            d = 2;
        end
        % uzorkovanje f
        if (rand()>Pf_D(d))
            cnt_fp = cnt_fp + 1; 
            % samo ako je f = f+, 
            % uzokrujemo e
            if (rand()<Pe_AB(a,b))
                % ako je e = e-
                cnt_emfp = cnt_emfp + 1;
            end
        end

    end

    p_emfp_b(k) = cnt_emfp/cnt_fp;

end

%%  1 zadatak c) metoda ponderisanja verodostojnoscu

% cnt = [#(*,*,d-,e-,f+) #(*,*,d+,e-,f+);   --> e-
%        #(*,*,d-,e+,f+) #(*,*,d+,e+,f+)]   --> e+

p_emfp_c = zeros(1,Nr);
cnt = zeros(2,2);

for k = 1:Nr
    for i = 1:N
        a = 1; b = 1; d = 1; e = 1;
        % uzorkovanje a
        if (rand()>Pa)
            a = 2;
        end
        % uzorkovanje b
        if (rand()>Pb_A(a))
            b = 2;
        end
        % uzorkovanje d
        if (rand()>Pd_A(a))
            d = 2;
        end  
        % uzorkovanje e
        if (rand()>Pe_AB(a,b))
            e = 2;
        end
        % brojanje
        cnt(e,d) = cnt(e,d) + 1;
    end

    p_emfp = cnt(1,1)*(1-Pf_D(1,1)) + cnt(1,2)*(1-Pf_D(2,1));
    p_epfp = cnt(2,1)*(1-Pf_D(1,1)) + cnt(2,2)*(1-Pf_D(2,1));

    alfa = 1/(p_emfp + p_epfp);
    p_emfp_c(k) = alfa*p_emfp;

end

%%  1 zadatak d) Gibbsovo uzorkovanje

% Pi_J =  [i-|j- i+|j-   -> j = j-
%          i-|j+ i+|j+]  -> j = j+


Pb_A = [0.7 0.3; 0.2 0.8];
Pc_A = [0.5 0.5; 0.7 0.3];
Pd_A = [0.3 0.7; 0.8 0.2];
Pg_A = [0.4 0.6; 0.6 0.4];
Pf_D = [0.8 0.2; 0.3 0.7];
Pe_AB = zeros(2,2,2);
Pe_AB(:,:,1) = [0.6 0.5; 0.3 0.8];
Pe_AB(:,:,2) = 1-Pe_AB(:,:,1);
Pa = 0.3;

f = 2;
line = randi(2,1,4);
a = line(1);
b = line(2);
d = line(3);
e = line(4);

p_emfp_d = zeros(1,Nr);


for k = 1:Nr
    cnt_em = 0;
    cnt_f = 0;
    for i = 1:N
        a = 1;
        p_am = Pa*Pb_A(1,b)*Pd_A(1,d)*Pe_AB(1,b,e);      % a = a-;
        p_ap = (1-Pa)*Pb_A(2,b)*Pd_A(2,d)*Pe_AB(2,b,e);  % a = a+;
        p_a = p_am/(p_am+p_ap);
        if (rand()>p_a)
            a = 2;
        end
        b = 1;
        p_bm = Pe_AB(a,1,e)*Pb_A(a,1);   % b = b-
        p_bp = Pe_AB(a,2,e)*Pb_A(a,2);   % b = b+
        p_b = p_bm/(p_bm+p_bp);
        if (rand()>p_b)
            b = 2;
        end
        d = 1;
        p_dm = Pf_D(1,f)*Pd_A(a,1);   % d = d-
        p_dp = Pf_D(2,f)*Pd_A(a,2);   % d = d+
        p_d = p_dm/(p_dm+p_dp);
        if (rand()>p_d)
            d = 2;
        end
        %if (rand()>Pf_D(d,f))
        %    cnt_f = cnt_f + 1;
            e = 2;
            if (rand()<Pe_AB(a,b,1))
                e = 1;
                cnt_em = cnt_em + 1;
            end
        %end
    end

    p_emfp_d(k) = cnt_em/N;

end

[p_emfp_b,p_emfp_c,p_emfp_d];

%% histogrami 

correct_val = 0.6428;
mean_b = mean(p_emfp_b);
std_b = std(p_emfp_b');
figure
Hb = histogram(p_emfp_b); hold on; 
b_lim = max(Hb.Values) + 5;
plot([correct_val correct_val],[0,b_lim],'r'); hold on;
plot([mean_b mean_b],[0,b_lim],'r--'); hold on;
plot([mean_b + std_b, mean_b + std_b],[0,b_lim],'k--',[mean_b - std_b, mean_b - std_b],[0,b_lim],'k--');
title('Uzorkovanje sa odbacivanjem')
ylim([0,b_lim]); legend('histogram','tacna vrednost','mean','\sigma')

mean_c = mean(p_emfp_c);
std_c = std(p_emfp_c');
figure
Hc = histogram(p_emfp_c); hold all;
c_lim = max(Hc.Values) + 5;
plot([correct_val correct_val],[0,c_lim],'r'); hold on;
plot([mean_c mean_c],[0,c_lim],'r--'); hold on;
plot([mean_c + std_c, mean_c + std_c],[0,c_lim],'k--',[mean_c - std_c, mean_c - std_c],[0,c_lim],'k--');
title('Metoda ponderisanja verodostojnoscu')
ylim([0,c_lim]); legend('histogram','tacna vrednost','mean','\sigma')

mean_d = mean(p_emfp_d);
std_d = std(p_emfp_d');
figure
Hd = histogram(p_emfp_d); hold on;
d_lim = max(Hd.Values) + 5;
plot([correct_val correct_val],[0,d_lim],'r'); hold on;
plot([mean_d mean_d],[0,d_lim],'r--'); hold on;
plot([mean_d + std_d, mean_d + std_d],[0,d_lim],'k--',[mean_d - std_d, mean_d - std_d],[0,d_lim],'k--');
title('Gibbsovo uzorkovanje')
ylim([0,d_lim]); legend('histogram','tacna vrednost','mean','\sigma')



