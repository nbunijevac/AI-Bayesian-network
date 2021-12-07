clc
clear all
close all

Nr = 10;
N_v = floor(linspace(1,400,100));
error_1 = zeros(Nr+1,length(N_v));
error_2 = zeros(Nr+1,length(N_v));

T = 50;
sigma_n = 1;
sigma_a = 2;
x_t = zeros(1,T);
y_t = zeros(1,T);

for i = 1:length(N_v)
    N = N_v(i);

    for j = 1:Nr
        
        x_t = zeros(1,T);
        % GENERISANJE OPSERVACIJA
        x_curr = rand()*2-1;
        for t = 1:T
            x_t(t) = 0.5*x_curr + 25*x_curr/(1+x_curr^2) + 8*cos(1.2*t) + sqrt(sigma_a)*randn();
            y_t(t) = x_t(t)^2/20 + sqrt(sigma_n)*randn();
            x_curr = x_t(t);           
        end
        
        % INICIJALIZACIJA
        x_curr_p = rand(1,N)*2-1;
        x_next_p = zeros(1,N);
        error1_curr = 0;
        error2_curr = 0;
        w_v = ones(1,N)/N;
        
        for t = 1:T
            x_curr_p = 0.5*x_curr_p + 25*x_curr_p./(1+x_curr_p.^2) + 8*cos(1.2*t) + sqrt(sigma_a)*randn(1,N);
            
            % tezinjenje
            y_curr_p = x_curr_p.^2/20; 
            w_v = w_v.*exp(-(y_curr_p-y_t(t)).^2/(2*sigma_n^2));
            w_v = w_v/sum(w_v); 
            
            x_weighted = w_v.*x_curr_p;
            mean_sample = sum(x_weighted);
            std_sample = sqrt(sum(w_v.*((x_curr_p - mean_sample).^2)));
            error1_curr = error1_curr + (abs(mean_sample - x_t(t)) > 2*std_sample);
            error2_curr = error2_curr + (mean_sample - x_t(t))^2;
            
            [w_v,index] = sort(w_v,'descend');
            x_curr_p = x_curr_p(index);
            
                        
            % reuzorkovanje
            w_v = cumsum(w_v);
            for ind = 1:N
                num = rand();
                index_curr = N - nnz((w_v-num)>0) + 1;
                x_next_p(ind) = x_curr_p(index_curr);
            end
            x_curr_p = x_next_p;
            w_v = 1/N*ones(1,N);
            
            
        end
        
        error_1(j,i) = error1_curr/T;
        error_2(j,i) = sqrt(error2_curr/T);
        
        error1_curr = 0;
        error2_curr = 0;
        
    end

    error_1(Nr+1,i) = 100*sum(error_1(1:Nr,i))/Nr;
    error_2(Nr+1,i) = (sum(error_2(1:Nr,i)))/Nr;
end


[~,ind1] = min(error_1(Nr+1,:)); [~, ind2] = min(error_2(Nr+1,:));
%N_v(ind1), N_v(ind2);

figure
plot(N_v,error_1(Nr+1,:),'m','LineWidth',1)
xlabel('broj iteracija'); ylabel('% cestica za koje vazi da |x_t - x_{pred}|>2*\sigma')
title({'Cesticni filtar','\it Greska 1'})
figure
plot(N_v,error_2(Nr+1,:),'m','LineWidth',1)
xlabel('broj iteracija'); ylabel('koren srednje-kvadratne greske estimacije')
title({'Cesticni filtar','\it Greska 2'})

%% 

N = 100;
x_sa_r = zeros(T,N);
x_bez_r = zeros(T,N);
        
x_t = zeros(1,T);
% generisanje opservacija
x_curr = rand()*2-1;
for t = 1:T
    x_t(t) = 0.5*x_curr + 25*x_curr/(1+x_curr^2) + 8*cos(1.2*t) + sqrt(sigma_a)*randn();
    y_t(t) = x_t(t)^2/20 + sqrt(sigma_n)*randn();
    x_curr = x_t(t);           
end

% SA REUZORKOVANJEM
x_curr_p = rand(1,N)*2-1;
x_next_p = zeros(1,N);
error1_curr = 0;
error2_curr = 0;
w_v = ones(1,N)/N;
x_pred_r = zeros(1,T);
std_r = zeros(1,T);

for t = 1:T
    x_curr_p = 0.5*x_curr_p + 25*x_curr_p./(1+x_curr_p.^2) + 8*cos(1.2*t) + sqrt(sigma_a)*randn(1,N);
    x_sa_r(t,:) = x_curr_p;
    % tezinjenje
    y_curr_p = x_curr_p.^2/20; %+ sqrt(sigma_n)*randn(1,N);
    w_v = w_v.*exp(-(y_curr_p-y_t(t)).^2/(2*sigma_n^2));
    w_v = w_v/sum(w_v); 

    x_weighted = w_v.*x_curr_p;
    mean_sample = sum(x_weighted);
    std_sample = sqrt(sum(w_v.*((x_curr_p - mean_sample).^2)));
    error1_curr = error1_curr + (abs(mean_sample - x_t(t)) > 2*std_sample);
    error2_curr = error2_curr + (mean_sample - x_t(t))^2;
    std_r(t) = std_sample;
    x_pred_r(t) = mean_sample;
    
    [w_v,index] = sort(w_v,'descend');
    x_curr_p = x_curr_p(index);


    % reuzorkovanje
    if (max(w_v)>0)
        w_v = cumsum(w_v);
        for ind = 1:N
            num = rand();
            index_curr = N - nnz((w_v-num)>0) + 1;
            x_next_p(ind) = x_curr_p(index_curr);% + sqrt(sigma_a)*randn(); 
        end 
        x_curr_p = x_next_p;
        w_v = 1/N*ones(1,N);
    end


end

error_1_r = 100*error1_curr/T;
error_2_r = sqrt(error2_curr/T);

% BEZ REUZORKOVANJA
x_curr_p = rand(1,N)*2-1;
x_next_p = zeros(1,N);
error1_curr = 0;
error2_curr = 0;
w_v = ones(1,N)/N;
x_pred_bez_r = zeros(1,T);
std_bez_r = zeros(1,T);

for t = 1:T
    x_curr_p = 0.5*x_curr_p + 25*x_curr_p./(1+x_curr_p.^2) + 8*cos(1.2*t) + sqrt(sigma_a)*randn(1,N);
    x_bez_r(t,:) = x_curr_p;
    % tezinjenje
    y_curr_p = x_curr_p.^2/20; %+ sqrt(sigma_n)*randn(1,N);
    w_v = w_v.*exp(-(y_curr_p-y_t(t)).^2/(2*sigma_n^2));
    w_v = w_v/sum(w_v); 

    x_weighted = w_v.*x_curr_p;
    mean_sample = sum(x_weighted);
    std_sample = sqrt(sum(w_v.*((x_curr_p - mean_sample).^2)));
    error1_curr = error1_curr + (abs(mean_sample - x_t(t)) > 2*std_sample);
    error2_curr = error2_curr + (mean_sample - x_t(t))^2;
    std_bez_r(t) = std_sample;
    x_pred_bez_r(t) = mean_sample;

end

error_1_bez_r = 100*error1_curr/T;
error_2_bez_r = sqrt(error2_curr/T);


%%

t = linspace(1,T,T);
figure
plot(t,x_t,'k','LineWidth',2); hold on;
plot(t,x_pred_r,'r','LineWidth',2); hold on;
plot(t,x_pred_r-2*std_r,'m-.','LineWidth',1); hold on;
plot(t,x_pred_r+2*std_r,'m-.','LineWidth',1); hold off;
legend('x_t','predikcije','2\sigma interval poverenja');
title('Cesticni filtar sa reuzorkovanjem')
xlabel('t')

figure
plot(t,x_t,'k','LineWidth',2); hold on;
plot(t,x_pred_bez_r,'r','LineWidth',2); hold on;
plot(t,x_pred_bez_r-2*std_bez_r,'m-.','LineWidth',1); hold on;
plot(t,x_pred_bez_r+2*std_bez_r,'m-.','LineWidth',1); hold off;
legend('x_t','predikcije','2\sigma interval poverenja');
title('Cesticni filtar bez reuzorkovanja')
xlabel('t')

