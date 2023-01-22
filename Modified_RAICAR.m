function [Outputt,N,Similarity] = RAICAR_Simple (data, repeat_num_RAICAR,ICA_Num_Comp,Selected_IC_Num_RAICAR)

for i=1:repeat_num_RAICAR
    fprintf('icaML_Run: %d \n',i);
    [SSS{1,i},AAA{1,i},UUU{1,i}]= icaML(data,ICA_Num_Comp, [] , 1);
    %     toc();
end

for i=1:repeat_num_RAICAR
    for j=1:repeat_num_RAICAR
        CRCM_cell{i,j}=double(abs(corr(SSS{i}',SSS{j}'))); % cross-realization correlation matrix
    end
end

% CRCM = cell2mat(CRCM_cell);

fi(1,1) = 1;
for fixx = 1:2    % to select Correct roww
    
    roww = fi(1,1);
    
    %% Sum mode
    for i = 1:repeat_num_RAICAR
        summ1(i) = round(sum(sum(CRCM_cell{roww,i})),3);
    end
    
    %% Max mode
    for i = 1:repeat_num_RAICAR
        [~,ir] = max(CRCM_cell{roww,i});
        [~,ic] = max(CRCM_cell{roww,i},[],2);
        ind_max{i} = [ir;ic'];
    end
    
    for i = 1:repeat_num_RAICAR
        for j = 1:repeat_num_RAICAR
            if(nnz(ind_max{i} == ind_max{j})==20)
                equl(i,j) = 1;
            end
        end
    end
    
    for i = 1:repeat_num_RAICAR
        summ2(i) = sum(find(equl(i,:)));
    end
    
    summ = summ1 + summ2;
    
    %% Union & find index
    summ_nozero = summ(logical(summ));
    N = histcounts(summ_nozero,'BinWidth',0.001); N = N(logical(N)); % number of each CELL_category
    [N, ind_N] = sort(N,'descend'); % Sort descend
    summm = union(summ_nozero,summ_nozero);  % Index_Value of each CELL_category
    for i = 1:size(summm,2)
        fi_cell{ind_N(i),1} = [find(summ == summm(i)),zeros(1,100 - size(find(summ == summm(i)),2))];
        CRCM_cell_Uni{ind_N(i)} =  CRCM_cell{roww,fi_cell{ind_N(i),1}(1)};
        SSS_Uni{ind_N(i),1} =  SSS{1,fi_cell{ind_N(i),1}(1)};
        AAA_Uni{ind_N(i),1} =  AAA{1,fi_cell{ind_N(i),1}(1)};
        UUU_Uni{ind_N(i),1} =  UUU{1,fi_cell{ind_N(i),1}(1)};
    end
    
    CRCM_Uni = cell2mat(CRCM_cell_Uni);
    Source_Infomax_Uni = cell2mat(SSS_Uni);
    fi = cell2mat(fi_cell);
    
end

Outputt.S = SSS_Uni{1,1};
Outputt.A = AAA_Uni{1,1};
Outputt.U = UUU_Uni{1,1};
Similarity = CRCM_cell_Uni;

end