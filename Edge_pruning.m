function [thresholded,numl_FC,numl_FC_sel] = Matrix_thresholding_justFC(Input,thresh_sel)
%%Initial thresholds
thresh_FC = 0.1:0.1:4; 
% thresh_SC = 0.1:0.1:4;
num_IC=10;
%% U * A * S
S = Input.S;
A = Input.A;
U = Input.U;
UA = U(1:num_IC,1:num_IC) * A;
Y = UA*S;

size_S=size(S,2);

Y_FC = Y(:,1:size_S);
% Y_SC = Y(:,8256+1:end);
S_FC = S(:,1:size_S);
% S_SC = S(:,8256+1:end);

%% Main of Code  ->  FC
NumofConn = size(Y_FC,2);
for jj = 1:NumofConn
    temp = Y_FC(2:num_IC,jj);
    temp = temp .^2;
    YY_FC(jj) = sqrt(sum(temp));    
end

for ind = 1:numel(thresh_FC)
    for row = 1:num_IC
        Coef_FC = max(abs(UA(2:num_IC,row)))./ YY_FC ;
        Mat_FC(row,:) = S_FC(row,:) .* Coef_FC;
        check_FC(row,:) = abs(Mat_FC(row,:)) > thresh_FC(ind);
    end
    
    for row = 1:10
        numl_FC(ind,row) = round(numel(find(check_FC(row,:) == 1))./NumofConn.*100);
    end
    
    if ind == thresh_sel
        check_FC_sel = check_FC;
        numl_FC_sel = numl_FC(ind,:);
    end
end

%% Thresholded
thresholded_FC = S_FC;
% thresholded_SC = S_SC;
thresholded_FC(find(~check_FC_sel)) = 0;
% thresholded_SC(find(~check_SC_sel)) = 0;
thresholded = [thresholded_FC];


