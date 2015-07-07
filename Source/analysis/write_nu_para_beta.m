%% write_nu_para_beta.m
% eta = # sheet
%   A     B     C     D     E     F     G     H
%1 alpha  #          eta    #
%2            below                   above
%3  q     -     =     +           -     =     +
%4  1     #     #     #           #     #     #
%5  2     #     #     #           #     #     #
%6  3     #     #     #           #     #     #
%
% Pb sheet
%       A     B     C     D     E     F     G     H
%1      q     1     2     3           1     2     3
%2 nu_paras       below                   above
%3      #     #     #     #           #     #     #
%
% -ADS 12*10*13

file_nu_para_beta = ['nu_para_beta_' generalize_base_name(base_name)];

%% eta sheet
eta_sheet = ['eta = ' num2str(eta)];
text_np = cell(3,10);
text_np{1,1} = 'alpha'; text_np{1,4} = 'eta'; 
  text_np{1,7} = {'below min_gen'}; text_np{1,10} = {'below max_gen'};
text_np{2,7} = {'above min_gen'}; text_np{2,10} = {'above max_gen'};
text_np{2,3} = 'below'; text_np{2,7} = 'above';
text_np{3,1} = 'q'; text_np{3,2} = '-'; text_np{3,3} = '='; text_np{3,4} = '+';
  text_np{3,6} = '-'; text_np{3,7} = '='; text_np{3,8} = '+';
xlswrite(file_nu_para_beta,text_np,eta_sheet,'A1:H3');
%write alpha
xlswrite(file_nu_para_beta,alpha,eta_sheet,'B1:B1');
%write eta
xlswrite(file_nu_para_beta,eta,eta_sheet,'E1:E1');
%write gen span
xlswrite(file_nu_para_beta,[limitedb(1);limiteda(end)],eta_sheet,'H1:H2');
xlswrite(file_nu_para_beta,[limitedb(end);limiteda(end)],eta_sheet,'K1:K2');
%write q
q_offset = int2str(Q+3);
xlswrite(file_nu_para_beta,[1:Q]',eta_sheet,['A4:A' q_offset]);
%write best nu_paras
xlswrite(file_nu_para_beta,best_npb,eta_sheet,['B4:D' q_offset]);
xlswrite(file_nu_para_beta,best_npa,eta_sheet,['F4:H' q_offset]);

%% Pb sheet
Pb_sheet = 'Pb';
text_Pb = cell(2,8);
text_Pb{1,1} = 'q';
text_Pb{2,1} = 'nu_para';
text_Pb{2,3} = 'below';
text_Pb{2,7} = 'above';
%write Pb for nu_para
xlswrite(file_nu_para_beta,text_Pb,'Pb','A1:H2');
xlswrite(file_nu_para_beta,1:Q,Pb_sheet,'B1:D1');
xlswrite(file_nu_para_beta,1:Q,Pb_sheet,'F1:H1');
xlswrite(file_nu_para_beta,nu_paras','Pb',['A3:A' int2str(Nnp+2)]);
xlswrite(file_nu_para_beta,Pbb','Pb',['B3:D' int2str(Nnp+2)]);
xlswrite(file_nu_para_beta,Pba','Pb',['F3:H' int2str(Nnp+1)]);