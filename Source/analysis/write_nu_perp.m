%% write_nu_para_beta.m
% eta = # sheet
%   A     B     C     D     E     F     G     H     I     J     K
%1 alpha  #          eta    #        min_gen  #        max_gen  #
%2            below         
%3  q     -     =     +     
%4  1     #     #     #     
%5  2     #     #     #     
%6  3     #     #     #     
%
% Pb sheet
%       A     B     C     D  
%1      q     1     2     3  
%2 nu_perps       below      
%3      #     #     #     #  
%
% -ADS 12*10*13

file_nu_perp = ['nu_perp_' generalize_base_name(base_name)];

%% eta sheet
eta_sheet = ['eta = ' num2str(eta)];
text_nr = cell(3,10);
text_nr{1,1} = 'alpha'; text_nr{1,4} = 'eta'; text_nr{1,7} = 'min_gen'; text_nr{1,10} = 'max_gen';
text_nr{2,3} = 'below';
text_nr{3,1} = 'q'; text_nr{3,2} = '-'; text_nr{3,3} = '='; text_nr{3,4} = '+';
xlswrite(file_nu_perp,text_nr,eta_sheet,'A1:D3');
%write alpha
xlswrite(file_nu_perp,alpha,eta_sheet,'B1:B1');
%write eta
xlswrite(file_nu_perp,eta,eta_sheet,'E1:E1');
%write gen span
xlswrite(file_nu_perp,limited(1),eta_sheet,'H1:H1');
xlswrite(file_nu_perp,limited(end),eta_sheet,'K1:K1');
%write q
xlswrite(file_nu_perp,[1:Q]',eta_sheet,['A4:A' int2str(Q+3)]);
%write best nu_paras
xlswrite(file_nu_perp,nr,eta_sheet,['B4:D' int2str(Q+3)]);

%% Pb sheet
text_Pb = cell(2,3);
text_Pb{1,1} = 'q';
text_Pb{2,1} = 'nu_perp';
text_Pb{2,3} = 'below';
%write Pb for nu_para
xlswrite(file_nu_perp,text_Pb,'Pb','A1:C2');
xlswrite(file_nu_perp,1:Q,'Pb','B1:D1');
xlswrite(file_nu_perp,nu_perps','Pb',['A3:A' int2str(Nnr+2)]);
xlswrite(file_nu_perp,Pb','Pb',['B3:D' int2str(Nnr+2)]);