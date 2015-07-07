%% print_passes
% [] = print_passes(outs)
% prints each line of output results
%

try,  
  print_report_pass(pout);
  print_report_pass(txout);
  print_report_pass(tyout);
  print_report_pass(tcsout);
  print_report_pass(sdout);
  print_report_pass(kout);
  print_report_pass(rout);
  print_report_pass(parout);
  print_report_pass(ncout);
  print_report_pass(oncout);
  print_report_pass(tcout);
  print_report_pass(cxout);
  print_report_pass(cyout);
  print_report_pass(grout);
  print_report_pass(divout);
  print_report_pass(diaout);
  print_report_pass(plout);
  print_report_pass(ndout);
  print_report_pass(ndcout);
  print_report_pass(dcout);
  print_report_pass(ncpout);
  print_report_pass(cpout);
  print_report_pass(ncfout);
  print_report_pass(cfout);
catch,  
  fprintf('No output strings for pass');
end