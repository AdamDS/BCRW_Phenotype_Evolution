%% print_fails
% [] = print_fails(outs)
% prints each line of output results
%
try,  
  print_report_fail(pout,do_DNE,do_only_DNE);
  print_report_fail(txout,do_DNE,do_only_DNE);
  print_report_fail(tyout,do_DNE,do_only_DNE);
  print_report_fail(tcsout,do_DNE,do_only_DNE);
  print_report_fail(sdout,do_DNE,do_only_DNE);
  print_report_fail(kout,do_DNE,do_only_DNE);
  print_report_fail(rout,do_DNE,do_only_DNE);
  print_report_fail(parout,do_DNE,do_only_DNE);
  print_report_fail(ncout,do_DNE,do_only_DNE);
  print_report_fail(oncout,do_DNE,do_only_DNE);
  print_report_fail(tcout,do_DNE,do_only_DNE);
  print_report_fail(cxout,do_DNE,do_only_DNE);
  print_report_fail(cyout,do_DNE,do_only_DNE);
  print_report_fail(grout,do_DNE,do_only_DNE);
  print_report_fail(divout,do_DNE,do_only_DNE);
  print_report_fail(diaout,do_DNE,do_only_DNE);
  print_report_fail(plout,do_DNE,do_only_DNE);
  print_report_fail(ndout,do_DNE,do_only_DNE);
  print_report_fail(ndcout,do_DNE,do_only_DNE);
  print_report_fail(dcout,do_DNE,do_only_DNE);
  print_report_fail(ncpout,do_DNE,do_only_DNE);
  print_report_fail(cpout,do_DNE,do_only_DNE);
  print_report_fail(ncfout,do_DNE,do_only_DNE);
  print_report_fail(cfout,do_DNE,do_only_DNE);
catch,  
  fprintf('No output strings for fail \n');
end