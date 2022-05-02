function write_to_excel(file, value, sheet, range)
% almost the same as xlswrite, except for an additional waiting time

xlswrite(file, value, sheet, range);
% waiting a bit because otherwise the excel file might still be open
% and you might get errors in matlab because of file still in use..
pause(1);

end