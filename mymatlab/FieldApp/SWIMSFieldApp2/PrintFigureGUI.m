function SW=PrintFigureGUI(SW)
if SW.printfilenum<10
    printname=strcat(SW.printfilename,'0',num2str(SW.printfilenum));
else
    printname=strcat(SW.printfilename,num2str(SW.printfilenum));
end
figure(1)
WriteEPS(printname,SW.prpath)
SW.printfilenum=SW.printfilenum+1;
