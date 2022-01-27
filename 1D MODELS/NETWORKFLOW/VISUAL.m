function []=VISUAL(FNAME)

    datafile = fopen(FNAME,'r');
    [NBR,ELASTA,LONG,T_period,NCycle,DT1,Msave] = data(datafile);
    fclose(datafile);
    
    tevol = savetime(T_period,DT1,NCycle,Msave);
    
    % Visualisations of the results
    mmHg = 1334;
    ok=1;
    while (ok==1)
        dis = sprintf('Which branch do you want to explore? [0 to %i] (if 0 no visualisation)  ',NBR);
        branch_visu = input(dis);
        if ((branch_visu>0)&(branch_visu<=NBR))
             output_inlet = fopen('A_U_P_DEB_T_inlet.dat','r');
             %output_pt4 = fopen('A_U_P_DEB_T_pt4.dat','r');
             output_outlet = fopen('A_U_P_DEB_T_outlet.dat','r');

            % Entry of the branch
            for i=1:length(tevol)
                fscanf(output_inlet,'%f',5*(branch_visu-1));
                scan = fscanf(output_inlet,'%f',5);
                fscanf(output_inlet,'%f',5*(NBR-branch_visu));
                section(i) = scan(1);
                velocity(i) = scan(2);
                pressure(i) = scan(3)/mmHg;
                debit(i) = scan(4);
                tau(i) = scan(5);
            end

            title1p = strcat('\bf Section evolution at the entry of the branch N°',num2str(branch_visu));
            title1 = strcat(title1p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,section);title(title1);xlabel('\bf time');grid on;saveas(gcf, title1p(4:length(title1p)), 'fig');
            title2p = strcat('\bf Velocity evolution at the entry of the branch N°',num2str(branch_visu));
            title2 = strcat(title2p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,velocity);title(title2);xlabel('\bf time');grid on;saveas(gcf, title2p(4:length(title2p)), 'fig');
            title3p = strcat('\bf Pressure evolution at the entry of the branch N°',num2str(branch_visu));
            title3 = strcat(title3p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,pressure);title(title3);xlabel('\bf time');grid on;saveas(gcf, title3p(4:length(title3p)), 'fig');
            title4p = strcat('\bf Tau evolution at the entry of the branch N°',num2str(branch_visu));
            title4 = strcat(title4p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,tau);title(title4);xlabel('\bf time');grid on;saveas(gcf, title4p(4:length(title4p)), 'fig');
            title5p = strcat('\bf Debit evolution at the entry of the branch N°',num2str(branch_visu));
            title5 = strcat(title5p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,debit);title(title5);xlabel('\bf time');grid on;saveas(gcf, title5p(4:length(title5p)), 'fig');
    
            % Output of the branch
            for i=1:length(tevol)
                fscanf(output_outlet,'%f',5*(branch_visu-1));
                scan = fscanf(output_outlet,'%f',5);
                fscanf(output_outlet,'%f',5*(NBR-branch_visu));
                section(i) = scan(1);
                velocity(i) = scan(2);
                pressure(i) = scan(3)/mmHg;
                debit(i) = scan(4);
                tau(i) = scan(5);
            end

            title1p = strcat('\bf Section evolution at the output of the branch N°',num2str(branch_visu));
            title1 = strcat(title1p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,section);title(title1);xlabel('\bf time');grid on;saveas(gcf, title1p(4:length(title1p)), 'fig');
            title2p = strcat('\bf Velocity evolution at the output of the branch N°',num2str(branch_visu));
            title2 = strcat(title2p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,velocity);title(title2);xlabel('\bf time');grid on;saveas(gcf, title2p(4:length(title2p)), 'fig');
            title3p = strcat('\bf Pressure evolution at the output of the branch N°',num2str(branch_visu));
            title3 = strcat(title3p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,pressure);title(title3);xlabel('\bf time');grid on;saveas(gcf, title3p(4:length(title3p)), 'fig');
            title4p = strcat('\bf Tau evolution at the output of the branch N°',num2str(branch_visu));
            title4 = strcat(title4p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,tau);title(title4);xlabel('\bf time');grid on;saveas(gcf, title4p(4:length(title4p)), 'fig');
            title5p = strcat('\bf Debit evolution at the output of the branch N°',num2str(branch_visu));
            title5 = strcat(title5p,'; Longueur=',num2str(LONG(branch_visu)),'; Elastance=',num2str(ELASTA(branch_visu)));
            figure;plot(tevol,debit);title(title5);xlabel('\bf time');grid on;saveas(gcf, title5p(4:length(title5p)), 'fig');

        end
        fclose(output_inlet);
        %fclose(output_pt4);
        fclose(output_outlet);

        again = input('Do you want to explore another branch? [y]es or [n]o   ','s');
        if ((again~='y')&(again~='Y'))
            ok=0;
        end

    end

end

% Vessels data
function [NBR,ELASTA,LONG,T_period,NCycle,DT1,Msave] = data(datafile)

    fscanf(datafile,'%s',5);
    T_period = fscanf(datafile,'%f',1);
    fscanf(datafile,'%f',1);
    NCycle = fscanf(datafile,'%i',1);
    DT1 = fscanf(datafile,'%f',1);
    Msave = fscanf(datafile,'%i',1);
    fscanf(datafile,'%s',4);
    NBR = fscanf(datafile,'%i');
    fscanf(datafile,'%s',3);
    howmany=fscanf(datafile,'%i',1);
    fscanf(datafile,'%i',howmany);
    fscanf(datafile,'%s',3);
    fscanf(datafile,'%i',NBR);
    fscanf(datafile,'%s',3);
    fscanf(datafile,'%f',NBR);
    fscanf(datafile,'%s',3);
    ELASTA = fscanf(datafile,'%f',NBR);
    unitelasta = fscanf(datafile,'%f',1);
    ELASTA = unitelasta*ELASTA;
    fscanf(datafile,'%s',3);
    LONG = fscanf(datafile,'%f',NBR);
    fscanf(datafile,'%s',5);
    fscanf(datafile,'%i',NBR);
    fscanf(datafile,'%s',4);
    fscanf(datafile,'%i',NBR);

    ELASTA=ELASTA';
    LONG=LONG';

end

function [t] = savetime(T,DT,NCycle,Msave)

    e=ceil(T/DT)-fix(T/DT);
    t = [0:Msave:fix(T/DT+e)];
    for i=2:NCycle
        t = [t [(i-1)*fix(T/DT+e)+1:Msave:i*fix(T/DT+e)]];
    end
    t = t*DT;

end
