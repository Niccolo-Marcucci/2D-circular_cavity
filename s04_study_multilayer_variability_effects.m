% Copyright 2020 Niccolò Marcucci <niccolo.marcucci@polito.it>
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('Lumerical-Objects/multilayer_design/functions');

design_name = "TM_Descrovi";
design_file = strcat("Lumerical-Objects/multilayer_design/designs/design_",...
    design_name,".mat");
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='';               % either 'buried' or empty

n1=idx_layers;

lambda = 570*1e-9 +[-30 -15 0 15 30]*1e-9;

theta = real(asin(n1(end)/n1(1))/pi*180)+linspace(0,15,1e3);
        % start from critical angle

folder="SIM04_variability_study/";
f_ID=zeros(1,length(lambda));
for l=1:length(lambda)
    results_file  =strcat(folder,design_name,'_eff_index_values-LAMBDA',...
                    num2str(lambda(l)*1e9),'.csv');
    f_ID(l)=fopen(results_file,'w');
    fprintf(f_ID(l),"ID,neff_L',neff_L'',R_L,Q_L,neff_H',neff_H'',R_H,Q_H\n");
end
% fclose('all');

thickness_file=strcat(folder,design_name,'_thickness_values','.csv');
f2_ID=fopen(thickness_file,'w');
fprintf(f2_ID,"ID,thickness\n");
str_format=['%03i',repmat(',%+.2f',1,length(n1)),'\n'];
% fclose(f_ID);


N_tests = 2e3;
tic
for j = 1:N_tests
    
    % insert uncertainty in the thickness
    delta_thickness=1e-9*normrnd(0,sqrt(5),length(n1),1);
    d1=d_layers+delta_thickness;
    for l = 1:length(lambda)
%         f_ID=fopen(results_file(l),'w');
        neff_vect=[j,zeros(1,8)]; 
        for k = 1:2
            if strcmp(design_type,'buried')
                n = n1;
                d = d1;
                n(end-2)=n1(end-k);
            else
                n = [n1(1:end-k) ; n1(end)];
                d = [d1(1:end-k) ; d1(end)];
            end
            
            [dr,nr,~,~] = prepare_multilayer(d,n);
            
            R = reflectivity_par(lambda(l),theta,dr,nr,pol);
            
            [pks,idxs] = findpeaks(1-R);
            [~,pk_ix] = max(pks);
            idx = idxs(pk_ix);
            
            peak    = R(idx);
            min_s   = min(abs(R(1:idx)- peak/2));
            min_idx = find((abs(R(1:idx)- peak/2) == min_s));
            max_s   = min(abs(R(idx:end)- peak/2));
            max_idx = idx + find((abs(R(idx:end)- peak/2) == max_s));
            FWHM=abs(theta(min_idx) -theta(max_idx));
            Q = theta(idx)/FWHM;
            
            neff=sin(theta(idx)*pi/180)*n(1);
            PL=lambda(l)/neff/(pi*FWHM);
            
            neff_vect(2+4*(k-1))=neff;
            neff_vect(3+4*(k-1))=FWHM/2*neff;
            neff_vect(4+4*(k-1))=peak;
            neff_vect(5+4*(k-1))=Q;
%             col=['r','b'];
%             plot(theta,R,col(k));
%             drawnow
%             hold on
        end
        fprintf(f_ID(l),'%03i,%.4f,%.5f,%.4f,%05.0f,%.4f,%.5f,%.4f,%05.0f,\n',...
                neff_vect);
%         fclose(f_ID(l));
    end
%     pause(2)
%     hold off
    fprintf(f2_ID,str_format,j,delta_thickness*1e9);
    %     dlmwrite(results_file,neff_vect,'-append')
end
fclose('all');
toc
fclose('all');
%%