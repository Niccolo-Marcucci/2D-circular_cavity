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

design_name = "TM_Descrovi";%gd3_buriedDBR";
design_file = strcat("Lumerical-Objects/multilayer_design/designs/design_",...
    design_name,".mat");
load(design_file);
pol='p';                            % polarisation: 'p' or 's'
design_type='';               % either 'buried' or empty

n1=idx_layers;

lambda = 570*1e-9 +[-30 -15 0 15 30]*1e-9;

if strcmp(design_type,'buried')
    theta = linspace(50,70,5e4);
else
    theta = real(asin(n1(end)/n1(1))/pi*180)+linspace(0,20,1e3);
    % start from critical angle
end
folder="SIM04_variability_study/";
f_ID=zeros(1,length(lambda));
for i=1:length(lambda)
    results_file  =strcat(folder,design_name,'_eff_index_values-LAMBDA',...
                    num2str(lambda(i)*1e9),'.csv');
    f_ID(i)=fopen(results_file,'w');
    fprintf(f_ID(i),"ID,neff_L',neff_L'',R_L,Q_L,neff_H',neff_H'',R_H,Q_H\n");
end
% fclose('all');

thickness_file=strcat(folder,design_name,'_thickness_values','.csv');
f2_ID=fopen(thickness_file,'w');
fprintf(f2_ID,"ID,thickness\n");
str_format=['%03i',repmat(',%+.2f',1,length(n1)),'\n'];
% fclose(f_ID);


N_tests = 2e3;
t=0;
j=0;
while  j < N_tests
    tic
    j = j+1;
    
    % insert uncertainty in the thickness
    delta_thickness=1e-9*normrnd(0,sqrt(5),length(n1),1);
    d1=d_layers+delta_thickness;
    for i = 1:length(lambda)
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
            
            R = reflectivity_par(lambda(i),theta,dr,nr,pol);
            
            [pks,idxs] = findpeaks(1-R);
            % in buried design the pseudo-BSW peak might be close to some 
            % bulk modes, hence in order to find the on relative to the BSW
            % we check the field distribution. Idf the macimum is close to
            % the external interface, it cannot be a bulky mode
            for  peak_idx = idxs 
                [z1, nz, P] = field_distribution(lambda(i),theta(peak_idx),d,n,pol);
                [~, max_loc] = max(P);
%                 plot(z1,P,z1,nz/max(nz)*max(P))
%                 hold on;
%                 plot(z1(max_loc)*[1,1],[0 max(P)],'k')
%                 plot((sum(d(2:length(d1)-3))-d(length(d1)-3)/3)*[1,1],[0 max(P)],'r')
                if z1(max_loc) > (sum(d(2:length(d1)-3))-d(length(d1)-3)/3)
                    break 
                else
                    peak_idx = NaN;
                end
            end
%             [~,pk_ix] = max(pks);
%             idx = idxs(pk_ix);
            
            if isnan(peak_idx)
                j=j-1;
                break
            end
            col=['r','b'];
            plot(theta,R,col(k),theta(peak_idx),R(peak_idx),'k+');
            drawnow
            hold on
            peak = R(peak_idx);
%             min_s   = min(abs(R(1:idx)- peak/2));
%             min_idx = find((abs(R(1:idx)- peak/2) == min_s));
%             max_s   = min(abs(R(idx:end)- peak/2));
%             max_idx = idx + find((abs(R(idx:end)- peak/2) == max_s));
            min_s = abs(diff(R(1:peak_idx)-(peak+(1-peak)/2) > 0)); 
                    %equal to one when R(1:idx)- peak/2  changes sign
            min_idx = find(min_s);
            max_s = abs(diff(R(peak_idx:end)-(peak+(1-peak)/2) > 0));
            max_idx = find(max_s);
            
            if isempty(min_idx) || isempty(max_idx)
                % repeat this j if the mode is not found
                j=j-1;
                break
            end
            FWHM=abs(theta(min_idx(end))-theta(peak_idx+max_idx(1)));
            Q = theta(peak_idx)/FWHM;
            
            neff=sin(theta(peak_idx)*pi/180)*n(1);
%             PL=lambda(i)/neff/(pi*FWHM);
            
            neff_vect(2+4*(k-1))=neff;
            neff_vect(3+4*(k-1))=FWHM/2*neff;
            neff_vect(4+4*(k-1))=peak;
            neff_vect(5+4*(k-1))=Q;
        end
        if isnan(peak_idx) || isempty(min_idx) || isempty(max_idx)
            % if i exited from previous loop with break, then repeat this j
            break
        end
        fprintf(f_ID(i),'%03i,%.4f,%.5f,%.4f,%05.0f,%.4f,%.5f,%.4f,%05.0f,\n',...
                            neff_vect);
    end
    if isnan(peak_idx) || isempty(min_idx) || isempty(max_idx)
        % if i exited from previous loop with break, then repeat this j
        continue
    end
    fprintf(f2_ID,str_format,j,delta_thickness*1e9);
    %     dlmwrite(results_file,neff_vect,'-append')
    t = t + toc;
    if mod(j,10) == 0
        display(['index j = ',num2str(j)])
        display(['time left is ',num2str(t/j*(N_tests-j)/60)])
    end
end
fclose('all');

fclose('all');
%%