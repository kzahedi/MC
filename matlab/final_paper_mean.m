%%

file_path = '/Users/zahedi/bwSyncAndShare/Article_QuantMorphComp/Model_Results_Data/';

suffix = '_noPerturbation';
files(1).name = [file_path '/NoPerturbation/' 'results_motor' suffix];
files(2).name = [file_path '/NoPerturbation/' 'results_muscle_linFv_constFl_FFB' suffix];
files(3).name = [file_path '/NoPerturbation/' 'results_muscle_HillFv_HillFl_FFB' suffix];
% files(4).name = [file_path '/NoPerturbation/' 'results_muscle_linFv_linFl_FFB' suffix];
% files(5).name = [file_path '/NoPerturbation/' 'results_muscle_HillFv_constFl_FFB' suffix];
% files(6).name = [file_path '/MuscleTendonModel/' 'results_muscle_MTC_FFB' suffix];

%% 1. Step: Extract min/max values for each column that will be used
p_min  = 0; p_max = 0;
ac_min = 0; ac_max = 0;
v_min  = 0; v_max = 0;
a_min  = 0; a_max = 0;
s_min  = 0; s_max = 0;

domain_string = '';


for file_index = 1:length(files)
    load(files(file_index).name,'SimData');
    
    position      = SimData(:,2);
    velocity      = SimData(:,3);
    accelaration  = SimData(:,4);
    action        = SimData(:,10);
    muscle_sensor = SimData(:,5);
    
    if file_index == 1
        p_min = min(position);
        p_max = max(position);
        
        v_min = min(velocity);
        v_max = max(velocity);
        
        a_min = min(action);
        a_max = max(action);
        
        s_min = min(muscle_sensor);
        s_max = max(muscle_sensor);
        
        ac_min = min(accelaration);
        ac_max = max(accelaration);
    else
        p_min = min(p_min, min(position));
        p_max = max(p_max, max(position));
        
        v_min = min(v_min, min(velocity));
        v_max = max(v_max, max(velocity));
        
        a_min = min(a_min, min(action));
        a_max = max(a_max, max(action));
        
        ac_min = min(ac_min, min(accelaration));
        ac_max = max(ac_max, max(accelaration));
        
        s_min = min(s_min, min(muscle_sensor));
        s_max = max(s_max, max(muscle_sensor));
    end
end
domain_string = [domain_string 'Domains\n' ...
    sprintf('  Position:        %f, %f\n', p_min, p_max) ...
    sprintf('  Velocity:        %f, %f\n', v_min, v_max) ...
    sprintf('  Accelaration:    %f, %f\n', ac_min, ac_max) ...
    sprintf('  Actuator signal: %f, %f\n', a_min, a_max)];

%%
for bins = [100]
    
    % 2. Step: Discretise data
    w_bins = bins;
    a_bins = bins;
    s_bins = bins;
    
    %%
    %fprintf('Binning the data with |W| = %d, |A| = %d, |S| = %d\n', w_bins, a_bins, s_bins);
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        % fprintf('Discretising file %s\n', name);
        load(files(file_index).name, 'SimData');
        
        position      = SimData(:,2);
        velocity      = SimData(:,3);
        accelaration  = SimData(:,4);
        action        = SimData(:,10);
        muscle_sensor = SimData(:,5);
        
        files(file_index).d_position     = discretiseMatrix(position, p_min, p_max, w_bins);
        files(file_index).d_velocity     = discretiseMatrix(velocity, v_min, v_max, w_bins);
        files(file_index).d_accelaration = discretiseMatrix(accelaration, ac_min, ac_max, w_bins);
        
        files(file_index).a              = discretiseMatrix(action, a_min, a_max, a_bins);
        files(file_index).w              = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity, files(file_index).d_accelaration]);
        
        if isempty(strfind(files(file_index).name,'muscle'))
            fprintf('%s is not a muscle\n', name);
            files(file_index).s          = combineAndRelabelBinnedMatrix([files(file_index).d_position, files(file_index).d_velocity]);
        else
            fprintf('%s is a muscle\n', name);
            files(file_index).s          = discretiseMatrix(muscle_sensor, s_min, s_max, s_bins);
            
        end
    end
    
    %% calculate the measures
    
    %fprintf('Starting calculations\n')
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        %fprintf('Working on file %s\n', name);
        w2 = files(file_index).w(2:end,:);
        w1 = files(file_index).w(1:end-1,:);
        a1 = files(file_index).a(1:end-1,:);
        s1 = files(file_index).s(1:end-1,:);
        
        
        fprintf('Calculating MC_MI\n')
        files(file_index).mcmi = MC_MI(w2, w1, s1, a1);
        fprintf('  Result %f\n', files(file_index).mcmi);
        fprintf('Calculating MC_W\n')
        files(file_index).mcw = MC_W(w2, w1, a1);
        fprintf('  Result %f\n', files(file_index).mcw);
        %fprintf('Calculating MC_C\n')
        %files(file_index).mcc = MC_C(w2, w1, s1, a1);
        %fprintf('  Result %f\n', files(file_index).mcc);
    end
    
    %% log to console
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        fprintf('Filename %s\n', name);
        fprintf('  MC_W:  %f\n', files(file_index).mcw);
        fprintf('  MC_MI: %f\n', files(file_index).mcmi);
        %fprintf('  MC_C:  %f\n', files(file_index).mcc);
    end
    
    %% log to file
    %filename = sprintf('/Users/zahedi/bwSyncAndShare/Article_QuantMorphComp/results/no pertubation/global_binning_results_w%d_a%d_s%d%s_%dbins.txt', w_bins, a_bins, s_bins, suffix, bins);
    filename = sprintf('/Users/zahedi/projects/QuantMorphComp/data/global_binning_results_w%d_a%d_s%d%s_%dbins.txt', w_bins, a_bins, s_bins, suffix, bins);
    fprintf('Writing results to %s\n', filename);
    
    fileID = fopen(filename,'w');
    fprintf(fileID, domain_string);
    fprintf(fileID, 'Bins %d\n', bins);
    for file_index = 1:length(files)
        [pathstr,name,ext] = fileparts(files(file_index).name);
        fprintf(fileID, 'Filename %s\n', name);
        fprintf(fileID, '  MC_W:  %f\n', files(file_index).mcw);
        fprintf(fileID, '  MC_MI: %f\n', files(file_index).mcmi);
        %fprintf(fileID, '  MC_C:  %f\n', files(file_index).mcc);
    end
    fclose(fileID);
    fprintf('done.\n')
end