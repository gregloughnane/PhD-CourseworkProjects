clear all; close all; clc; format compact
% Force ~ Unif(375,625) lbs
% aparamF = 375; bF = 625;
aparamF = 300; bF = 700;
xF = 500;
PDF = unifpdf(xF,aparamF,bF);
CDF = unifcdf(xF,aparamF,bF);
sigmaF = normpdf(norminv(CDF)) / PDF;
muF = xF - norminv(CDF) * sigmaF; %% this is the equivalent normal for a and b
% Thickness ~ Unif(0.4,0.6) in
% aparamt = 0.4; bt = 0.6;
aparamt = 0.25; bt = 0.75;
xt = 0.5;
PDF = unifpdf(xt,aparamt,bt);
CDF = unifcdf(xt,aparamt,bt);
sigmat = normpdf(norminv(CDF)) / PDF;
mut = xt - norminv(CDF) * sigmat; %% this is the equivalent normal for a and b
% E (backboard,tempered glass) ~ N(10,1) Mpsi
muEb = 10e6 ; sigmaEb = 1e6;
% nu (backboard,tempered glass)~ N(0.22,0.11)
muvb = 0.22 ; sigmavb = 0.1;
% E (rim,steel) ~ LN(29,0.725) Mpsi
muEs = 29e6 ; sigmaEs = 0.7225e6;
% nu (rim,steel)~ N(0.29,0.145)
muvs = 0.29 ; sigmavs = 0.1;

%
limitstateconstraint = 26.8e3;
OrigInputFileName = 'Job-34.inp';
perturb = 0.01;

% Perturb Input Files
InpFilePerturbt(OrigInputFileName,'Job-34Perturbt.inp',perturb)
InpFilePerturbF(OrigInputFileName,'Job-34PerturbF.inp',perturb)
InpFilePerturbEb(OrigInputFileName,'Job-34PerturbEb.inp',perturb)
InpFilePerturbvb(OrigInputFileName,'Job-34Perturbvb.inp',perturb)
InpFilePerturbEs(OrigInputFileName,'Job-34PerturbEs.inp',perturb)
InpFilePerturbvs(OrigInputFileName,'Job-34Perturbvs.inp',perturb)

% Run Simulations
RunAbaqus('Job-34')
RunAbaqus('Job-34PerturbF')
RunAbaqus('Job-34PerturbEb')
RunAbaqus('Job-34Perturbvb')
RunAbaqus('Job-34PerturbEs')
RunAbaqus('Job-34Perturbvs'); pause(60)

% Write Stresses
WriteStressFromAbaqus = 'MacroWriteMaxStress.txt'; MakeNewMacro = 'NewMacro.py'; %% Define Stress Writing Macro Files
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34','Stresses')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34Perturbt','Stresses_t')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34PerturbF','Stresses_F')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34PerturbEb','Stresses_Eb')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34Perturbvb','Stresses_vb')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34PerturbEs','Stresses_Es')
system('abaqus cae noGUI=NewMacro.py'); pause(1)
PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34Perturbvs','Stresses_vs')
system('abaqus cae noGUI=NewMacro.py'); pause(1)

% Get FD Stress Derivatives
[MaxStress] = RptFileReader('Stresses.rpt')
[MaxStresst] = RptFileReader('Stresses_t.rpt')
[MaxStressF] = RptFileReader('Stresses_F.rpt')
[MaxStressEb] = RptFileReader('Stresses_Eb.rpt')
[MaxStressvb] = RptFileReader('Stresses_vb.rpt')
[MaxStressEs] = RptFileReader('Stresses_Es.rpt')
[MaxStressvs] = RptFileReader('Stresses_vs.rpt')
LS0 = limitstateconstraint - abs(MaxStress)
LSpt0 = limitstateconstraint - abs(MaxStresst)
LSpF0 = limitstateconstraint - abs(MaxStressF)
LSpEb0 = limitstateconstraint - abs(MaxStressEb)
LSpvb0 = limitstateconstraint - abs(MaxStressvb)
LSpEs0 = limitstateconstraint - abs(MaxStressEs)
LSpvs0 = limitstateconstraint - abs(MaxStressvs)
dLSdt0 = (LSpt0 - LS0) / abs((perturb*mut))
dLSdF0 = (LSpF0 - LS0) / abs((perturb*muF))
dLSdEb0 = (LSpEb0 - LS0) / abs((perturb*muEb))
dLSdvb0 = (LSpvb0 - LS0) / abs((perturb*muvb))
dLSdEs0 = (LSpEs0 - LS0) / abs((perturb*muEs))
dLSdvs0 = (LSpvs0 - LS0) / abs((perturb*muvs))

% Calculate initial beta, alpha, and new design points
sigmahat = sqrt( (dLSdt0*sigmat)^2 + (dLSdF0*sigmaF)^2 + (dLSdEb0*sigmaEb)^2 ...
                + (dLSdvb0*sigmavb)^2 + (dLSdEs0*sigmaEs)^2 + (dLSdvs0*sigmavs)^2)
b0 = LS0 / sigmahat
at0 = -dLSdt0*sigmat / sigmahat
aF0 = -dLSdF0*sigmaF / sigmahat
aEb0 = -dLSdEb0*sigmaEb / sigmahat
avb0 = -dLSdvb0*sigmavb / sigmahat
aEs0 = -dLSdEs0*sigmaEs / sigmahat
avs0 = -dLSdvs0*sigmavs / sigmahat
xt0 = mut + b0 * at0 * sigmat
xF0 = muF + b0 * aF0 * sigmaF
xEb0 = muEb + b0 * aEb0 * sigmaEb
xvb0 = muvb + b0 * avb0 * sigmavb
xEs0 = muEs + b0 * aEs0 * sigmaEs
xvs0 = muvs + b0 * avs0 * sigmavs
ut0 = (xt0 - mut) / sigmat
uF0 = (xF0 - muF) / sigmaF
uEb0 = (xEb0 - muEb) / sigmaEb
uvb0 = (xvb0 - muvb) / sigmavb
uEs0 = (xEs0 - muEs) / sigmaEs
uvs0 = (xvs0 - muvs) / sigmavs

for i = 1:20
    
    if i == 1
        ut(i) = ut0
        uF(i) = uF0
        uEb(i) = uEb0
        uvb(i) = uvb0
        uEs(i) = uEs0
        uvs(i) = uvs0
        
        xt(i) = xt0
        xF(i) = xF0
        xEb(i) = xEb0
        xvb(i) = xvb0
        xEs(i) = xEs0
        xvs(i) = xvs0
    else
        ut(i) = b(i-1) * at(i-1)
        uF(i) = b(i-1) * aF(i-1)
        uEb(i) = b(i-1) * aEb(i-1)
        uvb(i) = b(i-1) * avb(i-1)
        uEs(i) = b(i-1) * aEs(i-1)
        uvs(i) = b(i-1) * avs(i-1)
        
        xt(i) = mut(i-1) + ut(i) * sigmat(i-1)
        xF(i) = muF(i-1) + uF(i) * sigmaF(i-1)
        xEb(i) = muEb + uEb(i) * sigmaEb
        xvb(i) = muvb + uvb(i) * sigmavb
        xEs(i) = muEs + uEs(i) * sigmaEs
        xvs(i) = muvs + uvs(i) * sigmavs
    end
    
    % Force ~ Unif(375,625) lbs
    PDF = unifpdf(xF(i),aparamF,bF);
    CDF = unifcdf(xF(i),aparamF,bF);
    sigmaF(i) = normpdf(norminv(CDF)) / PDF;
    muF(i) = xF(i) - norminv(CDF) * sigmaF(i); %% this is the equivalent normal for a and b
    % Thickness ~ Unif(0.4,0.6) in
    PDF = unifpdf(xt(i),aparamt,bt);
    CDF = unifcdf(xt(i),aparamt,bt);
    sigmat(i) = normpdf(norminv(CDF)) / PDF;
    mut(i) = xt(i) - norminv(CDF) * sigmat(i); %% this is the equivalent normal for a and b
    
    % Create New Input File & Perturb
    InpFileCreator('Job-34.inp','Job-34New.inp',xt(i),-xF(i),xEb(i),xvb(i),xEs(i),xvs(i))
    InpFilePerturbtLOOP('Job-34New.inp','PerttNew.inp',perturb)
    InpFilePerturbF('Job-34New.inp','PertFNew.inp',perturb)
    InpFilePerturbEb('Job-34New.inp','PertEbNew.inp',perturb)
    InpFilePerturbvb('Job-34New.inp','PertvbNew.inp',perturb)
    InpFilePerturbEs('Job-34New.inp','PertEsNew.inp',perturb)
    InpFilePerturbvs('Job-34New.inp','PertvsNew.inp',perturb)
    % Run Simulations
    RunAbaqus('Job-34New');
    RunAbaqus('PerttNew'); 
    RunAbaqus('PertFNew');
    RunAbaqus('PertEbNew'); 
    RunAbaqus('PertvbNew'); 
    RunAbaqus('PertEsNew'); 
    RunAbaqus('PertvsNew'); pause(60)
    
    % Write Stresses
    WriteStressFromAbaqus = 'MacroWriteMaxStress.txt'; MakeNewMacro = 'NewMacro.py';
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'Job-34New','Stresses')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PerttNew','Stresses_t')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PertFNew','Stresses_F')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PertEbNew','Stresses_Eb')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PertvbNew','Stresses_vb')
    system('abaqus cae noGUI=NewMacro.py'); pause(1);     
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PertEsNew','Stresses_Es')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,'PertvsNew','Stresses_vs')
    system('abaqus cae noGUI=NewMacro.py'); pause(1); 
    
    % Get FD Stress Derivatives
    [MaxStress(i)] = RptFileReader('Stresses.rpt')
    [MaxStresst(i)] = RptFileReader('Stresses_t.rpt')
    [MaxStressF(i)] = RptFileReader('Stresses_F.rpt')
    [MaxStressEb(i)] = RptFileReader('Stresses_Eb.rpt')
    [MaxStressvb(i)] = RptFileReader('Stresses_vb.rpt')
    [MaxStressEs(i)] = RptFileReader('Stresses_Es.rpt')
    [MaxStressvs(i)] = RptFileReader('Stresses_vs.rpt')
    LS(i) = limitstateconstraint - abs(MaxStress(i))
    LSpt(i) = limitstateconstraint - abs(MaxStresst(i))
    LSpF(i) = limitstateconstraint - abs(MaxStressF(i))
    LSpEb(i) = limitstateconstraint - abs(MaxStressEb(i))
    LSpvb(i) = limitstateconstraint - abs(MaxStressvb(i))
    LSpEs(i) = limitstateconstraint - abs(MaxStressEs(i))
    LSpvs(i) = limitstateconstraint - abs(MaxStressvs(i))
    dLSdt(i) = (LSpt(i) - LS(i)) / abs((perturb*xt(i)))
    dLSdF(i) = (LSpF(i) - LS(i)) / abs((perturb*xF(i)))
    dLSdEb(i) = (LSpEb(i) - LS(i)) / abs((perturb*xEb(i)))
    dLSdvb(i) = (LSpvb(i) - LS(i)) / abs((perturb*xvb(i)))
    dLSdEs(i) = (LSpEs(i) - LS(i)) / abs((perturb*xEs(i)))
    dLSdvs(i) = (LSpvs(i) - LS(i)) / abs((perturb*xvs(i)))
    
    % Calculate alpha & beta values, and normalized design points    
    sigmahat(i) = sqrt( (dLSdt(i)*sigmat(i))^2 + (dLSdF(i)*sigmaF(i))^2 + (dLSdEb(i)*sigmaEb)^2 ...
                + (dLSdvb(i)*sigmavb)^2 + (dLSdEs(i)*sigmaEs)^2 + (dLSdvs(i)*sigmavs)^2)
    b(i) = (LS(i) - ( dLSdt(i)*sigmat(i)*ut(i) + dLSdF(i)*sigmaF(i)*uF(i) + dLSdEb(i)*sigmaEb*uEb(i) ...
        + dLSdvb(i)*sigmavb*uvb(i) + dLSdEs(i)*sigmaEs*uEs(i) + dLSdvs(i)*sigmavs*uvs(i)))/ sigmahat(i)
    at(i) = -dLSdt(i)*sigmat(i) / sigmahat(i)
    aF(i) = -dLSdF(i)*sigmaF(i) / sigmahat(i)
    aEb(i) = -dLSdEb(i)*sigmaEb / sigmahat(i)
    avb(i) = -dLSdvb(i)*sigmavb / sigmahat(i)
    aEs(i) = -dLSdEs(i)*sigmaEs / sigmahat(i)
    avs(i) = -dLSdvs(i)*sigmavs / sigmahat(i)
    
    % Calculate percent error
    if i == 1
        e(i) = abs(( b(i) - b0 ) / b0)
    else
        e(i) = abs(( b(i) - b(i-1) ) / b(i-1))
        if e(i) <= 1e-3
            disp('CONVERGED!')
            break
        end    
    end
    
end

LimitState = [LS0;LS']
b = [b0;b']
e = [0;e']
xt = [xt0;xt']
xF = [xF0;xF']
xEb = [xEb0;xEb']
xvb = [xvb0;xvb']
xEs = [xEs0;xEs']
xvs = [xvs0;xvs']
dLSdt = [dLSdt0;dLSdt']
dLSdF = [dLSdF0;dLSdF']
dLSdEb = [dLSdEb0;dLSdEb']
dLSdvb = [dLSdvb0;dLSdvb']
dLSdEs = [dLSdEs0;dLSdEs']
dLSdvs = [dLSdvs0;dLSdvs']
ut = [ut0;ut']
uF = [uF0;uF']
uEb = [uEb0;uEb']
uvb = [uvb0;uvb']
uEs = [uEs0;uEs']
uvs = [uvs0;uvs']
at = [at0;at']
aF = [aF0;aF']
aEb = [aEb0;aEb']
avb = [avb0;avb']
aEs = [aEs0;aEs']
avs = [avs0;avs']




