clear all; close all; clc; format compact
% Force ~ Unif(250,1000) lbs
aF = 250; bF = 750;
force = 500;
muF = 500; sigmaF = 100;
% Thickness ~ Unif(0.4,0.6) in
at = 0.4; bt = 0.6;
thickness = 0.5;
mut = 0.5; sigmat = 0.0798;
% E (backboard,tempered glass) ~ N(10,1) Mpsi
muEb = 10e6 ; sigmaEb = 1e6;
% nu (backboard,tempered glass)~ N(0.22,0.11)
muvb = 0.22 ; sigmavb = 0.11;
% E (rim,steel) ~ LN(29,0.725) Mpsi
muEs = 29e6 ; sigmaEs = 0.725e6;
% nu (rim,steel)~ N(0.29,0.145)
muvs = 0.29 ; sigmavs = 0.145;

%% Randomly Sample Inputs
% [f,E_bb,nu_bb,thickness,E_r,nu_r] = GenerateInputs();

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
LS = limitstateconstraint - abs(MaxStress)
LSpt = limitstateconstraint - abs(MaxStresst)
LSpF = limitstateconstraint - abs(MaxStressF)
LSpEb = limitstateconstraint - abs(MaxStressEb)
LSpvb = limitstateconstraint - abs(MaxStressvb)
LSpEs = limitstateconstraint - abs(MaxStressEs)
LSpvs = limitstateconstraint - abs(MaxStressvs)
dLSdt = (LSpt - LS) / abs((perturb*mut))
dLSdF = (LSpF - LS) / abs((perturb*muF))
dLSdEb = (LSpEb - LS) / abs((perturb*muEb))
dLSdvb = (LSpvb - LS) / abs((perturb*muvb))
dLSdEs = (LSpEs - LS) / abs((perturb*muEs))
dLSdvs = (LSpvs - LS) / abs((perturb*muvs))

% Calculate initial beta, alpha, and new design points
%%% Compute equivalent normal means and standard deviations
%%%% force
% PDF = unifpdf(init_force,aF,bF)
% CDF = 

sigmahat = sqrt( (dLSdt*sigmat)^2 + (dLSdF*sigmaF)^2 + (dLSdEb*sigmaEb)^2 ...
                + (dLSdvb*sigmavb)^2 + (dLSdEs*sigmaEs)^2 + (dLSdvs*sigmavs)^2)
b0 = LS / sigmahat
at0 = -dLSdt*sigmat / sigmahat
aF0 = -dLSdF*sigmaF / sigmahat
aEb0 = -dLSdEb*sigmaEb / sigmahat
avb0 = -dLSdvb*sigmavb / sigmahat
aEs0 = -dLSdEs*sigmaEs / sigmahat
avs0 = -dLSdvs*sigmavs / sigmahat
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

for i = 1:10
    
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
        
        xt(i) = mut + ut(i) * sigmat
        xF(i) = muF + uF(i) * sigmaF
        xEb(i) = muEb + uEb(i) * sigmaEb
        xvb(i) = muvb + uvb(i) * sigmavb
        xEs(i) = muEs + uEs(i) * sigmaEs
        xvs(i) = muvs + uvs(i) * sigmavs
    end
    
    % Create New Input File & Perturb
    InpFileCreator('Job-34.inp','Job-34New.inp',xt(i),-xF(i),xEb(i),xvb(i),xEs(i),xvs(i))
    InpFilePerturbt('Job-34New.inp','PerttNew.inp',perturb)
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
    sigmahat(i) = sqrt( (dLSdt(i)*sigmat)^2 + (dLSdF(i)*sigmaF)^2 + (dLSdEb(i)*sigmaEb)^2 ...
                + (dLSdvb(i)*sigmavb)^2 + (dLSdEs(i)*sigmaEs)^2 + (dLSdvs(i)*sigmavs)^2)
    b(i) = (LS(i) - ( dLSdt(i)*sigmat*ut(i) + dLSdF(i)*sigmaF*uF(i) + dLSdEb(i)*sigmaEb*uEb(i) ...
        + dLSdvb(i)*sigmavb*uvb(i) + dLSdEs(i)*sigmaEs*uEs(i) + dLSdvs(i)*sigmavs*uvs(i)))/ sigmahat(i)
    at(i) = -dLSdt(i)*sigmat / sigmahat(i)
    aF(i) = -dLSdF(i)*sigmaF / sigmahat(i)
    aEb(i) = -dLSdEb(i)*sigmaEb / sigmahat(i)
    avb(i) = -dLSdvb(i)*sigmavb / sigmahat(i)
    aEs(i) = -dLSdEs(i)*sigmaEs / sigmahat(i)
    avs(i) = -dLSdvs(i)*sigmavs / sigmahat(i)
    
    % Calculate percent error
    if i == 1
        e(i) = abs(( b(i) - b0 ) / b0)
    else
        e(i) = abs(( b(i) - b(i-1) ) / b(i-1))
        if e(i) <= 0.001
            disp('CONVERGED!')
            break
        end    
    end
    
end

LimitState = [LS;LS(i)']
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




