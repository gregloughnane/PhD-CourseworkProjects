clear all; close all; clc; format compact
% Force ~ Unif(250,1000) lbs
a = 250; b = 1000;
force = 500;
muF = 500; sigmaF = 25;
% Thickness ~ Unif(0.4,0.6) in
a = 0.4; b = 0.6;
thickness = 0.5;
mut = 0.5; sigmat = 0.025;
% E (backboard,tempered glass) ~ N(10,1) Mpsi
muEb = 10e6 ; sigmaEb = 1e6;
% nu (backboard,tempered glass)~ N(0.22,0.11)
muvb = 0.22 ; sigmavb = 0.11;
% E (rim,steel) ~ LN(29,0.725) Mpsi
muEs = 29e6 ; sigmaEs = 0.725e6;
% nu (rim,steel)~ N(0.29,0.145)
muvs = 0.29 ; sigmavs = 0.145;

perturb = 0.01;

% %% Randomly Sample Inputs
% % [f,E_bb,nu_bb,thickness,E_r,nu_r] = GenerateInputs();
% 
% %% Limit State Constraint
% limitstateconstraint = 26.8e3;
% 
% %% Define Stress Writing Macro Files
% WriteStressFromAbaqus = 'MacroWriteMaxStress.txt';
% MakeNewMacro = 'NewMacro.py';
% 
% %% Get Stress From Original Model
% JobFileName = 'Job-34';
% RptFileName = 'Stresses';
% RptFile = 'Stresses.rpt';
% % RunAbaqus(JobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,JobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LS = limitstateconstraint - abs(MaxStress)
% 
% %% Get FD Stress Derivative Info
% OrigInputFileName = 'Job-34.inp';
% 
% % Perturb thickness
% RptFileName = 'Stresses_t';
% RptFile = 'Stresses_t.rpt';
% NewInputFile = 'Job-34Perturbt.inp';
% NewJobFileName = 'Job-34Perturbt';
% InpFilePerturbt(OrigInputFileName,NewInputFile,perturb,thickness)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpt = limitstateconstraint - abs(MaxStress)
% dLSdt = (LSpt - LS) / abs((perturb*thickness))
% 
% % Perturb Force
% RptFileName = 'Stresses_F';
% RptFile = 'Stresses_F.rpt';
% NewInputFile = 'Job-34PerturbF.inp';
% NewJobFileName = 'Job-34PerturbF';
% [force,newforce] = InpFilePerturbF(OrigInputFileName,NewInputFile,perturb)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpF = limitstateconstraint - abs(MaxStress)
% dLSdF = (LSpF - LS) / abs((perturb*force))
% 
% % Perturb E (backboard, tempered glass)
% RptFileName = 'Stresses_Eb';
% RptFile = 'Stresses_Eb.rpt';
% NewInputFile = 'Job-34PerturbEb.inp';
% NewJobFileName = 'Job-34PerturbEb';
% [Eb,newEb] = InpFilePerturbEb(OrigInputFileName,NewInputFile,perturb)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpEb = limitstateconstraint - abs(MaxStress)
% dLSdEb = (LSpEb - LS) / abs((perturb*Eb))
% 
% % Perturb nu (backboard, tempered glass)
% RptFileName = 'Stresses_vb';
% RptFile = 'Stresses_vb.rpt';
% NewInputFile = 'Job-34Perturbvb.inp';
% NewJobFileName = 'Job-34Perturbvb';
% [vb,newvb] = InpFilePerturbvb(OrigInputFileName,NewInputFile,perturb)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpvb = limitstateconstraint - abs(MaxStress)
% dLSdvb = (LSpvb - LS) / abs((perturb*vb))
% 
% % Perturb E (rim, steel)
% RptFileName = 'Stresses_Es';
% RptFile = 'Stresses_Es.rpt';
% NewInputFile = 'Job-34PerturbEs.inp';
% NewJobFileName = 'Job-34PerturbEs';
% [Es,newEs] = InpFilePerturbEs(OrigInputFileName,NewInputFile,perturb)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpEs = limitstateconstraint - abs(MaxStress)
% dLSdEs = (LSpEs - LS) / abs((perturb*Es))
% 
% % Perturb nu (rim, steel)
% RptFileName = 'Stresses_vs';
% RptFile = 'Stresses_vs.rpt';
% NewInputFile = 'Job-34Perturbvs.inp';
% NewJobFileName = 'Job-34Perturbvs';
% [vs,newvs] = InpFilePerturbvs(OrigInputFileName,NewInputFile,perturb)
% RunAbaqus(NewJobFileName); pause(45)
% PythonFileChanger(WriteStressFromAbaqus,MakeNewMacro,NewJobFileName,RptFileName)
% system('abaqus cae noGUI=NewMacro.py'); pause(5); 
% [MaxStress] = RptFileReader(RptFile)
% LSpvs = limitstateconstraint - abs(MaxStress)
% dLSdvs = (LSpvs - LS) / abs((perturb*vs))


% stressderivatives = [dLSdt dLSdF dLSdEb dLSdvb dLSdEs dLSdvs]'

load HL_init.mat
% Calculate initial beta, alpha, and new design points
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
% 
% 
% for i = 1:10
%     
%     if i == 1
%         xE(i) = xE0
%         xF(i) = xF0
%         xA(i) = xA0
%         xL(i) = xL0
%         xI(i) = xI0
%         
%         uE(i) = uE0
%         uF(i) = uF0
%         uA(i) = uA0
%         uL(i) = uL0
%         uI(i) = uI0
%     elseif i == 2
%         xE(i) = xE1
%         xF(i) = xF1
%         xA(i) = xA1
%         xL(i) = xL1
%         xI(i) = xI1
%         
%         uE(i) = x(1)
%         uF(i) = x(2)
%         uA(i) = x(3)
%         uL(i) = x(4)
%         uI(i) = x(5)
%     elseif i == 4
%         xE(i) = xE3
%         xF(i) = xF3
%         xA(i) = xA3
%         xL(i) = xL3
%         xI(i) = xI3
%         
%         uE(i) = x(1)
%         uF(i) = x(2)
%         uA(i) = x(3)
%         uL(i) = x(4)
%         uI(i) = x(5)
%     else
%         uE(i) = b(i-1) * aE(i-1)
%         uF(i) = b(i-1) * aF(i-1)
%         uA(i) = b(i-1) * aA(i-1)
%         uL(i) = b(i-1) * aL(i-1)
%         uI(i) = b(i-1) * aI(i-1)
%         xE(i) = muE + uE(i) * sigmaE
%         xF(i) = muF + uF(i) * sigmaF
%         xA(i) = muA + uA(i) * sigmaA
%         xL(i) = muL + uL(i) * sigmaL
%         xI(i) = muI + uI(i) * sigmaI
%     end
%     
%     % Create new input file, perturb variables and run abaqus FE simulations
%     InpFileCreator('6BarFrame.inp','6BarFrameNew.inp',xE(i),-xF(i),xA(i),xL(i),xI(i))
%     [E(i),newE(i)] = InpFilePerturbE('6BarFrameNew.inp','PertENew.inp',perturb)
%     [F(i),newF(i)] = InpFilePerturbF('6BarFrameNew.inp','PertFNew.inp',perturb)
%     [A(i),newA(i)] = InpFilePerturbA('6BarFrameNew.inp','PertANew.inp',perturb)
%     [L(i),newL(i)] = InpFilePerturbL('6BarFrameNew.inp','PertLNew.inp',perturb)
%     [I(i),newI(i)] = InpFilePerturbI('6BarFrameNew.inp','PertINew.inp',perturb)
%     RunAbaqus('6BarFrameNew')
%     RunAbaqus('PertENew'); 
%     RunAbaqus('PertFNew');
%     RunAbaqus('PertANew'); 
%     RunAbaqus('PertINew'); 
%     RunAbaqus('PertLNew'); pause(45)
%     % Calculate displacements, corresponding limit states, and derivatives 
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','6BarFrameNew','DisplacementatNode6New');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [gx(i)] = RptFileReader('DisplacementatNode6New.rpt')
%     LimitState(i) = limitstateconstraint - abs(gx(i))
%     % Perturb E
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','PertENew','DisplacementatNode6PertENew');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [d_pertE(i)] = RptFileReader('DisplacementatNode6PertENew.rpt')
%     LimitState_pertE(i) = limitstateconstraint - abs(d_pertE(i))
%     dDdE(i) = (LimitState_pertE(i) - LimitState(i)) / ((perturb-1)*E(i))
%     % Perturb F
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','PertFNew','DisplacementatNode6PertFNew');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [d_pertF(i)] = RptFileReader('DisplacementatNode6PertFNew.rpt')
%     LimitState_pertF(i) = limitstateconstraint - abs(d_pertF(i))
%     dDdF(i) = (LimitState_pertF(i) - LimitState(i)) / ((perturb-1)*-F(i))
%     % Perturb A
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','PertANew','DisplacementatNode6PertANew');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [d_pertA(i)] = RptFileReader('DisplacementatNode6PertANew.rpt')
%     LimitState_pertA(i) = limitstateconstraint - abs(d_pertA(i))
%     dDdA(i) = (LimitState_pertA(i) - LimitState(i)) / ((perturb-1)*A(i))
%     % Perturb L
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','PertLNew','DisplacementatNode6PertLNew');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [d_pertL(i)] = RptFileReader('DisplacementatNode6PertLNew.rpt')
%     LimitState_pertL(i) = limitstateconstraint - abs(d_pertL(i))
%     dDdL(i) = (LimitState_pertL(i) - LimitState(i)) / ((perturb-1)*L(i))
%     % Perturb I
%     PythonFileChanger('abaqusMacros.txt','abaqusChanged.py','PertINew','DisplacementatNode6PertINew');
%     system('abaqus cae noGUI=abaqusChanged.py'); pause(5)
%     [d_pertI(i)] = RptFileReader('DisplacementatNode6PertINew.rpt')
%     LimitState_pertI(i) = limitstateconstraint - abs(d_pertI(i))
%     dDdI(i) = (LimitState_pertI(i) - LimitState(i)) / ((perturb-1)*I(i))
%     
%     sigmahat(i) = sqrt( (dDdE(i)*sigmaE)^2 + (dDdF(i)*sigmaF)^2 + (dDdA(i)*sigmaA)^2 ...
%             + (dDdL(i)*sigmaL)^2 + (dDdI(i)*sigmaI)^2 )
% 
%     if i == 1
%         % TANA
%         % Solve Nonlinearity Index r
%         options = optimset('Display','off','TolFun',1e-10); x0 = 1;
%         f = @(r)(gx0 - (LimitState(i) +(1/r)*((xE(i)^(1-r)*dDdE(i)*(muE^r - xE(i)^r)) + ...
%                                               ((-xF(i))^(1-r)*dDdF(i)*(muF^r - (-xF(i))^r)) + ...
%                                               (xA(i)^(1-r)*dDdA(i)*(muA^r - xA(i)^r)) + ...
%                                               (xL(i)^(1-r)*dDdL(i)*(muL^r - xL(i)^r)) + ...
%                                               (xI(i)^(1-r)*dDdI(i)*(muI^r - xI(i)^r))) ) )
%         [r,fval] = fsolve(f,x0,options)
%         % Plug r into TANA expansion equation to determine g(x)
%         syms x1 x2 x3 x4 x5 gTANA
%         limitstate = sym(LimitState(i),'f')
%         xe = sym(xE(i),'f'); dgde = sym(dDdE(i),'f')
%         xf = sym(-xF(i),'f'); dgdf = sym(dDdF(i),'f')
%         xa = sym(xA(i),'f'); dgda = sym(dDdA(i),'f')
%         xl = sym(xL(i),'f'); dgdl = sym(dDdL(i),'f')
%         xi = sym(xI(i),'f'); dgdi = sym(dDdI(i),'f')
%         gTANA = symfun(limitstate + (1/r)*((xe^(1-r)*dgde*(x1^r - xe^r)) + ...
%                                            (xf^(1-r)*dgdf*(x2^r - xf^r)) + ...
%                                            (xa^(1-r)*dgda*(x3^r - xa^r)) + ...
%                                            (xl^(1-r)*dgdl*(x4^r - xl^r)) + ...
%                                            (xi^(1-r)*dgdi*(x5^r - xi^r)) ),...
%                                            [x1 x2 x3 x4 x5])
%         % Convert g(x) to g(u)
%         syms u1 u2 u3 u4 u5 gTANAu
%         x1n = muE + u1 * sigmaE
%         x2n = muF + u2 * sigmaF
%         x3n = muA + u3 * sigmaA
%         x4n = muL + u4 * sigmaL
%         x5n = muI + u5 * sigmaI
%         gTANAu = subs(gTANA,{x1,x2,x3,x4,x5},[x1n,x2n,x3n,x4n,x5n])
%         % Solve beta minimization on g(u)
%         x0 = [0 0 0 0 0]; LB = 2*[-1;-1;-1;-1;-1]; UB = 2*[1 ;1 ;1 ;1 ;1];
%         Amat = [ ]; B = [ ]; 
%         Aeq = []; Beq = []; options = optimset('Display','off');
%         [x,bTANA,exitflag,output] = fmincon(@myfun,x0,Amat,B,Aeq,Beq,LB,UB,@mycon,options);                           
%         xE1 = subs((muE + u1 * sigmaE),{u1,u2,u3,u4,u5},x)
%         xF1 = subs((muF + u2 * sigmaF),{u1,u2,u3,u4,u5},x)
%         xA1 = subs((muA + u3 * sigmaA),{u1,u2,u3,u4,u5},x)
%         xL1 = subs((muL + u4 * sigmaL),{u1,u2,u3,u4,u5},x)
%         xI1 = subs((muI + u5 * sigmaI),{u1,u2,u3,u4,u5},x)
%         b(i) = bTANA
%         aE(i) = x(1)/b(i)
%         aF(i) = x(2)/b(i)
%         aA(i) = x(3)/b(i)
%         aL(i) = x(4)/b(i)
%         aI(i) = x(5)/b(i)
%     elseif i == 3
%         % TANA
%         % Solve Nonlinearity Index r
%         options = optimset('Display','iter','TolFun',1e-5); x0 = 1;
%         f = @(r)(LimitState(i-1) - (LimitState(i) +(1/r)*((xE(i)^(1-r)*dDdE(i)*(xE(i-1)^r - xE(i)^r)) + ...
%                                                           ((-xF(i))^(1-r)*dDdF(i)*((-xF(i-1))^r - (-xF(i))^r)) + ...
%                                                           (xA(i)^(1-r)*dDdA(i)*(xA(i-1)^r - xA(i)^r)) + ...
%                                                           (xL(i)^(1-r)*dDdL(i)*(xL(i-1)^r - xL(i)^r)) + ...
%                                                           (xI(i)^(1-r)*dDdI(i)*(xI(i-1)^r - xI(i)^r))) ) )
%         [r,fval] = fsolve(f,x0,options)
%         % Plug r into TANA expansion equation to determine g(x)
%         syms x1 x2 x3 x4 x5 gTANA
%         limitstate = sym(LimitState(i),'f')
%         xe = sym(xE(i),'f'); dgde = sym(dDdE(i),'f')
%         xf = sym(-xF(i),'f'); dgdf = sym(dDdF(i),'f')
%         xa = sym(xA(i),'f'); dgda = sym(dDdA(i),'f')
%         xl = sym(xL(i),'f'); dgdl = sym(dDdL(i),'f')
%         xi = sym(xI(i),'f'); dgdi = sym(dDdI(i),'f')
%         gTANA = symfun(limitstate + (1/r)*((xe^(1-r)*dgde*(x1^r - xe^r)) + ...
%                                        (xf^(1-r)*dgdf*(x2^r - xf^r)) + ...
%                                        (xa^(1-r)*dgda*(x3^r - xa^r)) + ...
%                                        (xl^(1-r)*dgdl*(x4^r - xl^r)) + ...
%                                        (xi^(1-r)*dgdi*(x5^r - xi^r)) ),...
%                                        [x1 x2 x3 x4 x5])
%         % Convert g(x) to g(u)
%         syms u1 u2 u3 u4 u5 gTANAu
%         x1n = muE + u1 * sigmaE
%         x2n = muF + u2 * sigmaF
%         x3n = muA + u3 * sigmaA
%         x4n = muL + u4 * sigmaL
%         x5n = muI + u5 * sigmaI
%         gTANAu = subs(gTANA,{x1,x2,x3,x4,x5},[x1n,x2n,x3n,x4n,x5n])
%         % Solve beta minimization on g(u)
%         x0 = [0 0 0 0 0]; LB = 2*[-1;-1;-1;-1;-1]; UB = 2*[1 ;1 ;1 ;1 ;1];
%         Amat = [ ]; B = [ ]; 
%         Aeq = []; Beq = []; options = optimset('Display','off');
%         [x,bTANA,exitflag,output] = fmincon(@myfun,x0,Amat,B,Aeq,Beq,LB,UB,@mycon2,options);                           
%         xE3 = subs((muE + u1 * sigmaE),{u1,u2,u3,u4,u5},x)
%         xF3 = subs((muF + u2 * sigmaF),{u1,u2,u3,u4,u5},x)
%         xA3 = subs((muA + u3 * sigmaA),{u1,u2,u3,u4,u5},x)
%         xL3 = subs((muL + u4 * sigmaL),{u1,u2,u3,u4,u5},x)
%         xI3 = subs((muI + u5 * sigmaI),{u1,u2,u3,u4,u5},x)
%         b(i) = bTANA
%         aE(i) = x(1)/b(i)
%         aF(i) = x(2)/b(i)
%         aA(i) = x(3)/b(i)
%         aL(i) = x(4)/b(i)
%         aI(i) = x(5)/b(i)
%     else
%         % Calculate alpha & beta values, and normalized design points    
%         b(i) = (LimitState(i) - ( dDdE(i)*sigmaE*uE(i) + dDdF(i)*sigmaF*uF(i) ...
%             + dDdA(i)*sigmaA*uA(i) + dDdL(i)*sigmaL*uL(i) + dDdI(i)*sigmaI*uI(i)))/ sigmahat(i)
%         aE(i) = -dDdE(i)*sigmaE / sigmahat(i)
%         aF(i) = -dDdF(i)*sigmaF / sigmahat(i)
%         aA(i) = -dDdA(i)*sigmaA / sigmahat(i)
%         aL(i) = -dDdL(i)*sigmaL / sigmahat(i)
%         aI(i) = -dDdI(i)*sigmaI / sigmahat(i)
%     end
% 
%     % Calculate percent error
%     if i == 1
%         e(i) = abs(( b(i) - b0 ) / b0)
%     else
%         e(i) = abs(( b(i) - b(i-1) ) / b(i-1))
%         if e(i) <= 0.01
%             disp('CONVERGED!')
%             break
%         end    
%     end
%     
% end
% 
% LimitState = [gx0;LimitState']
% b = [b0;b']
% e = [0;e']
% xE = [xE0;xE']
% xF = [xF0;xF']
% xA = [xA0;xA']
% xL = [xL0;xL']
% xI = [xI0;xI']
% dDdE = [dDdE0;dDdE']
% dDdF = [dDdF0;dDdF']
% dDdA = [dDdA0;dDdA']
% dDdL = [dDdL0;dDdL']
% dDdI = [dDdI0;dDdI']
% uE = [uE0;uE']
% uF = [uF0;uF']
% uA = [uA0;uA']
% uL = [uL0;uL']
% uI = [uI0;uI']
% aE = [aE0;aE']
% aF = [aF0;aF']
% aA = [aA0;aA']
% aL = [aL0;aL']
% aI = [aI0;aI']



