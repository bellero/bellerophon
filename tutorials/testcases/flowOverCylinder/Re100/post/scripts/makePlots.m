clear all
close all


M=load('../results/overset.dat');
N=load('../results/single.dat');

n=size(M);
m=floor(n(1)/2)+1;

% Sampled Data
CDo = M(m:n(1),2);
CLo = M(m:n(1),3);
CDs = N(m:n(1),2);
CLs = N(m:n(1),3);

fs  = 1/M(1,1);                                % Abtastfrequenz
N   = length(CDo);                           % Laenge der Vektoren
f   = 0:1/N*fs:(N-1)/N*fs;                % Frequenzvektor

%Calculating FFTs
ftCDo = fft(CDo-mean(CDo));
ftCLo = fft(CLo-mean(CLo));
ftCDs = fft(CDs-mean(CDs));
ftCLs = fft(CLs-mean(CLs));

% PSDs
psdCDo = abs(ftCDo);
psdCDs = abs(ftCDs);
psdCLo = abs(ftCLo);
psdCLs = abs(ftCLs);

# figure(1);
# subplot(2,1,1);
# plot(f,psdCDo,'Color','red');
# hold on;
# plot(f,psdCDs,'Color','blue');
# axis([0.4,0.5,0,1])
# axis "autoy";
# hold off;
# subplot(2,1,2);
# plot(f,psdCLo,'Color','red');
# hold on;
# plot(f,psdCLs,'Color','blue');
# axis([0.1,0.3,0,1]);
# hold off;
# title('PSD of Lift')
# xlabel('f [Hz]');
# axis "autoy";

mCDo=mean(CDo);
mCDs=mean(CDs);
mCD=mCDo/mCDs*100-100;

sCDo=std(CDo);
sCDs=std(CDs);
sCD=sCDo/sCDs*100-100;

mCLo=mean(CLo);
mCLs=mean(CLs);
mCL=mCLo/mCLs*100-100;

sCLo=std(CLo);
sCLs=std(CLs);
sCL=sCLo/sCLs*100-100;
# fprintf('Average of CD:\n    %f (%f%%) \n',mCDo,mCD);
# fprintf('Standard Deviation of CD:\n    %f (%f%%) \n',sCDo,sCD);
# fprintf('Average of CL:\n    %f (single Grid: %f)\n',mCLo,mCLo-mCLs);
# fprintf('Standard Deviation of CL:\n    %f (%f%%)\n',sCLo,sCL);

dlmwrite("../results/spectra.dat",[f' psdCDo psdCDs psdCLo psdCLs]," ");

fid = fopen("../results/coeffs.tex", "w");
fputs(fid,"\\begin{tabular}{lccc}\n");
fputs(fid,"\\toprule\n");
fputs(fid," & Overset & Single & Deviation \\\\\n");
fputs(fid,"\\toprule\n");
fputs(fid,strcat("Avg. Drag & ",num2str(mCDo)," & ",num2str(mCDs)," & ",num2str(mCD)," \\% \\\\\n"));
fputs(fid,strcat("SD Drag & ",num2str(sCDo)," & ",num2str(sCDs)," & ",num2str(sCD)," \\% \\\\\n"));
fputs(fid,strcat("Avg. Lift & ",num2str(mCLo)," & ",num2str(mCLs)," & ",num2str(mCL)," \\% \\\\\n"));
fputs(fid,strcat("SD Lift & ",num2str(sCLo)," & ",num2str(sCLs)," & ",num2str(sCL)," \\% \\\\\n"));
fputs(fid,"\\toprule\n");
fputs(fid,"\\end{tabular}\n");
fclose(fid);
