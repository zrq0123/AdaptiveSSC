out_X=load('E:/Matlabspace/singlecell Data/Test_9_Yan.txt');
%out_X=log2(out_X+1);
%out_X=out_X-min(min(out_X));
true_label=load('E:/Matlabspace/singlecell Data/Test_9_Yan_label.txt');
n_space=length(unique(true_label));

out_X = FilterGenesZero(out_X);

ind=knnEstimate(out_X,5);

X = matrixNormalize(out_X');

Z = admmREnSC(X,ind,0.1,0.5);

%------esitmate the cluster num-----
%K2=Estimate_Number_of_Clusters_Gap(abs(Z) + abs(Z'),10)

clusters = SpectralClustering(abs(Z) + abs(Z'), n_space);

nmi = Cal_NMI(clusters, true_label);
ari = Cal_ARI(clusters, true_label);
disp(strcat(num2str(nmi),' ',num2str(ari)));