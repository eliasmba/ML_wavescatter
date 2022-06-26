clear
format compact

%% Initializing
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest,YTest,anglesTest] = digitTest4DArrayData;

disp(size(XTrain))
disp(size(YTrain))
disp(size(XTest))

Y = zeros(5000,10);
Ytest = zeros(5000,10);
labels = 0:9;
for n=1:5000
    Y(n,labels(YTrain(n))+1) = 1;
    Ytest(n,labels(YTest(n))+1) = 1;
end

U = reshape(XTrain,28*28,5000).';
Utest = reshape(XTest,28*28,5000).';
Util = Utest;

%% Getting Random matrix
% Note: Need while loop check and RAM split.
NRnd = 1000;
NDisks = 150;

MeanAccDATA = zeros(1,10);
gammaDATA = zeros(1,10);
varDATA = zeros(1,10);

for ii = 1:10
    W = get_scattermat(NRnd, 28*28, NDisks, 50, 2^6);
    UnbiasedW = (W-sum(W,'all')/(size(W,1)*size(W,2)));
    NormalizedRealW = real(UnbiasedW)/sqrt(var(real(UnbiasedW),0,'all'));
    NormalizedImagW = imag(UnbiasedW)/sqrt(var(imag(UnbiasedW),0,'all'));
    VarW = var(real(UnbiasedW(:))) + 1i*var(imag(UnbiasedW(:)));
    disp(['StdDevW = ',num2str(sqrt(real(VarW))+1i*sqrt(imag(VarW)))])
    varDATA(ii) =  VarW;
    Wnor = (NormalizedRealW+1i*NormalizedImagW)/(sqrt(NRnd));
    
    %% Machine Learning
    XMat = (abs(Wnor*U.')).';
    XtilMat = (abs(Wnor*Util.')).';
    
    gammaData = linspace(0.01, 1, 100);
    accuracy = zeros(1,length(labels)); %preload
    qualityData = zeros(1,length(gammaData));
    for g = 1:length(gammaData)
        disp(['g = ',num2str(g), ' / ', num2str(length(gammaData))])
        Ytil = XtilMat*((XMat.'*XMat+gammaData(g)*eye(size(Wnor,1)))\(XMat.'*Y));
        
        [~,estIdxLabel] = max(Ytil.');
        estLabel = (estIdxLabel-1).';
        for q=1:length(labels)
            accuracy(q) = 1-sum(estLabel((1:500)+(q-1)*500) ~= labels(q))/500;
        end
        qualityData(g) = mean(accuracy);
    end
    
    [~,bestIdx] = max(qualityData);
    Ytil = XtilMat*((XMat.'*XMat+gammaData(bestIdx)*eye(size(Wnor,1)))\(XMat.'*Y));
    [~,estIdxLabel] = max(Ytil.');
    estLabel = (estIdxLabel-1).';
    for q=1:length(labels)
        accuracy(q) = 1-sum(estLabel((1:500)+(q-1)*500) ~= labels(q))/500;
        disp(['Accuracy label : ',num2str(labels(q)), ' is ', num2str(accuracy(q))]);
    end
    disp(['Best gamma : ', num2str(gammaData(bestIdx))]);
    disp(['Accuracy : ', num2str(mean(accuracy))]);
    gammaDATA(ii) =  gammaData(bestIdx);
    MeanAccDATA(ii) =  mean(accuracy);
end
disp('end');