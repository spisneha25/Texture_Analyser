clear 
display('Quantization is set to 16');

%INPUT PARAMETERS
I = input('Enter Image Location');
N = input('Enter no. of training samples');
P = input('Enter no. of training features');
W = input('Enter window size');
A = input('Offset 0, 45, 90, 135');

% SET OFFSET BASED ON ANGLE
if A == 0
    O = [0 1];
elseif A == 45
    O = [-1 1];
elseif A == 90
    O = [-1 0];
else 
    O = [-1 -1];
end

% TRAINING THE CLASSIFIER
label = fopen('Label.txt', 'wt');
imgListFile = fopen('imgs/images_read.txt', 'r');
tline = fgetl(imgListFile); 
L = cell(N, 1);
F = zeros(N, P);
sf = zeros(1, P);
currentLine = 1;
while ischar(tline)        
    splittedLine = regexp(tline, ',[ ]*', 'split');
    imagePath = fullfile('imgs', splittedLine{1});
    I = imread(imagePath);
    GLCM = graycomatrix(I,'Offset', O, 'NumLevels', 16);
    stats = Glcm_f(GLCM, 0);
    sf = [stats.autoc stats.sosvh stats.savgh stats.senth stats.maxpr]; % SELECTING TOP 5 FEATURES OBTAINED THROUGH INFORMATION GAIN FUNCTION 
    F(currentLine, :) = sf;
    L{currentLine} = splittedLine{2};
    fprintf(label, strcat(splittedLine{2}, '\n'));
    tline = fgetl(imgListFile);
    currentLine = currentLine + 1;
end;   
fclose(imgListFile);
csvwrite('Cropclassification.csv', F);
  
% CREATING THE ARFF FILE TO TRAIN THE CLASSIFIER
t = fopen('Cropclassification.csv', 'r');
l = fopen('Label.txt', 'r');
f = fopen('Cropclassification.arff', 'wt');
Train = '@relation cropclassification\n';
fprintf(f,  Train);
Train = '@attribute autocorrelation real\n';
fprintf(f,  Train);
Train = '@attribute sumofsquares real\n';
fprintf(f,  Train);
Train = '@attribute sumofaverages real\n';
fprintf(f,  Train);
Train = '@attribute sumentropy real\n';
fprintf(f,  Train);
Train = '@attribute maximumprobability real\n';
fprintf(f,  Train);
Train = '@attribute class {water,urban,hill,crop}\n';
fprintf(f,  Train);
Train = '@data\n';
fprintf(f,  Train);
cL = 1;
Train = fgets(t);
Label = fgets(l);
while ischar(Train)
    Train = strcat(Train, ',', Label, '\n');
    fprintf(f,  Train);
    Train = fgets(t);
    Label = fgets(l);
end
fclose(f);
fclose(t);
fclose(l);
 
% GLCM OF THE IMAGE WITH ADJUSTABLE WINDOW SIZE
cntr = 1;
stst = zeros(1, P);
Test = imread('imgs/dates/maghh19.tif');
[tm, tn] = size(Test);
ws = W;
i = 1;
while i <= tn
    j = 1;
    while j <= tm
        imc = imcrop(Test, [i j ws ws]);
        GLCMt = graycomatrix(imc,'Offset', O, 'NumLevels', 16);
        statstst = Glcm_f(GLCMt, 0);
        stst = [statstst.autoc statstst.sosvh statstst.savgh statstst.senth statstst.maxpr];
        exc(cntr, :) = stst;
        cntr = cntr + 1;
        j = j+ws;
    end
    i = i+ws;
end
csvwrite('Croptest.csv', exc);
 
% CREATING THE TEST ARFF FILE
t = fopen('Croptest.csv', 'r');
f = fopen('Croptest.arff', 'wt');
Test = '@relation croptest\n';
fprintf(f,  Test);
Train = '@attribute autocorrelation real\n';
fprintf(f,  Train);
Train = '@attribute sumofsquares real\n';
fprintf(f,  Train);
Train = '@attribute sumofaverages real\n';
fprintf(f,  Train);
Train = '@attribute sumentropy real\n';
fprintf(f,  Train);
Train = '@attribute maximumprobability real\n';
fprintf(f,  Train);
Train = '@attribute class {water,urban,hill,crop}\n';
fprintf(f,  Train);
Train = '@data\n';
fprintf(f,  Train);
Test = fgets(t);
while ischar(Test)
    Test = strcat(Test, ',?\n');
    fprintf(f,  Test);
    Test = fgets(t);
end
fclose(f);
fclose(t);
         
% IMPORTING WEKA IN MATLAB
javaaddpath('C:/Weka-3-6/weka.jar');
classifier = weka.classifiers.trees.NBTree;
p = java.lang.String('-t ./Cropclassification.arff -T ./Croptest.arff -p 0');
param = p.split(' ');
r = weka.classifiers.Evaluation.evaluateModel(classifier, param);
a = r.split('   ');
[am, an] = size(a);

% CLASSIFYING THE IMAGE
Test = imread('imgs/dates/maghh19.tif');
[tm, tn] = size(Test);
ws = W;
classcntr = 6;
i = 6;
inc = 6;
while i <= tn
    j = 1;
    while j <= tm
        if classcntr == 6000
            classcntr = 5999;
            inc = 5;
        else
            class = a(classcntr).split(':');
            display(class);
            for tj = i:(i+ws)
                for ti = j:(j+ws)
                    if strcmp(class(2),'water')
                        Test(ti,tj) = 0;
                    elseif strcmp(class(2),'urban')
                        Test(ti,tj) = 255;
                    elseif strcmp(class(2),'hill')
                        Test(ti,tj) = 63;
                    elseif strcmp(class(2),'crop')
                        Test(ti,tj) = 127;
                    end
                end
            end
        end
        j = j+ws;
        classcntr = classcntr + inc;
    end
    i = i+ws;
end
figure, imshow(Test)

% SUPERVISED TESTING THE CLASSIFIER
S = input('Enter no. of testing samples');
label = fopen('Label.txt', 'wt');
imgListFile = fopen('imgs/test_read.txt', 'r');
tline = fgetl(imgListFile); 
L = cell(S, 1);
F = zeros(S, P); 
sf = zeros(1, P);
currentLine = 1;
while ischar(tline)        
    splittedLine = regexp(tline, ',[ ]*', 'split');
    imagePath = fullfile('imgs', splittedLine{1});
    I = imread(imagePath);
    GLCM = graycomatrix(I,'Offset', O, 'NumLevels', 16);
    stats = Glcm_f(GLCM, 0);
    sf = [statstst.autoc statstst.sosvh statstst.savgh statstst.senth statstst.maxpr];
    F(currentLine, :) = sf;
    L{currentLine} = splittedLine{2};
    fprintf(label, strcat(splittedLine{2}, '\n'));
    tline = fgetl(imgListFile);
    currentLine = currentLine + 1;
end;   
fclose(imgListFile);
csvwrite('Croptest.csv', F);
 
% CREATING THE TEST ARFF FILE
t = fopen('Croptest.csv', 'r');
l = fopen('Label.txt', 'r');
f = fopen('Croptest.arff', 'wt');
Test = '@relation croptest\n';
fprintf(f,  Test);
Train = '@attribute autocorrelation real\n';
fprintf(f,  Train);
Train = '@attribute sumofsquares real\n';
fprintf(f,  Train);
Train = '@attribute sumofaverages real\n';
fprintf(f,  Train);
Train = '@attribute sumentropy real\n';
fprintf(f,  Train);
Train = '@attribute maximumprobability real\n';
fprintf(f,  Train);
Train = '@attribute class {water,urban,hill,crop}\n';
fprintf(f,  Train);
Train = '@data\n';
fprintf(f,  Train);
cL = 1;
Train = fgets(t);
Label = fgets(l);
while ischar(Train)
    Train = strcat(Train, ',', Label, '\n');
    fprintf(f,  Train);
    Train = fgets(t);
    Label = fgets(l);
end
fclose(f);
fclose(t);
fclose(l);

%ALL CLASSIFIER COMMANDS
% classifier = weka.classifiers.trees.J48;
% classifier = weka.classifiers.bayes.NaiveBayes;
% classifier = weka.classifiers.trees.NBTree;
% classifier = weka.classifiers.functions.MultiLayerPerceptron;
% classifier = weka.classifiers.bayes.BayesNet;
% classifier = weka.classifiers.trees.RandomForest;
% classifier = weka.classifiers.functions.Logistic;
% classifier = weka.classifiers.functions.SMO;
% classifier = weka.classifiers.functions.LibSVM;
    %javaaddpath('C:/Weka-3-6/libsvm.jar');
    %p = java.lang.String('-t ./Cropclassification.arff -T ./Croptesting.arff -p 0 -K 0 -E 0.01 -Z -seed 4');

% % CREATING THE DATA ARFF FILE TO CALCULATE INFORMATION GAIN
% t = fopen('Alldata.csv', 'r');
% l = fopen('Label.txt', 'r');
% f = fopen('Alldata.arff', 'wt');
% Train = '@relation alldata\n';
% fprintf(f,  Train);
% Train = '@attribute autocorrelation real\n';
% fprintf(f,  Train);
% Train = '@attribute contrast real\n';
% fprintf(f,  Train);
% Train = '@attribute correlation real\n';
% fprintf(f,  Train);
% Train = '@attribute clusterprominence real\n';
% fprintf(f,  Train);
% Train = '@attribute clustershade real\n';
% fprintf(f,  Train);
% Train = '@attribute dissimilarity real\n';
% fprintf(f,  Train);
% Train = '@attribute energy real\n';
% fprintf(f,  Train);
% Train = '@attribute entropy real\n';
% fprintf(f,  Train);
% Train = '@attribute homogenity real\n';
% fprintf(f,  Train);
% Train = '@attribute maximumprobability real\n';
% fprintf(f,  Train);
% Train = '@attribute sumofsquares real\n';
% fprintf(f,  Train);
% Train = '@attribute sumofaverages real\n';
% fprintf(f,  Train);
% Train = '@attribute sumofvariances real\n';
% fprintf(f,  Train);
% Train = '@attribute sumentropy real\n';
% fprintf(f,  Train);
% Train = '@attribute diffvariance real\n';
% fprintf(f,  Train);
% Train = '@attribute diffentropy real\n';
% fprintf(f,  Train);
% Train = '@attribute infogain real\n';
% fprintf(f,  Train);
% Train = '@attribute inversediff real\n';
% fprintf(f,  Train);
% Train = '@attribute inversemoment real\n';
% fprintf(f,  Train);
% Train = '@attribute class {water,urban,hill,crop}\n';
% fprintf(f,  Train);
% Train = '@data\n';
% fprintf(f,  Train);
% cL = 1;
% Train = fgets(t);
% Label = fgets(l);
% while ischar(Train)
%     Train = strcat(Train, ',', Label, '\n');
%     fprintf(f,  Train);
%     Train = fgets(t);
%     Label = fgets(l);
% end
% fclose(f);
% fclose(t);
% fclose(l);