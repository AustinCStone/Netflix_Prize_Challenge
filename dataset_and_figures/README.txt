Terms of Use:

This dataset is a subset of the Netflix prize dataset.  Unauthorized distribution of the dataset is prohibited.

Format:

On matlab console, type:
load Netflix.mat

This will load following matrices: Ratings (sparse format), trR (sparse format), trIdx, testIdx, cvSet.

Ratings is the matrix of all known ratings for 1978 users and 4635 movies. If an element of the matrix is 0, it implies that rating is unknown for that particular element. 

trR is a matrix of size same as Ratings matrix but some of the entries have been set to be zero. This matrix is to be used as the Training matrix, where entries corresponding to the Test data has been set to be zero or unknown. 

trIdx is the indices of Ratings matrix which we use as training data, i.e. trR=Ratings(trIdx).

testIdx is the indices of Ratings matrix which we have set to be zero. Once you calculate U and M, then calculate PredictedRatings=U*M^T. Now your RMSE will be
sqrt(sum(sum((PredictedRatings(testIdx)-Ratings(testIdx)).^2))/length(testIdx)), i.e. you take mean squared error of your predictions for indices corresponding to testIdx with the known ratings for testIdx stored in matrix Ratings. 

We provide you with a matrix cvSet which contains indices for 10-fold Cross-validation set. Row cvSet(i, :) contains the indices for first cross-validation set. Thus, if you want to leave first set out and perform training on 2-10 sets, then you will type:
trR1=trR;
trR1(cvSet(1,:))=0;

i.e. you set elements of your training data corresponding to the first set to be unknown. Now you train using trR1 data and after you find out the solution matrix U1 and M1, calculate RMSE using:
sqrt(sum(sum((PredictedRatings1(cvSet(1,:))-trR(cvSet(1,:))).^2))/length(cvSet(1,:))).

