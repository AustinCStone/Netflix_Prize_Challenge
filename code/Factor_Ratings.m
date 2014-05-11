
function [] = Factor_Ratings(Ratings,cvSet,testIdx,trR,k,iterations)
%All parameters come directly out of the provided data.


%Ratings is the matrix of all known ratings for 1978 users and 4635 movies.
%If an element of the matrix is 0, it implies that rating is unknown
%for that particular element.

%cvSet contains indices for 10-fold Cross-validation set. Row cvSet(i, :)
%contains the indices for first cross-validation set.

%testIdx are the indices in the Ratings matrix not in the cvSet on which to
%test the predictions of the algorithm

%trR is a matrix of size same as Ratings matrix but some of the entries
%have been set to be zero. This matrix is to be used as the Training
%matrix, where entries corresponding to the Test data has been set to be zero
%or unknown.

%k is the number of "features" each row of U and M should have. U is
%1978 by k and M is 4635 by k.

%Iterations is the number of iterations to solve for U and M. It should be
%high enough for U and M to approach a steady state. Roughly 30 iterations
%seems to produce good results.


%The optimal lamda (based on the graph produced by this algorithm) is .4.

%For the optimal lamda, what is the RMSE obtained for the test set?
%The optimal RMSE is about 1.07.

lambda = .05; %can be changed to start out at a different lamda variable
lambda_add = .05; %what we add to lambda on each iteration
lambda_iterations=20; %number of iterations to solve for best lambda
lambda_start = lambda; %just since we change the value of lamda below,
%so we can know what the original value of lamda was
u = size(Ratings,1);
m = size(Ratings,2);
U = rand(size(Ratings,1),k);
M = rand(size(Ratings,2),k);
lambda_iteration = 1;
calc_iteration = 1;
display = zeros(lambda_iterations,1);
%loop through all different values of lambda
while(lambda_iteration<=lambda_iterations)
    disp(lambda_iteration);
    test_index = 1;
    while(test_index<=size(cvSet,1))
        trR_temp=trR;
        trR_temp(cvSet(test_index,:))=0;
        while(calc_iteration<=iterations)
   
            %loop through and solve for all columns of M via linear regression
            m_temp = 1;
            while(m_temp<=m)
                %find indices of all users who rated this movie
                indices = find(trR_temp(:,m_temp));
                M(m_temp,:)=((U(indices,:)'*U(indices,:)+lambda.*(eye(k)))^-1)*U(indices,:)'*trR_temp(indices,m_temp);
                m_temp = m_temp+1;
            end
            
            
            %loop through and solve for all columns of U via linear regression
            u_temp = 1;
            while(u_temp<=u)
                %find all indices of all movies rated by this user 
                indices = find(trR_temp(u_temp,:));
                U(u_temp,:)= ((M(indices,:)'*M(indices,:)+lambda.*(eye(k)))^-1)*M(indices,:)'*trR_temp(u_temp,indices)';
                u_temp = u_temp+1;
            end
            
            calc_iteration = calc_iteration+1;
            
        end
        
        predicted_ratings = U*M';
        RMSE = sqrt(sum(sum((predicted_ratings(cvSet(test_index,:))-trR(cvSet(test_index,:))).^2))/length(cvSet(test_index,:)));
        display(lambda_iteration) = display(lambda_iteration)+RMSE/size(cvSet,1);
        calc_iteration = 1;
        U = rand(size(Ratings,1),k);
        M = rand(size(Ratings,2),k);
        test_index=test_index+1;
    end
    lambda=lambda+lambda_add;
    lambda_iteration = lambda_iteration+1;
    
end

%disp(display);
%disp(size(display));
%disp(size(lambda_start:lambda_add:(lambda_iterations)*lambda_add));
h = scatter(lambda_start:lambda_add:(lambda_iterations)*lambda_add,display);
xlabel('value of lambda');
ylabel('RMSE');
title('RMSE vs. Lambda');
saveas(h,'0Lamda.fig');

%find the optimal value of lambda
[val index] = min(display);
optimal_lambda = index*lamda_add;


%now that optimal lambda has been found, use it to produce a factorization
%of the ratings matrix

calc_iteration = 1;
while(calc_iteration<=iterations)
    
     %loop through and solve for all columns of M via linear regression
    m_temp = 1;
    while(m_temp<=m)
        %find indices of all users who rated this movie
        indices = find(trR(:,m_temp));
        M(m_temp,:)=((U(indices,:)'*U(indices,:)+lambda.*(eye(k)))^-1)*U(indices,:)'*trR(indices,m_temp);
        m_temp = m_temp+1;
    end
    
    
    
    u_temp = 1;
    while(u_temp<=u)
        %find all indices of all movies rated by this user 
        indices = find(trR(u_temp,:));
        U(u_temp,:)= ((M(indices,:)'*M(indices,:)+lambda.*(eye(k)))^-1)*M(indices,:)'*trR(u_temp,indices)';
        u_temp = u_temp+1;
    end
    
    calc_iteration = calc_iteration+1;
    
end
predicted_ratings = U*M';
test_set_RMSE = sqrt(sum(sum((predicted_ratings(testIdx)-Ratings(testIdx)).^2))/length(testIdx));
fprintf('Optimal Lambda is: %d, RMSE with test set is %d', optimal_lambda,test_set_RMSE);
fprintf('\n');

end


