The R scripts for "Prediction of In-hospital Death in Pneumonic COPD Exacerbation: Multicenter Cohort Study using Machine Learning Algorithm" are uploaded in this repository.  
"Main analysis" containts "External validation of BAP-65 and CURB-65", "Development of the XGBoost model", and "Model comparisons between CURB-65, BAP-65 and the XGBoost model"  

<Grid search>  
Optimal hyperparameters were searched fixing the value of eta = 0.1, subsample = 0.8, colsample_bytree = 0.8, gamma = 0, lambda = 1, and alpha = 0. Candidates of max_depth, and min_child_weight are summarized below. The area under the receiver operating characteristic curves was maximized when max_depth = 2, min_child_weight = 2, and the maximum number of trees = 39.   
   
Hyperparameters	Candidates  
Eta 0.1  
Max_depth 2, 4, 6, 8, 10  
Min_child_weight 1, 2, 3, 4, 5  
Subsample 0.8  
Colsample_bytree 0.8  
Gamma	0  
Lambda 1  
Alpha	0  
