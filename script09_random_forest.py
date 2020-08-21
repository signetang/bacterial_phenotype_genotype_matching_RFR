import matplotlib
matplotlib.use('Agg')
import os
import numpy as np
import pandas as pd
from math import ceil, floor
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from matplotlib.lines import Line2D
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor, GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV, cross_val_score
from sklearn.metrics import mean_squared_error, explained_variance_score
from scipy.stats import pearsonr
from script02_color_scheme import load_color_scheme
import shap

from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette, cut_tree

random_seed = 42
np.random.seed(random_seed)


def read_data(genomic_matrix, phenotype_file, phenotype_column_variables, variable_columns, training_fraction=0.6, test_fraction=0.2):
    """Reads and structures the data. Splits it into training and test data.

    genomic_matrix: The matrix of gene or Pfam presence/absence or the Pfam count matrix.
    phenotypes: The phenotype recorded for the strains.
    phenotype_column_variables: A list of the name of the phenotype and which conditions it was meassured under, for when multiple different columns contain a target value.
    variable_columns: A list of column names of columns containing the conditions the phenotype was measured under.
    training_fraction: How big a part of the data set to use for training.
    test_fraction: How big a part of the data set to use for testing.
    The remaining part of the data is used for validation (1-training_fraction-test_fraction).
    """
    val_fraction = 1.0 - training_fraction - test_fraction
    if val_fraction < 0.0:
        print('\n\nIt is required that: training_fraction + test_fraction <= 1\nThis should be fixed, else it might cause errors!\n\n')

    # Load genome and phenotype data:
    target = pd.read_csv(phenotype_file, index_col=0, sep=';')
    genomic_matrix = pd.read_csv(genomic_matrix, sep=';', index_col=0)

    # Only keep the records for strains whihc are present both in the target and in the genomic_matrix:
    strain_list = list(set([i for i in target.index.values if i in genomic_matrix.index.values]))
    if len(strain_list) == 0:
        raise Exception('No strain IDs are present both in the genomic matrix and in the phenotypic matrix.')
    target = target.loc[strain_list,:]
    genomic_matrix = genomic_matrix.loc[strain_list,:]
    #if len(target) != len(genomic_matrix):
    #    raise Exception('There is not the same number of strain genomes as phenotypes... Possibly duplicate phenotype information?')

    # Make all the phenotype values go in one column and the conditions, under which it was measured, in each their own column:
    target = rearrange_target(target, phenotype_column_variables, variable_columns)

    # Merge target values (phenotypes), conditions, and genomic info:
    data_set = target.merge(genomic_matrix, on='Strains')

    # Make sure the there is no structured order behind which strains go into training and which go into test:
    np.random.shuffle(strain_list)

    # Split the data set into test and training data:
    n_train_strains = int(ceil(training_fraction*len(strain_list)))
    n_test_strains = int(floor(test_fraction*len(strain_list)))
    n_val_strains = len(strain_list) - n_train_strains - n_test_strains

    print n_train_strains, 'training strains,', n_test_strains,'test strains, and',n_val_strains,'validation strains.'

    # List which strains go in which set:
    train_strains = strain_list[:n_train_strains]
    test_strains = strain_list[n_train_strains:n_train_strains+n_test_strains]
    if n_val_strains > 0:
        val_strains = strain_list[n_train_strains+n_test_strains:]
    else:
        val_strains = []

    # Make the separate dataframes:
    train_set = data_set.loc[train_strains,:]
    test_set = data_set.loc[test_strains,:]
    val_set = data_set.loc[val_strains,:]

    # RandomizedSearchCV will re-indexes the data, so for the CV-index-lists to be correct,
    # the data_set index must be changed but the strain order remembered!
    train_set = train_set.reset_index()
    test_set = test_set.reset_index()
    val_set = val_set.reset_index()

    return train_set, test_set, val_set


def rearrange_target(target, phenotype_column_variables, variable_columns):
    """If multiple phenotypes have been recorded and each column is the same type of measurement of a phenotype under different conditions,
    then this function makes the conditions into variables (each in their own column) and puts all the recorded phenotypes into another column.
    If additional columns in the input contain conditions, these are also kept.

    Example 1: The phenotype is the pH recorded at 20 degrees on day 2. Then the column in the original phenotype_file is called pH_20_2 and
    the data will be rearranged into 3 columns: temperature, day, and pH.
    The target name should be the first one in the underscore-separated column names of the input. The other variables in the column names must be numerical.

    Example 2: One column contains the temperature, another column contains the type of media, and 5 others contain the phenotype at different time points.
    Then the script puts all phenotype measurements into one column, and has 'temperature', 'media', and 'time' as other columns.
    """
    # Make all the phenotype values go in one column and a description of the phenotype in another column:
    target['Strains'] = target.index
    target = pd.melt(target, id_vars=['Strains']+variable_columns, var_name="Variables", value_name="Value") #the column names becomes the values of the new "Variable" coluumn

    # Split the phenotype description into different columns:
    for i in range(len(phenotype_column_variables)):
        target[phenotype_column_variables[i]] = target['Variables'].str.split('_').str[i]

    # Move the condition-columns to after the target-column (phenotype-column).
    target = target[[col for col in target.columns if col not in variable_columns] + variable_columns]

    # Make 'Strains' the index and remove unnecessary columns and rows:
    target = target.set_index('Strains')
    target = target.drop(['Variables'], axis=1)
    target = target.dropna()

    return target


def get_best_params(data_set, model, k, param_grid, n_iter, output_prefix):
    """Uses randomized grid search with k-fold cross validation to establish the best parameters from a dictionary of suggestions.
    Returns the best parameters as key word arguments (a dictionary that can be passed to the model with **).

    data_set: The data to use for testing the parameters.
    model: The model to use.
    k: The number of folds to use in k-fold cross-validation.
    params: A dictionary of parameters to test in a grid search.
    """

    # Get indexes to split the data_set into the k folds of the cross-validation
    # (data from the same strain is not split into different folds):
    CV = split_data(data_set, k)

    # Define which parameter combinations to test:
    random_grid_search = RandomizedSearchCV(estimator = model, param_distributions=param_grid, n_iter=n_iter, cv=CV, iid=False, verbose=1, return_train_score=True) #, n_jobs=2)

    # Perform the randomized grid-search to get the  best parameters:
    X,y,taxonomy,strains = data_set_X_y(data_set)
    random_grid_search.fit(X,y)

    # Update the model with the best parameters and save them in a text file:
    best_params = random_grid_search.best_params_

    # Save the best parameters in a text file:
    with open(output_prefix+'_best_parameters.txt', 'w') as f:
        f.write(str(best_params))

    # Plot the mean test and training scores for each iteration:
    plot_cross_val_params(random_grid_search.cv_results_, output_prefix)

    return best_params


def split_data(data_set, k):
    """Splits data for k-fold cross-validation such that multiple measurements from the same strain ends in the same partition.

    data_set: The data to split into k folds.
    k: The number of folds to split the data into.
    """
    folds = [] # add tupples of index-arrays for each CV iteration. Tupple example: (train_index, test_index)

    strain_list = list(set(data_set['Strains'].values))
    n = len(strain_list) #total number of strains
    max_group_size = int(ceil(n/float(k))) #max number of strains in each partition
    group_index = np.tile(np.arange(k),max_group_size)[:n] #which partition each strain belongs to
    np.random.shuffle(group_index) #randomize the group assignment

    strains = np.array(strain_list)

    # Make each k-fold:
    for i in range(k):
        test = data_set.loc[data_set['Strains'].isin(strains[group_index==i])].index.values #indexes for strains in the i'th group
        train = data_set.loc[data_set['Strains'].isin(strains[group_index!=i])].index.values #indexes for strains in the other groups
        print 'In partition number',i,'the test size is', len(test), 'and training size is', len(train)
        folds.append((train,test))
    return folds


def data_set_X_y(data_set):
    """Splits the data_set into an input matrix, X, and a target-vector, y.
    Assumes the 0th column contains the Strains (the ID), the 1st column contains the target,]
    and the 2nd column contains the name of the type of measurement.
    """
    taxonomy = data_set['Taxonomy'].values
    strains = data_set.iloc[:,0].values
    y = data_set.iloc[:,1].astype(float).values
    X = data_set.iloc[:,3:] #the 2nd column contains the name of the phenotype.. That column is dropped here.
    X = X.drop(['Taxonomy'], axis=1) #the taxonomy column comes from the genomic_matrix file, and is therefor not in one of the first columns.
    X = X.astype(float).values

    return X,y,taxonomy,strains


def plot_cross_val_params(cv_results, output_prefix):
    """Plots the """
    cv_results = pd.DataFrame(cv_results)
    cv_results_params = cv_results.filter(regex='param_').columns #the parameter values

    train_test_legend = np.array([Line2D([0], [0], markerfacecolor='blue', marker='o', lw=0), Line2D([0], [0], markerfacecolor='red', marker='^', lw=0)])

    # Plot how the mean train and test scores depend on each parameter:
    plot_number = 1
    rows = ceil(len(cv_results_params)/2.)
    plt.figure(figsize=(12,rows*5))
    for param in cv_results_params:
        plt.subplot(rows,2,plot_number)
        # Numeric parameters:
        try:
            cv_results[param].astype(float) #test if the values are numerical
            plt.scatter(cv_results['mean_train_score'],cv_results[param],color='blue')
            plt.scatter(cv_results['mean_test_score'],cv_results[param],color='red',marker='^')
        # Non-numeric parameters:
        except:
            i = 0
            non_numeric_params = {}
            y_tick_labels = []
            for non_numeric in set(cv_results[param].values):
                non_numeric_params[non_numeric] = i
                y_tick_labels.append(non_numeric)
                i += 1
            plt.scatter(cv_results['mean_train_score'],[non_numeric_params[p] for p in cv_results[param]],color='blue')
            plt.scatter(cv_results['mean_test_score'],[non_numeric_params[p] for p in cv_results[param]],color='red',marker='^')
            plt.yticks(np.arange(i),y_tick_labels)
        plt.ylabel(param.split('param_')[1])
        if plot_number == 1: #put the legend above the first plot
            plt.legend(train_test_legend,['Mean train score','Mean test score'],bbox_to_anchor=(0, 1.5))
        if plot_number > len(cv_results_params)-2: #the last two subplots
            plt.xlabel('Mean score')
        plot_number += 1
    plt.suptitle('Mean scores for the different parameters in the grid search', fontsize=20)
    plt.savefig(output_prefix+'_parameter_cross_val_grid_search.png', bbox_inches='tight')


def train_model(data_set, model):
    """Trains the model on the data_set.
    """
    X,y,taxonomy,strains = data_set_X_y(data_set)
    trained_model = model.fit(X,y)

    return trained_model


def eval_model(data_set, trained_model, colors, name, time_series=False, series='Time', color_by='Temp'):
    """Uses the trained model to predict the values for the target phenotypes of the data_set.
    Compares the predicted values to the actual values with a plot, the Pearson Correlation score, the Root Mean Square Error, and the Explained Variance.

    data_set: The data set to predict values for.
    model: The trained model to use for the prediction.
    colors: A dictionary of colors, where the key is the subpopulation.
    name: A string to use as the title of the plot.
    """

    # Get input matrix and target vector:
    X,y,taxonomy,strains = data_set_X_y(data_set)

    # Predict the target:
    y_predict = trained_model.predict(X)

    # Evaluate the prediction with root mean square error (RMSE) and Pearson correlation (PC):
    RMSE = np.sqrt(mean_squared_error(y,y_predict))
    PC = pearsonr(y,y_predict)
    EV = explained_variance_score(y,y_predict)

    # Make plot of the predicted target value, y_predict, against the actual target value, y, and add distribution plots to both axes:
    fig = sns.JointGrid(y_predict, y)
    for tax in colors:
        if len(y_predict[taxonomy==tax]) > 0:
            fig.ax_joint.scatter(y_predict[taxonomy==tax], y[taxonomy==tax], c=colors[tax], linewidth=0, edgecolor=None, alpha=0.6, label=tax) #scatter plot
            sns.kdeplot(y_predict[taxonomy==tax], ax=fig.ax_marg_x, legend=False, c=colors[tax]) #x-axis distribution
            sns.kdeplot(y[taxonomy==tax], ax=fig.ax_marg_y, vertical=True, legend=False, c=colors[tax]) #y-axis distribution
    fig.ax_joint.legend(bbox_to_anchor=(1.795, 1.35)) #loc=3)
    plt.subplots_adjust(top=0.9)
    fig.fig.suptitle(name, fontsize=16)
    fig.fig.text(0.85,0.8,'Explained Variance score='+'{:.2f}'.format(EV)+'\nPearson Correlation='+'{:.2f}'.format(PC[0])+', p-val='+'{:.2e}'.format(PC[1])+'\nRoot Mean Square Error='+'{:.2f}'.format(RMSE))
    fig.ax_joint.set_xlabel('Predicted value')
    fig.ax_joint.set_ylabel('Actual value')


    percentile90 = np.percentile(y-y_predict, 90.0)

    # Figure of the prediction errors:
    fig2 = plt.figure(figsize=(14,7))
    plt.subplot(1, 2, 1)
    for tax in colors:
        plt.hist(y[taxonomy==tax]-y_predict[taxonomy==tax], bins=30, color=colors[tax], alpha=0.4, label=tax, density=True)
    plt.xlabel('Actual value - predicted value',fontsize=14)
    plt.ylabel('Percentage',fontsize=14)
    #plt.ylabel('Number of datapoints',fontsize=14)
    plt.title('Prediction errors',fontsize=16)
    #plt.legend(loc=1)

    plt.subplot(1, 2, 2)
    for tax in colors:
        plt.hist(abs(y[taxonomy==tax]-y_predict[taxonomy==tax]), bins=30, color=colors[tax], alpha=0.4, label=tax, density=True)
    plt.axvline(percentile90, color='k', linestyle='dashed', linewidth=1, label='90th percentile')
    plt.xlabel('|Actual value - predicted value|',fontsize=14)
    plt.ylabel('Percentage',fontsize=14)
    #plt.ylabel('Number of datapoints',fontsize=14)
    plt.title('Absolute prediction errors',fontsize=16)
    plt.legend(loc=1)



    # Plot actual time series and the predicted time series (one column in the data_set must be 'Time'):
    if time_series != False:
        plot_series(data_set, y, y_predict, colors, time_series+name.replace(' ', '_'), series, color_by)

    return RMSE, PC, EV, fig, fig2, percentile90


def plot_series(data_set, y, y_predict, colors, output_name, series='Time', color_by='Temp'):
    """ Plots time series for each strain in the data_set in a different plot.
    Saves the plots by output_name and the strain name.
    Plots the actual time series, y, at the top and the predicted time series, y_predict, at the bottom.
    Colors the time series by the column color_by in the data_set. Uses the taxonomy to decide the color palette from.
    """

    data = data_set.loc[:,['Strains','Taxonomy',color_by,series]]
    data = data.set_index('Strains')
    data['y'] = y
    data['y_predict'] = y_predict

    # Plot series for each strain:
    for strain in set(data.index.values):
        fig = plt.figure(figsize=(6,6))

        # Make color-palette that depends on the taxonomy and where light colors indicate lower values of the color_by varaible:
        tax = list(set(data.loc[strain,'Taxonomy']))[0]
        color = colors[tax]
        n_colors = len(set(data.loc[strain,color_by]))
        color_palette = sns.dark_palette(color, input='rgb', n_colors=n_colors, reverse=True)

        # Plot actual values:
        plt.subplot(2,2,1)
        g = sns.scatterplot(data=data.loc[strain,:], y='y', x=series, hue=color_by, palette=color_palette)
        plt.ylim([0,1.5])
        plt.xlim([0,20])
        g.set(xticklabels=[])
        plt.ylabel('Actual values')
        plt.xlabel('')
        plt.title(strain,fontsize=14)
        plt.legend(bbox_to_anchor=(2, 1))

        # Plot predicted values:
        plt.subplot(2,2,3)
        sns.scatterplot(data=data.loc[strain,:], y='y_predict', x=series, hue=color_by, palette=color_palette)
        plt.ylim([0,1.5])
        plt.xlim([0,20])
        plt.ylabel('Predicted values')
        plt.legend(bbox_to_anchor=(2, 1))
        plt.savefig(output_name+strain+'.png')
        plt.close()

    return


def eval_permuted(data_set, trained_model, genomic_matrix, output_file):
    """Trains the model 1000 times on permuted input and compares the performance to the performance of the non-permuted model.
    Compares the predicted values to the actual values with Pearson Correlation scores and Root Mean Square Errors.
    """
    # Get input matrix and target vector:
    X,y,taxonomy,strains = data_set_X_y(data_set)

    # Predict the target:
    y_predict = trained_model.predict(X)

    # Evaluate the prediction with root mean square error (RMSE) and Pearson correlation scores:
    RMSE = np.sqrt(mean_squared_error(y,y_predict))
    PC = pearsonr(y,y_predict)[0]
    EV = explained_variance_score(y,y_predict)

    # Get the length of the list of all genomic features:
    genomic_matrix = pd.read_csv(genomic_matrix, sep=';', index_col=0)
    genomic_matrix = genomic_matrix.drop('Taxonomy',axis=1) #already have the taxonomic information in an array
    n_genomic_features = len(genomic_matrix.columns.values)

    # New matrix which will contain the permuted genomic content:
    permuted_input = pd.DataFrame(index=strains, data=X[:,:X.shape[1]-n_genomic_features]) #the non-genomic data
    permuted_input = permuted_input.rename_axis('Strains')

    # Save results to:
    f = open(output_file, 'w')

    # lists to keep the scores in:
    RMSE_list = []
    PC_list = []
    EV_list = []

    # 1000 times: Permute the profiles of each Pfam domain and see how well the model predicts on that:
    for i in range(1000):
        # Each strain has the same genomic content for every instance after the shuffeling:
        input = permuted_input.merge(genomic_matrix.apply(np.random.permutation, axis=0), on='Strains')
        y_predict_randomized = trained_model.predict(input)
        RMSE_list.append(np.sqrt(mean_squared_error(y,y_predict_randomized)))
        PC_list.append(pearsonr(y,y_predict_randomized)[0])
        EV_list.append(explained_variance_score(y,y_predict_randomized))

    f.write('0-hypothesis: "The model predicts equally well when the genetic values are randomly permuted for each feature."\n')
    f.write('Alternative hypothesis: "The model predicts better than when the genetic values are randomly permuted for each feature."\n')
    f.write('If the Pearson Correlation of '+str(PC)+' is above the interval ['+str(np.mean(PC_list)-2*np.std(PC_list))+', '+str(np.mean(PC_list)+2*np.std(PC_list))+'], then the 0-hypothesis can be rejected on a 2.5% significance level.\n')
    f.write('If the Root Mean Square Error of '+str(RMSE)+' is below the interval ['+str(np.mean(RMSE_list)-2*np.std(RMSE_list))+', '+str(np.mean(RMSE_list)+2*np.std(RMSE_list))+'], then the 0-hypothesis can be rejected on a 2.5% significance level.\n')
    f.write('If the Explained Variance score of '+str(EV)+' is above the interval ['+str(np.mean(EV_list)-2*np.std(EV_list))+', '+str(np.mean(EV_list)+2*np.std(EV_list))+'],  then the 0-hypothesis can be rejeced on a 2.5% significance level.\n')

    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.hist(RMSE_list, bins=50)
    plt.title('Root Mean Square Error')
    plt.xlabel('Non-permuted RMSE: '+'{:.3f}'.format(RMSE))
    plt.subplot(1,3,2)
    plt.hist(PC_list, bins=50)
    plt.title('Pearson Correlation')
    plt.xlabel('Non-permuted PC: '+'{:.3f}'.format(PC))
    plt.subplot(1,3,3)
    plt.hist(EV_list, bins=50)
    plt.title('Explained Variance')
    plt.xlabel('Non-permuted EV: '+'{:.3f}'.format(EV))
    plt.subplots_adjust(top=0.8)
    plt.suptitle('Scores for the model when predicting on input where the\ngenetic values are randomly permuted for each feature', fontsize=16)
    plt.savefig(output_file+'eval_permuted_histograms.png') # bbox_inches='tight')
    plt.close()

    # lists to keep the scores in:
    RMSE_list = []
    PC_list = []
    EV_list = []

    # 1000 times: Switch the genomic content between the strains (each strain gets another strains' genomic profile) and see how well the model predicts on that:
    for i in range(1000):
        # Each strain has the same genomic content for every instance after the shuffeling:
        index = genomic_matrix.index
        input = genomic_matrix.sample(frac=1)
        input.index = index
        input = permuted_input.merge(input, on='Strains')

        y_predict_randomized = trained_model.predict(input)
        RMSE_list.append(np.sqrt(mean_squared_error(y,y_predict_randomized)))
        PC_list.append(pearsonr(y,y_predict_randomized)[0])
        EV_list.append(explained_variance_score(y,y_predict_randomized))

    f.write('0-hypothesis: "The model predicts equally well when the genetic values are randomly switched between strains."\n')
    f.write('Alternative hypothesis: "The model predicts better than when the genetic values are randomly switched between strains."\n')
    f.write('If the Pearson Correlation of '+str(PC)+' is above the interval ['+str(np.mean(PC_list)-2*np.std(PC_list))+', '+str(np.mean(PC_list)+2*np.std(PC_list))+'], then the 0-hypothesis can be rejected on a 2.5% significance level.\n')
    f.write('If the Root Mean Square Error of '+str(RMSE)+' is below the interval ['+str(np.mean(RMSE_list)-2*np.std(RMSE_list))+', '+str(np.mean(RMSE_list)+2*np.std(RMSE_list))+'], then the 0-hypothesis can be rejected on a 2.5% significance level.\n')
    f.write('If the Explained Variance score of '+str(EV)+' is above the interval ['+str(np.mean(EV_list)-2*np.std(EV_list))+', '+str(np.mean(EV_list)+2*np.std(EV_list))+'],  then the 0-hypothesis can be rejeced on a 2.5% significance level.\n')

    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.hist(RMSE_list, bins=50)
    plt.title('Root Mean Square Error')
    plt.xlabel('Non-permuted RMSE: '+'{:.3f}'.format(RMSE))
    plt.subplot(1,3,2)
    plt.hist(PC_list, bins=50)
    plt.title('Pearson Correlation')
    plt.xlabel('Non-permuted PC: '+'{:.3f}'.format(PC))
    plt.subplot(1,3,3)
    plt.hist(EV_list, bins=50)
    plt.title('Explained Variance')
    plt.xlabel('Non-permuted EV: '+'{:.3f}'.format(EV))
    plt.subplots_adjust(top=0.8)
    plt.suptitle('Scores for the model when predicting on input where the\ngenetic values are switched between strains', fontsize=16)
    plt.savefig(output_file+'eval_switched_histograms.png') # bbox_inches='tight')
    plt.close()

    f.close()


def feature_importance(data_set, trained_model, output_prefix, non_genetic_features, name_dictionary=None, method='impurity_decrease', score='EV', n_features=30, iterations_per_feature=10):
    """Calculates the n_features most important features based on either impurity_decrease or permutation strategy.
    Saves a figure of the feature importances with standard deviations.
    Can either use the featue names directly or add more descrption by use of a dictionary."""
    output_prefix += '_feature_importances'

    variables = data_set.drop('Taxonomy', axis=1).columns.values[3:] #not the 'Strain' column or the 2nd column, that contains the name of the phenotype or the target value column.

    if name_dictionary != None:
        if output_prefix.split('/')[-1][:4] == 'gene':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            #variables = np.array([','.join([v.split('_')[1] if v.split('_')[1] != 'GCF' else v for v in var.split(',')])+': '+','.join([name_dictionary[v] for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
            variables = np.array([', '.join([name_dictionary[v]+' ('+v.split('_')[1]+')' if v.split('_')[1] != 'GCF' else name_dictionary[v]+' ('+v+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
        elif output_prefix.split('/')[-1][:4] == 'pfam':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            variables = np.array([', '.join([v+' ('+name_dictionary[v]+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])

    if method == 'impurity_decrease':
        # Next three lines are from http://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html
        importances = trained_model.feature_importances_
        importance_std = np.std([tree.feature_importances_ for tree in trained_model.estimators_], axis=0)
        index = np.argsort(importances)[::-1]
        output_prefix += '_impurity_decrease'

    elif method == 'permutation':
        # Get input matrix and target vector:
        X,y,taxonomy,strains = data_set_X_y(data_set)

        importances, importance_std, index = permutationStrategy(X, y, trained_model, len(variables), score, iterations_per_feature)
        output_prefix += '_permutation_'+score

    #outfile = open(output_prefix+'.txt', 'w')
    #outfile.close()

    # Plot the n most important features and their relative importance:
    fig = plt.figure(figsize=(15,15))
    plt.bar(range(n_features), importances[index][:n_features], yerr=importance_std[index][:n_features], align="center")
    fig.subplots_adjust(bottom=0.60)
    plt.xticks(range(n_features), variables[index][:n_features], rotation='vertical')
    plt.xlim([-1, n_features])
    plt.title('Feature importances')
    fig.savefig(output_prefix+'.png')
    plt.close(fig)

    return variables[index][:n_features]


def iterativePermutation(train_data, test_data, model, n_iterations, non_genetic_features, outputfile, name_dictionary=None, features_to_remove='default', score='EV', iterations_per_feature=10):
    """Iteratively remove features with lowest importance and build a new model based on the smaller set of features.
    n_iterations is the number of iterations to perform.
    features_to_remove is the number of features to remove per iteration. The default is calculated from n_iterations and the total number of features.
    iterationsations_per_feature is the number of times to re-run the permutationStrategy.
    """
    variables = train_data.drop('Taxonomy', axis=1).columns.values[3:] #not the 'Strain' column or the 2nd column, that contains the name of the phenotype or the target value column.

    if name_dictionary != None:
        if output_prefix.split('/')[-1][:4] == 'gene':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            #variables = np.array([','.join([v.split('_')[1] if v.split('_')[1] != 'GCF' else v for v in var.split(',')])+': '+','.join([name_dictionary[v] for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
            variables = np.array([', '.join([name_dictionary[v]+' ('+v.split('_')[1]+')' if v.split('_')[1] != 'GCF' else name_dictionary[v]+' ('+v+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
        elif output_prefix.split('/')[-1][:4] == 'pfam':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            variables = np.array([', '.join([v+' ('+name_dictionary[v]+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])

    if features_to_remove == 'default':
        features_to_remove = int(floor(len(variables)/(n_iterations+1.0))) #number of features to remove in each iteration
        print 'Features removed per iteration:', features_to_remove
    # else use whatever value is given...

    X_train,y_train,taxonomy_train,strains_train = data_set_X_y(train_data)
    X_test,y_test,taxonomy_test,strains_test = data_set_X_y(test_data)

    trained_model = model.fit(X_train, y_train)

    old_PC = 0

    f = open(outputfile+'.txt','w')

    for i in range(n_iterations):
        # Get the importances:
        importances, importance_std, index = permutationStrategy(X_train, y_train, trained_model, len(variables), score, iterations_per_feature)

        # Remove the least important features and retrain the model:
        X_train = X_train[:,index[:-features_to_remove]]
        trained_model.fit(X_train, y_train)

        # Predict the test set with the re-trained model (now using fewer features):
        X_test = X_test[:,index[:-features_to_remove]]
        y_predict = trained_model.predict(X_test)

        # Evaluate the prediction with root mean square error (RMSE), Pearson correlation (PC), and Explained Variance (EV):
        RMSE = np.sqrt(mean_squared_error(y_test,y_predict))
        PC = pearsonr(y_test,y_predict)[0]
        EV = explained_variance_score(y_test,y_predict)

        print 'iteration:',i,', Scores (RMSE, PC, EV):', RMSE, PC, EV
        f.write('iteration: '+str(i)+', Scores (RMSE, PC, EV): '+str(RMSE)+', '+str(PC)+', '+str(EV)+'\n')

        # Remove the least important features:
        variables = variables[index[:-features_to_remove]]

        if PC >= old_PC:
            best_features = [i, RMSE, PC, EV, variables, importances[index[:-features_to_remove]], importance_std[index[:-features_to_remove]]]
            old_PC = PC

        f.write(', '.join(variables[:20])+'\n')


    f.write(', '.join(variables)+'\n')
    f.write(', '.join([str(im) for im in importances[index[:-features_to_remove]]])+'\n')
    f.write(', '.join([str(im) for im in importance_std[index[:-features_to_remove]]])+'\n')
    f.write('BEST:\nIteration: '+str(best_features[0])+', Scores (RMSE, PC, EV):'+str(best_features[1])+', '+str(best_features[2])+', '+str(best_features[3])+'\n')
    f.write('Features:\n'+', '.join(best_features[4])+'\n')
    f.write('Importances:\n'+', '.join([str(im) for im in best_features[5]])+'\n')
    f.write('Std:\n'+', '.join([str(im) for im in best_features[6]]))
    f.close()

    print variables, importances[index[:-features_to_remove]], importance_std[index[:-features_to_remove]]
    print 'Best:', best_features


def permutationStrategy(X, y, trained_model, n_variables, score='EV', iterations_per_feature=100):

    # Model prediction:
    y_predict = trained_model.predict(X)

    # Evaluate the prediction with root mean square error (RMSE), Pearson correlation (PC), and Explained Variance (EV):
    RMSE = np.sqrt(mean_squared_error(y,y_predict))
    PC = pearsonr(y,y_predict)[0]
    EV = explained_variance_score(y,y_predict)

    # lists to hold the feature importances:
    importances = np.zeros(n_variables)
    importance_std = np.zeros(n_variables)

    # Permute each feature 100 times, get the average of the model predicitons and calculate the drop/rise in scores from the 'real' score:
    for feature in range(n_variables):

        # List to keep the scores of all the permutations of a feature:
        importances_feature = np.zeros(iterations_per_feature)
        for i in range(len(importances_feature)):
            # Permute feature:
            X_permuted = X.copy()
            X_permuted[:,feature] = np.random.permutation(X_permuted[:,feature])

            # Predict the target for permuted input:
            y_predict = trained_model.predict(X_permuted)

            # Evaluate predicitons:
            if score == 'RMSE':
                importances_feature[i] = np.sqrt(mean_squared_error(y,y_predict))
            elif score == 'PC':
                importances_feature[i] = pearsonr(y,y_predict)[0]
            elif score == 'EV':
                importances_feature[i] = explained_variance_score(y,y_predict)

        # The feature importances (change in prediction score):
        if score == 'RMSE':
            importances[feature] = np.mean(importances_feature) - RMSE #lower scores expected with the real input
        elif score == 'PC':
            importances[feature] = PC - np.mean(importances_feature) #higher score expected with the real input
        elif score == 'EV':
            importances[feature] = EV - np.mean(importances_feature) #higher scores expected with the real input

        # Standard deviation of prediction scores for permuted input:
        importance_std[feature] = np.std(importances_feature)

    index = np.argsort(importances)[::-1]

    return importances, importance_std, index


def shapImportances(test_data, trained_model, output_prefix, non_genetic_features, name_dictionary=None, n_features=30, individual_sample_plots=False):
    """Get the most important features, as identified by SHAP and make SHAP feature importance plots.
    """
    variables = test_data.drop('Taxonomy', axis=1).columns.values[3:] #not the 'Strain' column or the 2nd column, that contains the name of the phenotype or the target value column.

    if name_dictionary != None:
        if output_prefix.split('/')[-1][:4] == 'gene':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            #variables = np.array([','.join([v.split('_')[1] if v.split('_')[1] != 'GCF' else v for v in var.split(',')])+': '+','.join([name_dictionary[v] for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
            variables = np.array([', '.join([name_dictionary[v]+' ('+v.split('_')[1]+')' if v.split('_')[1] != 'GCF' else name_dictionary[v]+' ('+v+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])
        elif output_prefix.split('/')[-1][:4] == 'pfam':
            # Use a dictionary to convert the genetic IDs to more 'readable' names:
            variables = np.array([', '.join([v+' ('+name_dictionary[v]+')' for v in var.split(',')]) if var not in non_genetic_features else var for var in variables])

    X,y,taxonomy,strains = data_set_X_y(test_data)
    X = pd.DataFrame(X)
    X.columns = variables

    # Explain the model's predictions with SHAP values:
    explainer = shap.TreeExplainer(trained_model)
    shap_values = explainer.shap_values(X)

    # Plot the effects of the features:
    shap.summary_plot(shap_values, X, max_display=n_features)
    plt.savefig(output_prefix+'_Shap_feature_importances.png', bbox_inches='tight')
    plt.close()

    if individual_sample_plots:
        for i in range(len(X)):
            shap.force_plot(explainer.expected_value, shap_values[i,:], X.iloc[i,:], matplotlib=True)
            plt.xlabel(strains[i]+', '+taxonomy[i]+', '+', '.join([str(v) for v in X.iloc[i,:4]]))
            plt.savefig(output_prefix+'_Shap_feature_importances_'+strains[i]+'_'+'_'.join([str(v) for v in X.iloc[i,:4]])+'.png', bbox_inches='tight')
            plt.close()


def pfam_dict(path):
    if path[-1] != '/':
        path += '/'

    pfam_dict = {}

    for file in os.listdir(path):
        strain_pfams = pd.read_csv(path+file, sep='\s+', skiprows=28, header=None, dtype='str', usecols=[5,6])
        #strain_pfams = strain_pfams.drop_duplicates() #faster to use this line?
        for i in range(len(strain_pfams)):
            pfam_dict[strain_pfams.iloc[i,0]] = strain_pfams.iloc[i,1]

    return pfam_dict


def gene_dict(pa_file):
    gene_dict = {}

    gene_names = pd.read_csv(pa_file, sep=',', dtype='str', usecols=[0,1,2])
    gene_names = gene_names.fillna("")
    for i in range(len(gene_names)):
        if gene_names.iloc[i,1] != "" and len(gene_names.iloc[i,1]) < 6:
            gene_dict[gene_names.iloc[i,0]] = gene_names.iloc[i,1] # short name
        else:
            gene_dict[gene_names.iloc[i,0]] = gene_names.iloc[i,2] # description

    return gene_dict


def necessary_number_of_training_strains(train_data, test_data, model, colors, output_prefix):
    # For printing the feature names (of most important features):
    variables = train_data.drop('Taxonomy', axis=1).columns.values[3:] #not the 'Strain' column or the 2nd column, that contains the name of the phenotype or the target value column.

    # Reading the data:
    X_train,y_train,taxonomy_train,strains_train = data_set_X_y(train_data)
    X_test,y_test,taxonomy_test,strains_test = data_set_X_y(test_data)

    # Output files:
    evaluation_scores_file = open(output_prefix+'_necessary_number_of_training_strains_Evaluation_scores.txt','w')
    evaluation_scores_file.write('Explained Variance\tPearson Correlation\tp-val\tRoot Mean Square Error\n')
    feature_importances_file = open(output_prefix+'_necessary_number_of_training_strains_Feature_importances.txt','w')

    for n_strains in range(10, len(set(strains_train)), 5):
        print 'Number of strains included in the training: '+str(n_strains)+'\n'
        for subset in range(5): # get 5 different data subsets
            strain_subset = np.random.choice(list(set(strains_train)), n_strains, replace=False) # select randomly n_strains strains as a subset of the training data.
            subset_index = np.array([i for i in range(len(strains_train)) if strains_train[i] in strain_subset])

            # train model:
            trained_model = model.fit(X_train[subset_index], y_train[subset_index])

            # feature importance:
            importances = trained_model.feature_importances_
            importance_index = np.argsort(importances)[::-1]
            print variables[importance_index][:30], '\n'
            feature_importances_file.write('\t'.join(variables[importance_index][:30])+'\n')

            # test model on separate test data:
            y_predict = trained_model.predict(X_test)

            # Evaluate the prediction with root mean square error (RMSE) and Pearson correlation (PC):
            RMSE = np.sqrt(mean_squared_error(y_test,y_predict))
            PC = pearsonr(y_test,y_predict)
            EV = explained_variance_score(y_test,y_predict)
            scores = '{:.2f}'.format(EV)+'\t'+'{:.2f}'.format(PC[0])+'\t'+'{:.2e}'.format(PC[1])+'\t'+'{:.2f}'.format(RMSE)+'\n'
            print scores
            evaluation_scores_file.write(scores)

            # For one of the models with a subset of n_strains, plot the predictions:
            # Make plot of the predicted target value, y_predict, against the actual target value, y, and add distribution plots to both axes:
            fig = sns.JointGrid(y_predict, y_test)
            for tax in colors:
                if len(y_predict[taxonomy_test==tax]) > 0:
                    fig.ax_joint.scatter(y_predict[taxonomy_test==tax], y_test[taxonomy_test==tax], c=colors[tax], linewidth=0, edgecolor=None, alpha=0.6, label=tax) #scatter plot
                    sns.kdeplot(y_predict[taxonomy_test==tax], ax=fig.ax_marg_x, legend=False, c=colors[tax]) #x-axis distribution
                    sns.kdeplot(y_test[taxonomy_test==tax], ax=fig.ax_marg_y, vertical=True, legend=False, c=colors[tax]) #y-axis distribution
            fig.ax_joint.legend(bbox_to_anchor=(1.795, 1.35)) #loc=3)
            plt.subplots_adjust(top=0.9)
            fig.fig.suptitle(str(n_strains)+' strains used for training, subset no. '+str(subset), fontsize=16)
            fig.fig.text(0.85,0.8,'Explained Variance score='+'{:.2f}'.format(EV)+'\nPearson Correlation='+'{:.2f}'.format(PC[0])+', p-val='+'{:.2e}'.format(PC[1])+'\nRoot Mean Square Error='+'{:.2f}'.format(RMSE))
            fig.ax_joint.set_xlabel('Predicted value')
            fig.ax_joint.set_ylabel('Actual value')
            fig.savefig(output_prefix+'_'+str(n_strains)+'strains_subset'+str(subset)+'.png', bbox_inches='tight')
            plt.close('all')


    evaluation_scores_file.close()
    feature_importances_file.close()

    df = pd.read_csv(output_prefix+'_necessary_number_of_training_strains_Evaluation_scores.txt', sep='\t')
    n_train_strains = np.repeat(np.arange(10,len(set(strains_train)),5),5)

    plt.figure()
    ax = plt.subplot(111)
    plt.scatter(n_train_strains,df['Pearson Correlation'], color='deepskyblue', label='Pearson Correlation', marker='*')
    plt.scatter(n_train_strains,df['Explained Variance'], color='mediumvioletred', label='Explained Variance', marker='d')
    plt.scatter(n_train_strains,df['Root Mean Square Error'], color='mediumslateblue', label='Root Mean Square Error')
    plt.xlim([0,len(set(strains_train))])
    plt.ylim([0,1])
    plt.xlabel('Number of strains in the training set',fontsize=13)
    plt.ylabel('Evaluation score',fontsize=13)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.16), ncol=3)
    plt.suptitle('Gene presence/absence model', fontsize=20)
    plt.title('Evaluation of test set predictions for increasing size of training set',fontsize=14)
    plt.savefig(output_prefix+'_necessary_number_of_training_strains_Evaluation.png', bbox_inches='tight')
    plt.close()




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomic_matrix', help='The gene or Pfam presence/absence matrix or the Pfam count matrix. First column must be the "Strains" ID.')
    parser.add_argument('--phenotype_file', help='The phenotype data. First column must be the "Strains" ID.')
    parser.add_argument('--colors', help='The color dictionary to use.')
    parser.add_argument('--output_prefix', help='A prefix to give the output (can be just the path).')
    parser.add_argument('--phenotype_column_variables', help="List of the conditions that the phenotype was collected under, when multiple different columns contain target values. First item in the list must be the phenotype name. Example: 'Measurement,Volume,Temperature,Salt,Yeast', or for a time series: 'Time,Value'.")
    parser.add_argument('--variable_columns', help="List of the columns that contain conditions, which the phenotype was measured under. Example: 'Temperature,Media'")
    args = parser.parse_args()

    # Load the color dictionary:
    colors = load_color_scheme(args.colors)

    # Read and structure the data, and split it into a trainig set and a test set:
    phenotype_column_variables = args.phenotype_column_variables.split(',') if (len(args.phenotype_column_variables) > 0) else []
    variable_columns = args.variable_columns.split(',') if (len(args.variable_columns) > 0) else []
    training_fraction = 0.75
    test_fraction = 0.25
    train_set, test_set, val_set = read_data(args.genomic_matrix, args.phenotype_file, phenotype_column_variables, variable_columns, training_fraction=training_fraction, test_fraction=test_fraction)
    non_genetic_features = phenotype_column_variables+variable_columns


    # Define which model to use:
    model = RandomForestRegressor()

    # Get the best parameters for the model by randomized grid-search with k-fold cross-validation:
    # Parameter options to test:
    n_estimators = [10]+list(np.arange(50, 800, 50))
    max_features = ['auto'] #, 'sqrt', 'log2']
    max_depth = list(np.arange(10, 160, 10))
    max_depth.append(None)
    min_samples_split = list(np.arange(2, 17, 2))
    min_samples_leaf = list(np.arange(1, 10, 1))
    oob_score = [True, False]
    param_grid = {'random_state': [random_seed],\
                  'n_estimators': n_estimators,\
                  'max_features': max_features,\
                  'max_depth': max_depth,\
                  'min_samples_split': min_samples_split,\
                  'min_samples_leaf': min_samples_leaf,\
                  'oob_score': oob_score  } #,\ 'criterion': criterion  }

    # Update model with the best parameters based on randomized grid-search with k-fold cross-validation:
    n_iter = 40 #iterations of the randomized grid-search
    k = 3 #k-fold cross-validation

    print '\nSearching for the best parameters using randomized grid-search with '+str(k)+'-fold cross-validation...'
    best_param_model = model.set_params(**get_best_params(train_set, model, k, param_grid, n_iter, args.output_prefix))

    #For testing: Naive model/my best guess:
    #best_param_model = model.set_params(**{'n_estimators':10, 'min_samples_leaf':6}) #faster than getting the best params. Only for tests.
    test_model = False

    print best_param_model.get_params()

    # Train the model on the training data:
    print '\nTraining the model using the best parameters...'
    trained_model = train_model(train_set, best_param_model)

    # Evaluate the model on training data and on test data:
    print '\nEvaluating the model on the training data...'
    RMSE_train, PC_train, EV_train, figure_train, error_figure_train, percent90_train = eval_model(train_set, trained_model, colors, name='Training data') #, time_series=args.output_prefix)
    print '\nTraining Root Mean Squared Error:',RMSE_train,'\nTraining Pearson Correlation scores:',PC_train,'\nTraining Explained Variance score:',EV_train
    print '\n90th percentile of the absolute error: ',percent90_train
    figure_train.savefig(args.output_prefix+'_train_figure.png', bbox_inches='tight')
    plt.close('all')

    print '\nEvaluating the model on the test data...'
    RMSE_test, PC_test, EV_test, figure_test, error_figure_test, percent90_test = eval_model(test_set, trained_model, colors, name='Test data') #, time_series=args.output_prefix)
    print '\nTest Root Mean Squared Error:',RMSE_test,'\nTest Pearson Correlation scores:',PC_test,'\nTest Explained Variance score:',EV_test
    print '\n90th percentile of the absolute error: ',percent90_test
    figure_test.savefig(args.output_prefix+'_test_figure.png', bbox_inches='tight')
    error_figure_test.savefig(args.output_prefix+'_test_errors.png', bbox_inches='tight')
    plt.close('all')

    f = open(args.output_prefix+'_Train_and_test_scores.txt','w')
    f.write('Training Root Mean Squared Error: '+str(RMSE_train)+'\nTraining Pearson Correlation scores: '+','.join([str(p) for p in PC_train])+'\nTraining Explained Variance score:'+str(EV_train)+'\n')
    f.write('Test Root Mean Squared Error: '+str(RMSE_test)+'\nTest Pearson Correlation scores: '+','.join([str(p) for p in PC_test])+'\nTest Explained Variance score:'+str(EV_train))
    f.close()

    print '\nEvaluating the model on permuted test data...'
    eval_permuted(test_set, trained_model, args.genomic_matrix, output_file=args.output_prefix+'_model_evaluations_of_permuted_inputs.txt')

    if args.genomic_matrix[-25:] in ['pruned_pfam_pa_matrix.txt', 'ned_pfam_count_matrix.txt']:
        print 'Getting the names of the Pfam IDs...'
        name_dictionary = pfam_dict('../data/Pfam/')
    elif args.genomic_matrix[-25:] in ['pruned_gene_pa_matrix.txt','ned_gene_count_matrix.txt']:
        print 'Getting the names of the gene IDs...'
        name_dictionary = gene_dict('../data/roary/new_gene_presence_absence.csv')

    print '\nFinding the most important features with Shap...'
    shapImportances(test_set, trained_model, args.output_prefix, non_genetic_features, individual_sample_plots=False, name_dictionary=name_dictionary, n_features=40)
    #print '\nFinding the most important features with impurity decrease...'
    #feature_importance(test_set, trained_model, args.output_prefix, name_dictionary=name_dictionary, non_genetic_features=non_genetic_features, method='impurity_decrease')
    #n_iterations, features_to_remove, iterations_per_feature, score = len(train_set.columns)-5, 1, 10, 'PC'
    #print '\nFinding the most important features with permutation strategy...'
    #feature_importance(test_set, trained_model, args.output_prefix, name_dictionary=name_dictionary, non_genetic_features=non_genetic_features, method='permutation', iterations_per_feature=iterations_per_feature)

    # Save the parameters that were used in a text file:
    with open(args.output_prefix+'_settings.txt', 'w') as f:
        f.write('\nOutput prefix = '+args.output_prefix)
        f.write('\nGenomic matrix = '+args.genomic_matrix)
        f.write('\nPhenotype file = '+args.phenotype_file)
        f.write('\nNon genetic features = ['+','.join(non_genetic_features)+']')
        f.write('\ntraining_fraction = '+str(training_fraction))
        f.write('\ntest_fraction = '+str(test_fraction))
        f.write('\nResulting validation fraction = '+str(1-training_fraction-test_fraction))
        f.write('\nk-fold cross-validation was used = '+str(not test_model))
        if not test_model:
            f.write('\n\tk = '+str(k))
            f.write('\n\tn_iter = '+str(n_iter))
            f.write('\n\tgrid = '+str(param_grid))
            f.write('\n\tOptimal parameters used = '+str(best_param_model.get_params()))
        else:
            f.write('\nTest parameters used = '+str(best_param_model.get_params()))
        f.write('\nFor permutation strategy feature importance analyses:')
        f.write('\n\tn_iterations = '+str(n_iterations))
        f.write('\n\tfeatures_to_remove = '+str(features_to_remove))
        f.write('\n\titerations_per_feature = '+str(iterations_per_feature))
        f.write('\n\tscore = '+score)
        f.write('\nVersion control = '+str(version))


    print 'Starting analysis on the necessary number of training strains...'
    necessary_number_of_training_strains(train_set, test_set, best_param_model, colors, args.output_prefix)
