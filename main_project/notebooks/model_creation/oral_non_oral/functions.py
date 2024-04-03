import pandas as pd
from sklearn.preprocessing import StandardScaler
from time import perf_counter, strftime, gmtime
import numpy as np

from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV



import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns

from time import perf_counter, gmtime, strftime

from sklearn.linear_model import LogisticRegression
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC

from sklearn.model_selection import RepeatedKFold, train_test_split, cross_val_score, GridSearchCV, RepeatedStratifiedKFold


from sklearn.metrics import roc_curve, roc_auc_score, confusion_matrix, ConfusionMatrixDisplay, make_scorer




def remove_high_corr(df: pd.DataFrame, threshold: float):
  """
    Remove highly correlated variables from a DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame containing variables to be checked for correlation.
    - threshold (float): The correlation threshold above which variables will be removed.

    Returns:
    - pd.DataFrame: A DataFrame with highly correlated variables removed.

    This function calculates the correlation matrix of the input DataFrame and removes
    variables that have a correlation coefficient greater than the specified threshold.
    """

  corr_matrix = df.corr().abs()

  # selecting upper triangle
  matrix_tri_sup = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k = 1).astype(bool))

  # selecting variables to be deleted
  remove = []

  for column in matrix_tri_sup.columns:
    if any(matrix_tri_sup[column] > threshold):
      remove.append(column)  
  
  print(f'Number of excluded variables: {len(remove)}')

  return df.drop(remove, axis = 1)


def scale_variables(df: pd.DataFrame):
  scaler = StandardScaler()
  scaled_data_array = scaler.fit_transform(df)

  # Merging the column name with the scale data array
  scaled_df = pd.DataFrame(scaled_data_array, columns=df.columns)

  return scaled_df

def best_params_grid(x, y, model_params, n_splits=5, n_repeats=5, scoring='roc_auc'):
    '''DataFrame, DataFrame, Dictionary, int, int --> DataFrame

    -----------------------------------------------------------------------------

    This function recieves the X and Y dataframes. It also recieves a dictionary, 
    with the model and parameters to be testes. As the function uses RepeatedStratifiedKFold,
    it can also recieve the n_splits and n_repeats parameters.
    Tha function will execute GridSearchCV on all parameters for the model.
    It returns two datasets, the first one is formatted to contain only the relevant information.
    The second one constains all information about the models, in case it's necessary'''

    t0 = perf_counter()
    np.random.seed(1428)

    cv = RepeatedStratifiedKFold(n_splits = n_splits, n_repeats = n_repeats)
    scores = []
    full_results = {}

    for model_name, mp in model_params.items():
        now = strftime("%H:%M", gmtime())
        print(f"Starting Grid Search for {model_name}: {now}")
        clf =  GridSearchCV(mp['model'], mp['params'], cv=cv,
                                return_train_score=True, scoring=scoring)
        clf.fit(x, y)
        scores.append({
            'model': model_name,
            'best_score': clf.best_score_,
            'best_params': clf.best_params_,
        })
        now = strftime("%H:%M", gmtime())
        print(f"Finished Grid Search for {model_name}: {now}")

        full_results[model_name] = pd.DataFrame(clf.cv_results_)
        
    important_result = pd.DataFrame(scores,columns=['model','best_score','best_params'])
    t1 = perf_counter() - t0
    print(f'Tempo de execução: {(t1/60): 0.1f} minutos')
    return important_result, full_results


def calculate_confusion_matrix(y_test, y_pred):
    """
    This function recieves a y_test and y_predictes arrays
    it calculates and return the misclassification rate
    and the confusion matrix
    """
    cm = confusion_matrix(y_test, y_pred)
    mis_rate = (cm[[1],[0]].flat[0] + cm[[0],[1]].flat[0])/len(y_test)
    print(f"Misclassification rate: {mis_rate:.4f}")
    disp = ConfusionMatrixDisplay(confusion_matrix=cm)

    return mis_rate, disp


def plot_roc_auc_curve(y_test, y_pred, model_name, title, save_path):
    """
    This function recieves the y_test and y_pred (in probability)
    and return the auc, and plots the AUC ROC curve and saves it
    """
    # calculating auc
    auc = roc_auc_score(y_test, y_pred)

    # creating diagonal line
    r_proba = [0 for _ in range(len(y_test))]
    r_fpr, r_tpr, _ = roc_curve(y_test, r_proba)
    r_auc = roc_auc_score(y_test, r_proba)

    # creating ROC
    fpr, tpr, _ = roc_curve(y_test, y_pred)

    # Creating Graph

    #FigSize
    fig = plt.figure(figsize=(10,8))

    #Plot
    plt.plot(fpr, tpr, '.-', label = f'{model_name} (AUC = %0.3f)' % auc)
    plt.plot(r_fpr, r_tpr, '--', label = '')

    #Title
    plt.title(title)
    #Axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    #Show legend
    plt.legend(loc=4)
    #Show plot
    plt.show()

    fig.savefig(save_path, dpi=300)

    return auc

def make_label(x):
    if x >= 0.5:
        return 1.0
    return 0.0

def bart_auc_scorer(y_true, y_pred):
    # create y_pred
    #y_pred = model.predict(x_train)
    
    # creating y_pred with labels
    make_label_v = np.vectorize(make_label)
    y_pred = make_label_v(y_pred)

    return roc_auc_score(y_true, y_pred)

def get_error_and_auc(model, x, y_true, transform_prob_into_label=False):

    # predict y_pred
    y_pred = model.predict(x)

    # transform prov into label, if necessary
    if transform_prob_into_label:
        make_label_v = np.vectorize(make_label)
        y_pred = make_label_v(y_pred)
    
    # calculate error
    cm = confusion_matrix(y_true, y_pred)
    mis_rate = (cm[[1],[0]].flat[0] + cm[[0],[1]].flat[0])/len(y_true)

    # calculate auc
    auc = roc_auc_score(y_true, y_pred)

    # printing results
    print(f"Training Misclassification Rate: {mis_rate:.4f}")
    print(f"Training AUC: {auc:.4f}")

    return mis_rate, auc