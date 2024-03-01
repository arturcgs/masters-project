import pandas as pd
from sklearn.preprocessing import StandardScaler

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