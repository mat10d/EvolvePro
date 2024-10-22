import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from evolvepro.src.data import load_experimental_data

pd.options.mode.chained_assignment = None  # default='warn'

def read_dms_data(directory, datasets, model, experiment, group_columns, aggregate_columns, 
              file_pattern="{dataset}_{model}_{experiment}.csv"):
    """
    Read and process data from multiple CSV files.
    
    Args:
    directory (str): Directory containing the CSV files
    datasets (list): List of dataset names
    model (str): Model name
    experiment (str): Experiment name
    group_columns (list): Columns to group by
    aggregate_columns (list): Columns to aggregate
    file_pattern (str): File name pattern for CSV files
    
    Returns:
    pd.DataFrame: Processed and concatenated data from all datasets
    """
    all_dfs = []
    
    for dataset in datasets:
        file_name = file_pattern.format(dataset=dataset, model=model, experiment=experiment)
        file_path = os.path.join(directory, file_name)
        
        try:
            df = pd.read_csv(file_path)
            df = process_dataframe(df, group_columns, aggregate_columns)
            df['dataset'] = dataset  
            df['model'] = model  
            df['experiment'] = experiment
            all_dfs.append(df)
        except FileNotFoundError:
            print(f"File {file_name} not found. Skipping...")
        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}")
    
    return pd.concat(all_dfs, ignore_index=True) if all_dfs else pd.DataFrame()

def read_exp_data(round_base_path, round_file_names_single, wt_fasta_path, round_file_names_multi=None):
    """
    Read and process experimental data from multiple files.

    Args:
    round_base_path (str): Base path to the data directory containing the excel files.
    round_file_names_single (list): List of single mutant round file names.
    wt_fasta_path (str): Path to the wild-type FASTA file.
    round_file_names_multi (list): List of multi mutant round file names.

    Returns:
    pd.DataFrame: Processed and concatenated experimental data
    """

    # Load experimental data
    all_experimental_data = []
    for round_file_name in round_file_names_single:
        experimental_data = load_experimental_data(round_base_path, round_file_name, wt_fasta_path, single_mutant=True)
        all_experimental_data.append(experimental_data)

    if round_file_names_multi is not None:
        for round_file_name in round_file_names_multi:
            experimental_data = load_experimental_data(round_base_path, round_file_name, wt_fasta_path, single_mutant=False)
            all_experimental_data.append(experimental_data)
    
    processed_dfs = []
    # Process each round's data
    for round_num, df in enumerate(all_experimental_data, start=1):
        df_copy = df.copy()
        
        # Set iteration for WT in first round, exclude WT from subsequent rounds
        if round_num == 1:
            df_copy.loc[df_copy['updated_variant'] == 'WT', 'iteration'] = 0
        else:
            df_copy = df_copy[df_copy['updated_variant'] != 'WT']
        
        df_copy.loc[df_copy['updated_variant'] != 'WT', 'iteration'] = round_num
        df_copy['iteration'] = df_copy['iteration'].astype(float)
        df_copy.rename(columns={'updated_variant': 'variant'}, inplace=True)
        
        processed_dfs.append(df_copy)

    # Combine all processed dataframes
    combined_df = pd.concat(processed_dfs, ignore_index=True)

    return combined_df

def process_dataframe(df, group_columns, aggregate_columns):
    """
    Process a dataframe by grouping and aggregating columns.

    Args:
    df (pd.DataFrame): Input dataframe
    group_columns (list): Columns to group by
    aggregate_columns (list): Columns to aggregate

    Returns:
    pd.DataFrame: Processed dataframe
    """
    
    df.replace("None", np.nan, inplace=True)
    df = df.dropna(subset=aggregate_columns, how='all')
    df[aggregate_columns] = df[aggregate_columns].apply(pd.to_numeric, errors='coerce')
    
    grouped = df.groupby(group_columns)
    stats = grouped[aggregate_columns].agg(['mean', 'std'])
    stats.columns = [f'{col}_{stat}' for col, stat in stats.columns]
    return stats.reset_index()

def filter_dataframe(df, conditions, output_dir=None, output_file=None):
    """
    Filter a dataframe based on conditions.

    Args:
    df (pd.DataFrame): Input dataframe
    conditions (dict): Dictionary of column-value pairs to filter on
    output_dir (str, optional): Output directory for the CSV file
    output_file (str, optional): Output file name    

    Returns:
    pd.DataFrame: Filtered dataframe
    """
    for column, value in conditions.items():
        if isinstance(value, list):
            filtered_df = df[df[column].isin(value)]
        else:
            filtered_df = df[df[column] == value]

    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        filtered_df.to_csv(os.path.join(output_dir, output_file), index=False)
        
    return filtered_df

def save_dataframe(df, output_dir=None, output_file=None):
    """
    Save a dataframe.

    Args:
    df (pd.DataFrame): Input dataframe
    output_dir (str, optional): Output directory for the CSV file
    output_file (str, optional): Output file name

    """
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(os.path.join(output_dir, output_file), index=False)
        
def apply_labels(df, column, prefix='', suffix='', value_column=None, format_string='{}'):
    """
    Apply labels to a DataFrame column.
    
    Args:
    df (pd.DataFrame): Input dataframe
    column (str): Column name to be created
    prefix (str, optional): Prefix for the label
    suffix (str, optional): Suffix for the label
    value_column (str, optional): Column to use for the label value
    format_string (str): Format string for the label value

    Returns:
    pd.DataFrame: DataFrame with the new column
    """   
    if value_column is None:
        df[column] = prefix + df.index.map(lambda x: format_string.format(x)) + suffix
    else:
        df[column] = prefix + df[value_column].map(lambda x: format_string.format(x)) + suffix
    return df

def load_external_data(file_path, label=None, rename_columns=None):
    """
    Load external data from a CSV file, optionally add a label column, and rename columns if specified.
    
    Args:
    file_path (str): Path to the CSV file
    label (str, optional): Label to be added as a new column
    rename_columns (dict, optional): Dictionary of column names to rename, e.g., {'old_name': 'new_name'}
    
    Returns:
    pd.DataFrame: Loaded and processed dataframe
    """
    df = pd.read_csv(file_path)
    
    if label:
        df['label'] = label

    if rename_columns:
        df = df.rename(columns=rename_columns)
    
    return df

def concatenate_dataframes(dataframes, output_dir = None, output_file = None):
    """
    Concatenate a list of dataframes and optionally save the result to a CSV file.
    
    Args:
    dataframes (list): List of dataframes to concatenate
    output_dir (str, optional): Output directory for the CSV file
    output_file (str, optional): Output file name
    
    Returns:
    pd.DataFrame: Concatenated dataframe
    """
    concatenated_df = pd.concat(dataframes, ignore_index=True)
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        concatenated_df.to_csv(os.path.join(output_dir, output_file), index=False)
    
    return concatenated_df

def plot_comparison(concatenated_df, palette=None, variable='fitness_binary_percentage_mean', title=None, output_dir=None, output_file=None):
    """
    Generate plots from a concatenated dataframe.
    
    Args:
    concatenated_df (pd.DataFrame): Dataframe containing the data to plot
    palette (dict, optional): Custom color palette for the plots
    variable (str): Name of the variable to plot on y-axis
    output_dir (str, optional): Directory to save the plots
    output_file (str, optional): Base name for the output files
    
    Returns:
    None
    """
    if palette is None:
        palette_colors = sns.color_palette("tab10")
    else:
        palette_colors = palette
    
    # Plot 1: Bar plot by dataset
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=concatenated_df, x='dataset', y=variable, 
                     hue='label', palette=palette_colors, errorbar=None, alpha=0.75)
    plt.xlabel('Dataset')
    plt.ylabel(f'{variable.replace("_", " ").title()}')
    plt.title(title)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_by_dataset.png"), dpi=300, bbox_inches='tight')
    
    plt.show()
    
    # Plot 2: Bar plot by label with swarm plot
    plt.figure(figsize=(7, 6))
    ax = sns.barplot(data=concatenated_df, x='label', y=variable,
                     palette=palette_colors, errorbar=None, alpha=0.75)
    sns.swarmplot(data=concatenated_df, x='label', y=variable, 
                  size=4, color="black")
    plt.xlabel('Model')
    plt.ylabel(f'{variable.replace("_", " ").title()}')
    plt.title(title)
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    if output_dir and output_file:
        plt.savefig(os.path.join(output_dir, f"{output_file}_by_model.png"), dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_grid_search_bar(df, variable='fitness_binary_percentage_mean', strategy_column=None, title=None, output_dir=None, output_file=None):
    """
    Generate plots from a dataframe to compare different strategies.
    
    Args:
    df (pd.DataFrame): Dataframe containing the data to plot
    variable (str): Name of the variable to plot on y-axis
    strategy_column (str): Name of the column containing different strategies
    title (str, optional): Title for the plot
    output_dir (str, optional): Directory to save the plots
    output_file (str, optional): Base name for the output files
    
    Returns:
    None
    """
    if strategy_column is None:
        raise ValueError("strategy_columns must be a list of exactly two column names")

    round_num = df['round_num'].iloc[0]
    grouped = df.groupby(['dataset', strategy_column])[variable].mean().unstack()

    # Plot
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(data=df, x='dataset', y=variable, hue=strategy_column, errorbar=None, alpha=0.75)
    
    if title is None:
        title = f'{variable.replace("_", " ").title()} by {strategy_column.replace("_", " ").title()} ({round_num} rounds)'
    
    ax.set_title(title)
    ax.set_xlabel('Dataset')
    ax.set_ylabel(variable.replace("_", " ").title())
    ax.legend(title=strategy_column.replace("_", " ").title(), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_grid_bar.png"), dpi=300, bbox_inches='tight')
    
    plt.show()

    # Count the occurrences of each strategy being the best for each dataset
    winning_counts = grouped.apply(lambda x: x.idxmax(), axis=1).value_counts()

    # Print counts of winning strategies
    print("\nCounts of Winning Strategies:")
    print(winning_counts)

def plot_grid_search_heatmap(df, variable='fitness_binary_percentage_mean', strategy_columns=None, title=None, output_dir=None, output_file=None):
    """
    Generate a heatmap from a dataframe to compare the intersection of two strategies.
    
    Args:
    df (pd.DataFrame): Dataframe containing the data to plot
    variable (str): Name of the variable to average and plot
    strategy_columns (list): List of two columns containing different strategies
    title (str, optional): Title for the plot
    output_dir (str, optional): Directory to save the plot
    output_file (str, optional): Base name for the output file
    
    Returns:
    None
    """
    if strategy_columns is None or len(strategy_columns) != 2:
        raise ValueError("strategy_columns must be a list of exactly two column names")

    # Extract round number
    round_num = df['round_num'].iloc[0]

    # Group by the two strategy columns and calculate the mean of the variable
    grouped = df.groupby(strategy_columns)[variable].mean().unstack()

    # Plot
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(grouped, cmap='viridis', annot=True, fmt=".2f", linewidths=0.5)
    
    if title is None:
        title = f'Average {variable.replace("_", " ").title()} by {strategy_columns[0].replace("_", " ").title()} and {strategy_columns[1].replace("_", " ").title()} ({round_num} rounds)'
    
    ax.set_title(title)
    ax.set_xlabel(strategy_columns[1].replace("_", " ").title())
    ax.set_ylabel(strategy_columns[0].replace("_", " ").title())
    plt.tight_layout()
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_grid_heatmap.png"), dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_by_round(df, variable='fitness_binary_percentage_mean', output_dir=None, output_file=None):
    """
    Plot round-by-round comparison for each dataset.
    
    Args:
    df (pd.DataFrame): Dataframe containing the data to plot
    variable (str): Name of the variable to plot on y-axis
    output_dir (str, optional): Directory to save the plot
    output_file (str, optional): Base name for the output file
    
    Returns:
    None
    """
    plt.figure(figsize=(10, 6))
    ax = plt.gca()

    # Generate a color palette for the unique datasets
    datasets = df['dataset'].unique()
    color_palette = sns.color_palette("tab10", n_colors=len(datasets))
    color_map = dict(zip(datasets, color_palette))

    for dataset in datasets:
        dataset_df = df[df['dataset'] == dataset]
        x_values = dataset_df['round_num']
        y_values = dataset_df[variable]
        color = color_map[dataset]
        ax = sns.lineplot(x=x_values, y=y_values, ax=ax, marker='o', color=color, label=dataset)
 
    ax.set_xlabel('Number of Iterations')
    ax.set_ylabel(variable.replace('_', ' ').title())
    ax.set_title(f'{variable.replace("_", " ").title()} by Iterations')
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_by_round.png"), dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_by_round_split(df, variable='fitness_binary_percentage_mean', split_variable='num_mutants_per_round', output_dir=None, output_file=None):
    """
    Plot round-by-round comparison with separate subplots for each dataset, with lines split by a specified variable.
    Includes a shared legend for all subplots.
    
    Args:
    df (pd.DataFrame): Dataframe containing the data to plot
    variable (str): Name of the variable to plot on y-axis
    split_variable (str): Name of the variable to split lines by within each subplot
    output_dir (str, optional): Directory to save the plot
    output_file (str, optional): Base name for the output file
    
    Returns:
    None
    """
    datasets = df['dataset'].unique()
    n_datasets = len(datasets)
    
    # Create color palette for split variable
    split_values = sorted(df[split_variable].unique())
    color_palette = sns.color_palette("tab10", n_colors=len(split_values))
    color_map = dict(zip(split_values, color_palette))
    
    # Calculate number of rows and columns for subplots
    n_cols = 3  # You can adjust this
    n_rows = (n_datasets + n_cols - 1) // n_cols
    
    # Create figure and subplots
    fig = plt.figure(figsize=(20, 4*n_rows))
    gs = fig.add_gridspec(n_rows, n_cols)
    
    # Create plots
    for idx, dataset in enumerate(datasets):
        row = idx // n_cols
        col = idx % n_cols
        ax = fig.add_subplot(gs[row, col])
        
        dataset_df = df[df['dataset'] == dataset]
        
        for split_value in split_values:
            subset_df = dataset_df[dataset_df[split_variable] == split_value]
            
            sns.lineplot(data=subset_df, x='round_num', y=variable, 
                        marker='o', ax=ax, color=color_map[split_value],
                        label=f'{split_variable}: {split_value}')
        
        ax.set_xlabel('Number of Iterations')
        ax.set_ylabel(variable.replace('_', ' ').title())
        ax.set_title(dataset)
        ax.legend().remove()  # Remove individual legends
        
    # Remove empty subplots if any
    for idx in range(len(datasets), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        fig.delaxes(plt.subplot(gs[row, col]))
    
    # Create shared legend
    lines = []
    labels = []
    for split_value in split_values:
        lines.append(plt.Line2D([0], [0], color=color_map[split_value], marker='o', linestyle='-'))
        labels.append(f'{split_variable}: {split_value}')
    
    fig.legend(lines, labels, loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    # Add overall title
    fig.suptitle(f'{variable.replace("_", " ").title()} by Iterations', y=1.02)
    
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.98, 1]) 
    
    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_by_round_split_{split_variable}.png"), 
                    dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_variants_by_iteration(df, fitness_column='fitness', output_dir=None, output_file=None):
    """
    Simple bar plot of variants grouped by iteration.
    
    Args:
    df: DataFrame with 'variant', 'iteration', and fitness column
    fitness_column: Column name containing fitness values
    output_dir: Directory to save the plot
    output_file: Filename for the saved plot
    """
    # sort the dataframe by iteration and the fitness column within each iteration
    df['iteration'] = df['iteration'].astype(int)
    df = df.sort_values(['iteration', fitness_column], ascending=[True, True])
    df = df.reset_index(drop=True)

    plt.figure(figsize=(12, 6))
    
    # Plot each variant in the order of the dataframe, colored by iteration
    for iteration, group in df.groupby('iteration'):
        plt.bar(group.index, group[fitness_column], label=f"Round {iteration}")

    # Customize
    plt.xticks(df.index, df['variant'], rotation=90)
    plt.ylabel(fitness_column.capitalize())
    plt.legend()
    
    plt.tight_layout()

    if output_dir and output_file:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f"{output_file}_by_iteration.png"), dpi=300, bbox_inches='tight')

    plt.show()
