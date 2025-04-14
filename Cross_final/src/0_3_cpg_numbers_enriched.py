# %%
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.expand_frame_repr', False)
# %%
def load_cpg_enrichment_data(file_path):
    """
    Loads the CpG enrichment data from the specified CSV file into a pandas DataFrame.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        pandas.DataFrame: A DataFrame containing the CpG enrichment data.
    """
    try:
        df = pd.read_csv(file_path)
        return df
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# %%
STUDY = "Neu"
base_dir = "D:/Github/SRF_MeCP2_cut_tag/iterative_alternative/results/no_dedup/cpg_enrichment"
file_path = f"{base_dir}/{STUDY}/broad/cpg_enrichment_1_rep_in_peaks/lists/up_enriched_signal_2_exo_over_20.csv"
cpg_data = load_cpg_enrichment_data(file_path)
print(f"Shape of cpg_data: {cpg_data.shape}")
print(f"Number of enriched CpGs: {len(cpg_data)}")
cpg_data.head()

# %%
STUDY = "NSC"
base_dir = "D:/Github/SRF_MeCP2_cut_tag/iterative_alternative/results/no_dedup/cpg_enrichment"
file_path = f"{base_dir}/{STUDY}/broad/cpg_enrichment_1_rep_in_peaks/lists/up_enriched_signal_2_exo_over_20.csv"
cpg_data = load_cpg_enrichment_data(file_path)
print(f"Shape of cpg_data: {cpg_data.shape}")
print(f"Number of enriched CpGs: {len(cpg_data)}")
cpg_data.head()
# %%
