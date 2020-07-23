"""Path definitions"""
from pathlib import Path

# Path to this repository
PROJECT_DIR = Path(__file__).parents[1]
ANALYSIS_DIR = PROJECT_DIR / "sra_experiment"
METADATA_DIR = ANALYSIS_DIR / "data" / "metadata"
SAMPLE_ANNOT = METADATA_DIR / "sample_groups.txt"
GENE_PAO1_ANNOT = METADATA_DIR / "PAO1_ID_2_PA14_ID_PAO1ref.csv"
GENE_PA14_ANNOT = METADATA_DIR / "PA14_ID_2_PAO1_ID_PA14ref.csv"

# Path to local directory where data files will be stored
LOCAL_DIR = Path.home() / "Documents" / "Data" / "Generic_expression_patterns_test"

# Location where data is stored
TEMPLATE_FILE = LOCAL_DIR / "recount2_template_data.tsv"
COMPENDIUM_FILE = LOCAL_DIR / "recount2_compendium_data.tsv"
NORMAL_COMPENDIUM_FILE = LOCAL_DIR / "normalized_recount2_compendium_data.tsv"
SHARED_GENES_FILE = LOCAL_DIR / "shared_genes_human_test.pickle"
SCALER_TRANSFORM_FILE = LOCAL_DIR / "scaler_transform_human_test.pickle"

# Location of reference gene rank file
REF_GENE_FILE = (
    "https://raw.githubusercontent.com/maggiecrow/DEprior/master/DE_Prior.txt"
)

