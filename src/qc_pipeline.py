import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from sklearn.ensemble import IsolationForest
def load_data(data):
    raw = pd.read_csv(data, sep="\t", index_col=0)  # genes as rows
    df_sample_counts = raw.T  # flip: cells become rows
    df_sample_counts = df_sample_counts.reset_index().rename(columns={"index": "cell_id"})
    print(df_sample_counts.head())
    print(df_sample_counts.shape)
    return df_sample_counts

def calculate_total_counts(df_sample_counts):
#get total counts per cell 
    gene_counts_df=df_sample_counts.drop(columns={"cell_id"})
    total_counts=gene_counts_df.sum(axis=1) 
    df_sample_counts["total_counts"]=total_counts
    print(df_sample_counts[["cell_id","total_counts"]])
    return df_sample_counts, gene_counts_df

def calculate_genes_detected(df_sample_counts,gene_counts_df):
#QC metric which is the number of genes detected per cell
    QC_metric=gene_counts_df>0 
    genes_detected=QC_metric.sum(axis=1) 
    df_sample_counts["genes_detected"]=genes_detected
    print(df_sample_counts[["cell_id","total_counts","genes_detected"]])
    return df_sample_counts 

def find_mitochondrial_genes(df_sample_counts):
# tells us what percent of the cell's total count count comes from mito genes
    mito_detect=[]
    for column in df_sample_counts.columns: 
        if column.startswith("MT-"):
            mito_detect.append(column)
    print(mito_detect)
    return mito_detect

def calculate_mito_percentage(df_sample_counts,mito_detect):
    mito_counts=df_sample_counts[mito_detect].sum(axis=1)
    df_sample_counts["mito_counts"]=mito_counts
    print(df_sample_counts[["cell_id","total_counts","mito_counts"]])
    percentage_mito=(df_sample_counts["mito_counts"]/df_sample_counts["total_counts"])*100
    df_sample_counts["percentage_mito"]=percentage_mito
    print(df_sample_counts[["cell_id","total_counts","percentage_mito"]])
    return df_sample_counts

def assign_qc_status(df_sample_counts, mito_thresh=20, gene_thresh=200, count_thresh=500):
    df_sample_counts["qc_status"] = "pass QC"
    df_sample_counts.loc[
        (df_sample_counts["percentage_mito"] > mito_thresh) |
        (df_sample_counts["genes_detected"] < gene_thresh) |
        (df_sample_counts["total_counts"] < count_thresh),
        "qc_status"
    ] = "fail QC"
    print(df_sample_counts[["cell_id", "total_counts", "qc_status"]])
    return df_sample_counts  
def assign_qc_status_isolation_forest(df_sample_counts, contamination=0.1, random_state=42):
    features = df_sample_counts[["total_counts", "genes_detected", "percentage_mito"]]
    iso_forest = IsolationForest(contamination=contamination, random_state=random_state)
    predictions = iso_forest.fit_predict(features)
    df_sample_counts["anomaly_score"] = iso_forest.decision_function(features)  # lower = more anomalous
    df_sample_counts["qc_status"] = ["fail QC" if p == -1 else "pass QC" for p in predictions]
    print(df_sample_counts[["cell_id", "total_counts", "genes_detected", "percentage_mito", "anomaly_score", "qc_status"]])
    return df_sample_counts

def save_results(df_sample_counts):
    qc_summary=df_sample_counts[["cell_id","total_counts","genes_detected","mito_counts","percentage_mito","qc_status"]]
    qc_summary.to_csv("results/QC_summary.csv", index=False)

    passed_cells=qc_summary.drop(qc_summary[qc_summary["qc_status"]=="fail QC"].index)
    print(passed_cells[["cell_id","qc_status"]])
    passed_cells.to_csv("results/passed_cells.csv", index=False)
    return passed_cells

def creating_plot(data,x_axis,y_axis,y_label,title,output_path):
    plt.figure(figsize=(5.5,5.5))
    sns.barplot(x=x_axis,y=y_axis,data=data)
    plt.xlabel('cell ID', fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.xticks(rotation=45)
    plt.title(title,fontsize=15,fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path,dpi=300)
    plt.show()

def plot_distribution(data, column, x_label, title, output_path, threshold=None):
    plt.figure(figsize=(6, 5))
    sns.histplot(data[column], bins=30)
    if threshold is not None:
        plt.axvline(threshold, color='red', linestyle='--', label=f'cutoff = {threshold}')
        plt.legend()
    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.ylabel("Number of Cells", fontsize=14, fontweight='bold')
    plt.title(title, fontsize=15, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

def create_qc_status_plot(df_sample_counts):
    qc_counts = df_sample_counts["qc_status"].value_counts()

    plt.figure(figsize=(5.5, 5.5))

    sns.barplot(
        x=qc_counts.index,
        y=qc_counts.values
    )

    plt.xlabel("QC Status", fontsize=14, fontweight="bold")
    plt.ylabel("Number of Cells", fontsize=14, fontweight="bold")
    plt.title("QC Status Counts", fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig("results/qc_status_counts.png", dpi=300)
    plt.show()

def create_all_plots(df_sample_counts):
    #so that all plotting happens and can be called in one spot
    plot_distribution(
        df_sample_counts,
        column="percentage_mito",
        x_label="Percentage mitochondrial counts",
        title="Mitochondrial Percentage Distribution",
        output_path="results/mitochondrial_percentage_plot.png",
        threshold=20
    )
    plot_distribution(
        df_sample_counts,
        column="total_counts",
        x_label="Total counts per cell",
        title="Total Counts Distribution",
        output_path="results/total_counts_plot.png",
        threshold=500
    )
    plot_distribution(
        df_sample_counts,
        column="genes_detected",
        x_label="Genes detected per cell",
        title="Genes Detected Distribution",
        output_path="results/genes_detected_plot.png",
        threshold=200
    )
    create_qc_status_plot(df_sample_counts)

def write_report(df_sample_counts,passed_cells):
#printing QC summary counts
    Total_cells=len(df_sample_counts) #this gets the number of columns
    passed_count=len(passed_cells)
    failed_count=Total_cells-passed_count
    avg_mito=df_sample_counts["percentage_mito"].mean()
    avg_total=df_sample_counts["total_counts"].mean()
    print(f"Total number of cells is,{Total_cells}")
    print(f"Total number of passed count is,{passed_count}")
    print(f"Total number of failed count is,{failed_count}")
    print(f"Average mitochondrial percentage is,{avg_mito}")
    print(f"Average total count is,{avg_total}")

    with open("results/qc_report.txt", "w") as file:
        file.write("Single-Cell RNA-seq QC Report\n\n")
        file.write(f"Total cells: {Total_cells}\n")
        file.write(f"Passed QC: {passed_count}\n")
        file.write(f"Failed QC: {failed_count}\n")
        file.write(f"Average mitochondrial percentage: {round(avg_mito, 2)}\n")
        file.write(f"Average total counts: {round(avg_total, 2)}\n")

    qc_passed_summary=df_sample_counts[["cell_id","total_counts","genes_detected","mito_counts","percentage_mito","qc_status"]]
    qc_passed_summary.to_csv("results/QC_passed_summary.csv", index=False)



def main():
    os.makedirs("results", exist_ok=True)
    df_sample_counts = load_data("data/CRC_subset.txt")
    df_sample_counts, gene_counts_df = calculate_total_counts(df_sample_counts)
    df_sample_counts = calculate_genes_detected(df_sample_counts, gene_counts_df)
    mito_detect = find_mitochondrial_genes(df_sample_counts)
    df_sample_counts = calculate_mito_percentage(df_sample_counts, mito_detect)
    df_sample_counts = assign_qc_status_isolation_forest(df_sample_counts, contamination=0.1)
    passed_cells = save_results(df_sample_counts)
    create_all_plots(df_sample_counts)
    write_report(df_sample_counts, passed_cells)

if __name__ == "__main__":
    main()





