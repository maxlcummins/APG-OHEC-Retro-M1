import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.utils import resample
from joblib import Parallel, delayed


data = pd.read_table("Pairwise_SNP_and_allelic_distances_with_source.txt.zip")


def create_network(data, threshold):
    # Initialize a graph
    G = nx.Graph()

    # Set to track added edges
    added_edges = set()

    # Add edges with weights
    for _, row in data.iterrows():
        if row['SNPs'] <= threshold:
            sample_pair = tuple(sorted([row['Sample 1'], row['Sample 2']]))
            if sample_pair not in added_edges:
                G.add_edge(row['Sample 1'], row['Sample 2'], weight=row['SNPs'])
                added_edges.add(sample_pair)

    # Remove isolated nodes (nodes with degree 0)
    G.remove_nodes_from(list(nx.isolates(G)))

    # Remove clusters with less than 2 members
    clusters = list(nx.connected_components(G))
    for cluster in clusters:
        if len(cluster) < 2:
            G.remove_nodes_from(cluster)

    return G

def plot_network(G, threshold, ax, cluster_labels=None, silhouette_avg=None):
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato")

    # Get edge weights and normalize them for width scaling
    weights = nx.get_edge_attributes(G, 'weight')
    max_weight = max(weights.values())
    min_weight = min(weights.values())
    widths = [5 * (weight - min_weight) / (max_weight - min_weight) + 0.1 for weight in weights.values()]

    # Assign colors to nodes based on cluster labels and silhouette scores
    if cluster_labels is not None and silhouette_avg is not None:
        unique_labels = set(cluster_labels)
        colors = [plt.cm.jet(float(i) / len(unique_labels)) for i in range(len(unique_labels))]
        color_map = dict(zip(unique_labels, colors))
        node_colors = [color_map[cluster_labels[i]] for i in range(len(G.nodes))]
        nx.draw_networkx_nodes(G, pos, node_size=50, node_color=node_colors, ax=ax)
    else:
        nx.draw_networkx_nodes(G, pos, node_size=50, node_color='skyblue', ax=ax)

    nx.draw_networkx_edges(G, pos, edge_color='gray', width=widths, ax=ax)
    ax.set_title(f'Network with SNP Distance < {threshold}')
    ax.axis('off')

# Vectorized function to calculate WCSS
def calculate_wcss(square_distance_matrix, cluster_labels):
    unique_clusters = np.unique(cluster_labels)
    wcss = 0
    for cluster in unique_clusters:
        cluster_points = square_distance_matrix[cluster_labels == cluster]
        centroid = np.mean(cluster_points, axis=0)
        wcss += np.sum((cluster_points - centroid) ** 2)
    return wcss

# Vectorized cohesion and separation calculation
def calculate_cohesion_separation(square_distance_matrix, cluster_labels):
    unique_clusters = np.unique(cluster_labels)
    centroids = np.array([np.mean(square_distance_matrix[cluster_labels == cluster], axis=0) for cluster in unique_clusters])

    cohesion = np.mean([np.sum(square_distance_matrix[cluster_labels == cluster]) for cluster in unique_clusters])
    separation = np.mean([np.linalg.norm(centroids[i] - centroids[j]) for i in range(len(centroids)) for j in range(i + 1, len(centroids))])

    return cohesion, separation

# Function to evaluate clustering quality
def evaluate_clustering(threshold, linkage_matrix, grp):
    cluster_labels = fcluster(linkage_matrix, threshold, criterion='distance')
    n_clusters = len(np.unique(cluster_labels))

    if n_clusters > 1:
        square_distance_matrix = squareform(grp)
        silhouette = silhouette_score(square_distance_matrix, cluster_labels, metric='precomputed')
        ch_score = calinski_harabasz_score(square_distance_matrix, cluster_labels)
        db_score = davies_bouldin_score(square_distance_matrix, cluster_labels)

        wcss = calculate_wcss(square_distance_matrix, cluster_labels)
        cohesion, separation = calculate_cohesion_separation(square_distance_matrix, cluster_labels)
    else:
        silhouette, ch_score, db_score, wcss, cohesion, separation = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    return silhouette, ch_score, db_score, wcss, cohesion, separation

# Parallelized permutation test
def parallel_compare_thresholds(condensed_distance_matrix, linkage_matrix, threshold_1, threshold_2, n_permutations=1000):
    # Helper function to evaluate clustering in parallel
    def evaluate_permutation(_):
        permuted_data = resample(condensed_distance_matrix)
        perm_linkage_matrix = linkage(permuted_data, method='ward')
        return evaluate_clustering(threshold_1, perm_linkage_matrix, permuted_data), evaluate_clustering(threshold_2, perm_linkage_matrix, permuted_data)

    score_1_sil, score_1_ch, score_1_db, score_1_wcss, score_1_cohesion, score_1_separation = evaluate_clustering(threshold_1, linkage_matrix, condensed_distance_matrix)
    score_2_sil, score_2_ch, score_2_db, score_2_wcss, score_2_cohesion, score_2_separation = evaluate_clustering(threshold_2, linkage_matrix, condensed_distance_matrix)

    observed_diff_sil = score_2_sil - score_1_sil
    observed_diff_ch = score_2_ch - score_1_ch
    observed_diff_db = score_1_db - score_2_db
    observed_diff_wcss = score_1_wcss - score_2_wcss
    observed_diff_cohesion = score_1_cohesion - score_2_cohesion
    observed_diff_separation = score_2_separation - score_1_separation

    # Run permutations in parallel
    results = Parallel(n_jobs=-1)(delayed(evaluate_permutation)(i) for i in range(n_permutations))

    # Accumulate permutation differences
    perm_diffs_sil, perm_diffs_ch, perm_diffs_db, perm_diffs_wcss, perm_diffs_cohesion, perm_diffs_separation = [], [], [], [], [], []
    for res_1, res_2 in results:
        score_1_perm_sil, score_1_perm_ch, score_1_perm_db, score_1_perm_wcss, score_1_perm_cohesion, score_1_perm_separation = res_1
        score_2_perm_sil, score_2_perm_ch, score_2_perm_db, score_2_perm_wcss, score_2_perm_cohesion, score_2_perm_separation = res_2

        if not np.isnan(score_1_perm_sil) and not np.isnan(score_2_perm_sil):
            perm_diffs_sil.append(score_2_perm_sil - score_1_perm_sil)
        if not np.isnan(score_1_perm_ch) and not np.isnan(score_2_perm_ch):
            perm_diffs_ch.append(score_2_perm_ch - score_1_perm_ch)
        if not np.isnan(score_1_perm_db) and not np.isnan(score_2_perm_db):
            perm_diffs_db.append(score_1_perm_db - score_2_perm_db)
        if not np.isnan(score_1_perm_wcss) and not np.isnan(score_2_perm_wcss):
            perm_diffs_wcss.append(score_1_perm_wcss - score_2_perm_wcss)
        if not np.isnan(score_1_perm_cohesion) and not np.isnan(score_2_perm_cohesion):
            perm_diffs_cohesion.append(score_1_perm_cohesion - score_2_perm_cohesion)
        if not np.isnan(score_1_perm_separation) and not np.isnan(score_2_perm_separation):
            perm_diffs_separation.append(score_2_perm_separation - score_1_perm_separation)

    def calculate_p_value(observed_diff, perm_diffs, reverse=False):
      if len(perm_diffs) > 0:
          if reverse:
              return np.sum(observed_diff <= perm_diffs) / len(perm_diffs)  # For lower-is-better metrics
          else:
              return np.sum(observed_diff >= perm_diffs) / len(perm_diffs)  # For higher-is-better metrics
      else:
          return np.nan


    # p-value calculation
    p_value_sil = calculate_p_value(observed_diff_sil, perm_diffs_sil)  # Higher is better
    p_value_ch = calculate_p_value(observed_diff_ch, perm_diffs_ch)     # Higher is better
    p_value_db = calculate_p_value(observed_diff_db, perm_diffs_db, reverse=True)  # Lower is better
    p_value_wcss = calculate_p_value(observed_diff_wcss, perm_diffs_wcss, reverse=True)  # Lower is better
    p_value_cohesion = calculate_p_value(observed_diff_cohesion, perm_diffs_cohesion, reverse=True)  # Lower is better
    p_value_separation = calculate_p_value(observed_diff_separation, perm_diffs_separation)  # Higher is better


    return {
        'observed_diff_sil': observed_diff_sil, 'p_value_sil': p_value_sil,
        'observed_diff_ch': observed_diff_ch, 'p_value_ch': p_value_ch,
        'observed_diff_db': observed_diff_db, 'p_value_db': p_value_db,
        'observed_diff_wcss': observed_diff_wcss, 'p_value_wcss': p_value_wcss,
        'observed_diff_cohesion': observed_diff_cohesion, 'p_value_cohesion': p_value_cohesion,
        'observed_diff_separation': observed_diff_separation, 'p_value_separation': p_value_separation
    }

# Example: Running the analysis
threshold_1 = 20
threshold_2 = 100

st_groups = data["Sequence Type - Sample 1"].unique()

for ST in st_groups:
    grp = data[(data['Sequence Type - Sample 1'] == ST)]
    samples = pd.concat([grp['Sample 1'], grp['Sample 2']]).unique()
    n_samples = len(samples)
    sample_index = {sample: idx for idx, sample in enumerate(samples)}

    # Initialize a distance matrix
    distance_matrix = np.zeros((n_samples, n_samples))

    for _, row in grp.iterrows():
        i = sample_index[row['Sample 1']]
        j = sample_index[row['Sample 2']]
        distance_matrix[i, j] = row['SNPs']
        distance_matrix[j, i] = row['SNPs']

    # Convert the distance matrix to condensed form
    condensed_distance_matrix = squareform(distance_matrix, checks=False)

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distance_matrix, method='ward')

    # Perform parallel permutation test with multiple metrics
    results = parallel_compare_thresholds(condensed_distance_matrix, linkage_matrix, threshold_1, threshold_2)

    print(ST, results)
      # Create the networks
    G1 = create_network(grp, threshold_1)
    G2 = create_network(grp, threshold_2)

    # Plot the networks
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [6, 1], 'width_ratios': [1, 1]})

    # Plot network 1
    plot_network(G1, threshold_1, axes[0, 0])

    # Plot network 2
    plot_network(G2, threshold_2, axes[0, 1])

    # Turn off the axis where the table will go
    axes[1, 0].axis('off')
    axes[1, 1].axis('off')

    # Create a table with the metrics
    metrics_data = [
        ["Metric", "Observed Difference", "p-value"],
        ["Silhouette", f"{results['observed_diff_sil']:.4f}", f"{results['p_value_sil']:.4f}"],
        ["Calinski-Harabasz", f"{results['observed_diff_ch']:.4f}", f"{results['p_value_ch']:.4f}"],
        ["Davies-Bouldin", f"{results['observed_diff_db']:.4f}", f"{results['p_value_db']:.4f}"],
        ["WCSS", f"{results['observed_diff_wcss']:.4f}", f"{results['p_value_wcss']:.4f}"],
        ["Cohesion", f"{results['observed_diff_cohesion']:.4f}", f"{results['p_value_cohesion']:.4f}"],
        ["Separation", f"{results['observed_diff_separation']:.4f}", f"{results['p_value_separation']:.4f}"]
    ]

    # Create the table and add it to the figure, spanning the last two axes
    table = axes[1, 1].table(cellText=metrics_data, colWidths=[0.3, 0.4, 0.4], loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.2)

    # Adjust the layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Set the figure title
    plt.suptitle(f"Cluster Quality Metrics for {ST}", fontsize=14, y=0.98)

    # Save and show the plot
    plt.savefig(f"{ST}_metrics_comparison.pdf", format="pdf", bbox_inches="tight")
    plt.show()
