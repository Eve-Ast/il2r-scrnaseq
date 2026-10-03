"""Shared helper functions for the IL-2 receptor scRNA-seq analysis
(loading, QC, normalization, dimensionality reduction, clustering,
cluster annotation, plotting and statistical tests)."""

# Imports directs
import gc
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.sparse import csc_matrix
from sklearn.decomposition import TruncatedSVD
from scipy.io import mmread
from scipy.stats import f_oneway, ttest_ind
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap

genes_of_interest = ["IL2RA", "IL2RB", "IL2RG"]  # Liste des gènes à tester

"""## Fonction pour l'import des données"""

def import_data(chemin_matrice, chemin_genes, chemin_barcodes):
    """
    Charge la matrice d'expression sparse et les métadonnées associées (noms des gènes et barcodes).

    Paramètres :
    - chemin_matrice : chemin vers le fichier .mtx.gz
    - chemin_genes : chemin vers le fichier des gènes (.txt.gz)
    - chemin_barcodes : chemin vers le fichier des barcodes (.txt.gz)

    Retour :
    - df : DataFrame sparse (gènes x cellules)
    - genes : tableau des noms de gènes
    - barcodes : tableau des noms de cellules
    - matrice : matrice CSC brute
    """

    # Charge la matrice sparse
    with gzip.open(chemin_matrice, 'rb') as f:
        matrice = mmread(f).tocsc()

    # Charge les gènes et les barcodes
    genes = pd.read_csv(chemin_genes, sep="\t", header=None)[1].values
    barcodes = pd.read_csv(chemin_barcodes, sep="\t", header=None)[0].values

    # Création d'un DataFrame sparse
    df = pd.DataFrame.sparse.from_spmatrix(matrice, index=genes, columns=barcodes)

    print(f"Data loaded: {df.shape}")
    return df, genes, barcodes, matrice

"""## Fonction pour le controle qualité"""

def QC_violin(df, genes):
    """
    Affiche des violin plots pour évaluer la qualité des cellules (QC) à partir d'une matrice d'expression.

    Paramètres :
    - df : DataFrame sparse (gènes x cellules)
    - genes : liste ou tableau des noms de gènes

    Retour :
    - genes_per_cell : nombre de gènes exprimés par cellule
    - umi_per_cell : nombre total d'UMI par cellule
    - mitochondrial_percent : pourcentage d'UMI mitochondriaux par cellule
    """

    # 1. Distribution du nombre de gènes par cellule
    print("Comptage des gènes exprimés par cellule...")
    genes_per_cell = df.astype(bool).sum(axis=0)  # Nombre de gènes exprimés par cellule
    print("Calcul du nombre de gènes exprimés par cellule terminé.")

    # 2. Distribution du nombre d'UMI par cellule
    print("Calcul du nombre d'UMI par cellule...")
    umi_per_cell = df.sum(axis=0)  # Nombre total d'UMI par cellule
    print("Calcul du nombre d'UMI par cellule terminé.")

    # 3. Distribution du nombre de gènes mitochondriaux par cellule
    print("Calcul du pourcentage de gènes mitochondriaux...")
    mitochondrial_genes = [gene for gene in genes if gene.startswith('MT-')]
    mitochondrial_counts = df.loc[mitochondrial_genes].sum(axis=0)  # Somme des gènes mitochondriaux par cellule
    mitochondrial_percent = (mitochondrial_counts / umi_per_cell) * 100
    print("Calcul du pourcentage de gènes mitochondriaux terminé.")

    # 4. Créer un DataFrame pour faciliter la visualisation avec Seaborn
    qc_data = pd.DataFrame({
        'Genes_par_cellule': genes_per_cell,
        'UMI_par_cellule': umi_per_cell,
        'Pourcentage_gènes_mitochondriaux': mitochondrial_percent
    })

    # 5. Tracer les violin plots
    plt.figure(figsize=(18, 6))

    # 5.1. Violin plot pour le nombre de gènes par cellule
    plt.subplot(1, 3, 1)
    sns.violinplot(y=qc_data['Genes_par_cellule'], color='skyblue')
    sns.stripplot(y=qc_data['Genes_par_cellule'], color='black', alpha=0.1, jitter=0.2, size=1)
    plt.title('Nombre de gènes par cellule')
    plt.ylabel('Nombre de gènes')

    # 5.2. Violin plot pour le nombre d'UMI par cellule
    plt.subplot(1, 3, 2)
    sns.violinplot(y=qc_data['UMI_par_cellule'], color='lightgreen')
    sns.stripplot(y=qc_data['UMI_par_cellule'], color='black', alpha=0.1, jitter=0.2, size=1)
    plt.title('Nombre d\'UMI par cellule')
    plt.ylabel('Nombre d\'UMI')

    # 5.3. Violin plot pour le pourcentage de gènes mitochondriaux
    plt.subplot(1, 3, 3)
    sns.violinplot(y=qc_data['Pourcentage_gènes_mitochondriaux'], color='salmon')
    sns.stripplot(y=qc_data['Pourcentage_gènes_mitochondriaux'], color='black', alpha=0.2, jitter=0.4, size=1)
    plt.title('Pourcentage de gènes mitochondriaux par cellule')
    plt.ylabel('Pourcentage de gènes mitochondriaux')

    plt.tight_layout()
    plt.show()

    return genes_per_cell, umi_per_cell, mitochondrial_percent

"""## Fonction de nettoyage"""

def filter_cells(df, genes, min_gene, max_gene, max_mito):
    """
    Filtre les cellules en fonction du nombre de gènes et du pourcentage de gènes mitochondriaux.

    Paramètres :
        df (pd.DataFrame) : DataFrame contenant les données d'expression (gènes en lignes, cellules en colonnes).
        genes (list) : Liste des noms des gènes (doit correspondre aux index de df).
        min_gene : Entier correspondant au nombre min de gènes dans une cellule, les cellules ayant moins de gène que ce nombre doivent etre supprimés
        max_gene : Entier correspondant au nombre max de gènes dans une cellule, les cellules ayant plus de gène que ce nombre doivent etre supprimés
        max_mito : Entier correspondant au pourcentage max de gènes mitochondriaux dans une cellule, les cellules ayant plus de gène mitochondriaux que
        ce nombre doivent etre supprimés
    Retourne :
        df_filtered (pd.DataFrame) : DataFrame filtré.
    """
    # 1. Calculer le nombre de gènes par cellule (somme des gènes exprimés par cellule)
    genes_per_cell = df.astype(bool).sum(axis=0)  # Nombre de gènes exprimés par cellule

    # 2. Calculer le pourcentage de gènes mitochondriaux par cellule
    mitochondrial_genes = [gene for gene in genes if gene.startswith('MT-')]  # Identifier les gènes mitochondriaux
    mitochondrial_counts = df.loc[mitochondrial_genes].sum(axis=0)  # Somme des gènes mitochondriaux par cellule
    mitochondrial_percent = (mitochondrial_counts / genes_per_cell) * 100  # Pourcentage de gènes mitochondriaux

    # 3. Appliquer les filtres
    mask = (genes_per_cell >= min_gene) & (genes_per_cell <= max_gene) & (mitochondrial_percent < max_mito)
    df_filtered = df.loc[:, mask]  # Filtrer les colonnes (cellules) qui respectent les critères

    # 4. Afficher le nombre de cellules filtrées
    print(f"Nombre de cellules avant filtrage : {df.shape[1]}")
    print(f"Nombre de cellules après filtrage : {df_filtered.shape[1]}")

    return df_filtered

"""## Fonction Normalisation"""

def normalize(df, target_sum=10_000):
    """
    On fait en sorte que la quantité d'expression des gènes dans chaque cellule soit la même et soit fixés à 10 000

    Paramètres :
        df (pd.DataFrame) : DataFrame avec les gènes en lignes et les cellules en colonnes.
        target_sum (int) : Valeur cible pour la normalisation (par défaut 10 000).

    Retourne :
        df_normalized (pd.DataFrame) : DataFrame normalisé
    """
    # 1. Normalisation par la taille de la bibliothèque
    # Calculer la somme des comptages par cellule (colonne)
    cell_sums = df.sum(axis=0)

    # Normaliser chaque cellule pour atteindre la somme cible
    df_normalized = df / cell_sums * target_sum

    return df_normalized

"""## Reduction de dimensionalité

- PCA
"""

def apply_PCA(data, n_components=20, random_state=42):
    """
    Applique une PCA (TruncatedSVD) adaptée aux données sparse après log-transformation.

    Paramètres :
    - data : matrice d'expression (cells x genes)
    - n_components : nombre de composantes principales (défaut=20)
    - random_state : graine pour la reproductibilité (défaut=42)

    Retour :
    - Matrice PCA (cells x n_components)
    """

    print("On applique la log transformation...")

    # On convertit en sparse au cas où ce n'est pas déjà le cas
    if not isinstance(data, csc_matrix):
        try:
            sparse_matrix = csc_matrix(data.sparse.to_coo())
        except AttributeError:
            sparse_matrix = csc_matrix(data)  # On convertit directement si on n'a pas .sparse
    else:
        sparse_matrix = data

    # On applique une transformation log uniquement sur les valeurs non nulles
    sparse_matrix.data = np.log1p(sparse_matrix.data)

    print("Log transformation complète.")
    print("Début de l'ACP")

    # On réduit les dimensions avec SVD car c'est plus adapté aux données sparse
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    pca_result = svd.fit_transform(sparse_matrix.T)

    # On libère la mémoire
    gc.collect()

    print(f"PCA complète, dimension = {pca_result.shape}")
    return pca_result

"""- UMAP"""

def apply_umap(pca_result, n_neighbors=15, min_dist=0.4):
    """
    Applique UMAP pour réduire la dimensionnalité après PCA.

    Paramètres :
    - pca_result : matrice PCA (cells x components)
    - n_neighbors : nombre de voisins (défaut=15)
    - min_dist : distance minimale entre les points (défaut=0.4)

    Retour :
    - Matrice UMAP (cells x 2)
    """

    umap_result = umap.UMAP(
        n_neighbors=n_neighbors,  # On peut réduire n_neighbors pour diminuer l'usage de RAM
        min_dist=min_dist,  # On peut augmenter min_dist pour diminuer la charge mémoire
        metric='euclidean',
        random_state=42
    ).fit_transform(pca_result)

    print(f"UMAP done: {umap_result.shape}")
    return umap_result

"""## Clustering

 - Kmeans
"""

def cluster_kmeans(pca_result, n_clusters):
    """
    Applique un clustering KMeans sur les données réduites (PCA).

    Paramètres :
    - pca_result : matrice de données après PCA (cells x components)
    - n_clusters : nombre de clusters à former

    Retour :
    - clusters : tableau des labels de cluster pour chaque cellule
    """

    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(pca_result)

    return clusters

"""## Visualisation

### Visualisation de projection
"""

def plot_2D_projection(data_2D, title="Projection 2D", x_label="Dim 1", y_label="Dim 2"):
    """
    Affiche une projection 2D (PCA, UMAP, etc.) sous forme de scatter plot.

    Paramètres :
    - data_2D : tableau numpy de forme (n_samples, 2)
    - title : titre du graphique
    - x_label : nom de l’axe X
    - y_label : nom de l’axe Y

    Retour :
    - Aucun (affiche le graphique)
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    plt.scatter(data_2D[:, 0], data_2D[:, 1], cmap='viridis', s=12, alpha=0.6)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.show()

"""### Visualisation d'UMAP

"""

def plot_umap(umap_result, labels):
    plt.figure(figsize=(10, 10))
    sns.scatterplot(
        x=umap_result[:, 0],
        y=umap_result[:, 1],
        hue=labels,
        palette="tab20",
        s=10,
        alpha=0.8,
        linewidth=0
    )
    plt.title("Projection de l'UMAP")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="cluster")
    plt.tight_layout()
    plt.show()

def plot_umap_with_annotations(umap_result, annotations):
    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        x=umap_result[:, 0], y=umap_result[:, 1],
        hue=annotations,
        palette="tab10", s=10, alpha=0.8, linewidth=0
    )
    plt.title("UMAP annotée selon les types cellulaires")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Type")
    plt.tight_layout()
    plt.show()

def plot_umap_cluster_labels(umap_result, clusters, cluster_annotation):
    """
    Affiche l'UMAP avec une couleur par cluster et les noms des types cellulaires au centroïde de chaque cluster,
    sauf pour ceux annotés comme "autre" ou "indéterminé".
    """
    df_umap = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
    df_umap["Cluster"] = clusters
    df_umap["Cell Type"] = df_umap["Cluster"].map(cluster_annotation)

    n_clusters = df_umap["Cluster"].nunique()
    palette = sns.color_palette("tab20", n_clusters)

    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        data=df_umap,
        x="UMAP1", y="UMAP2",
        hue="Cluster",
        palette=palette,
        s=10, alpha=0.7, linewidth=0, legend=False
    )

    # Ajouter les étiquettes sauf pour "autre"/"indéterminé"
    for cluster_id in df_umap["Cluster"].unique():
        label = cluster_annotation.get(cluster_id, f"Cluster {cluster_id}")
        if label != "Autres/Indeterminé":
            subset = df_umap[df_umap["Cluster"] == cluster_id]
            x_mean = subset["UMAP1"].mean()
            y_mean = subset["UMAP2"].mean()
            plt.text(x_mean, y_mean, label, fontsize=9, weight="bold", ha="center", va="center")

    plt.title("UMAP avec annotations par cluster (noms au centroïde)")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()
    plt.show()



def plot_gene_umap(umap_coords, expression_data, gene):
    """
    Affiche une projection UMAP colorée par l'expression du gène spécifié.

    Paramètres :
    - umap_coords : DataFrame contenant les coordonnées UMAP avec les colonnes 'UMAP1' et 'UMAP2'.
    - expression_data : DataFrame contenant les niveaux d'expression des gènes.
    - gene : Nom du gène à visualiser.
    """
    if gene not in expression_data.index:
        print(f"Le gène {gene} n'est pas présent dans les données.")
        return

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(
        umap_coords['UMAP1'], umap_coords['UMAP2'],
        c=expression_data.loc[gene], cmap='viridis', alpha=0.2, s= 10
    )
    plt.colorbar(sc, label=f'Expression de {gene}')
    plt.title(f'UMAP coloré par l\'expression de {gene}')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.show()

"""### Visualisation de Violin plots"""

def plot_violin(data, genes, clusters):
    """
    Affiche des violin plots avec jitter pour des gènes d’intérêt (ex. récepteurs IL-2), groupés par cluster.

    Paramètres :
    - data : DataFrame d'expression génique (gènes x cellules)
    - genes : liste des gènes à visualiser
    - clusters : tableau ou liste des labels de cluster par cellule

    Retour :
    - Aucun (affiche les graphiques)
    """

    # Conversion en DataFrame pour Seaborn
    expression_data = pd.DataFrame({gene: data.loc[gene] for gene in genes})
    expression_data["Cluster"] = clusters

    for gene in genes:
        plt.figure(figsize=(12, 6))

        # Violin plot
        sns.violinplot(x="Cluster", y=gene, data=expression_data, palette="muted", inner="quartile", alpha=0.8)

        # Stripplot pour ajouter du jitter
        sns.stripplot(x="Cluster", y=gene, data=expression_data, color="black", size=2, jitter=True, alpha=0.5)

        plt.title(f"Violin Plot de l'expression de {gene} par cluster avec jitter")
        plt.xlabel("Cluster")
        plt.ylabel("Expression")
        plt.show()



def plot_violin_with_significance(data, genes, clusters, p_value_table, threshold=0.05):
    """
    Affiche des violin plots avec jitter pour des gènes d’intérêt groupés par cluster,
    et des étoiles pour indiquer les différences significatives basées sur les p-values.

    Paramètres :
    - data : DataFrame d'expression génique (gènes x cellules)
    - genes : liste des gènes à visualiser
    - clusters : tableau ou liste des labels de cluster par cellule
    - p_value_table : tableau des p-values pour chaque gène et chaque cluster
    - threshold : seuil de significativité pour afficher les étoiles (par défaut 0.05)

    Retour :
    - Aucun (affiche les graphiques)
    """

    # Conversion en DataFrame pour Seaborn
    expression_data = pd.DataFrame({gene: data.loc[gene] for gene in genes})
    expression_data["Cluster"] = clusters

    for gene in genes:
        plt.figure(figsize=(12, 6))

        # Violin plot
        sns.violinplot(x="Cluster", y=gene, data=expression_data, palette="muted", inner="quartile", alpha=0.8)

        # Stripplot pour ajouter du jitter
        sns.stripplot(x="Cluster", y=gene, data=expression_data, color="black", size=2, jitter=True, alpha=0.5)

        # Ajout des étoiles pour les clusters significatifs
        clusters_unique = expression_data["Cluster"].unique()
        for cluster in clusters_unique:
            p_value = p_value_table.loc[cluster, gene]  # Récupérer la p-value pour ce gène et ce cluster

            if p_value is not None and p_value < threshold:
                # Placer l'étoile un peu au-dessus de la valeur d'expression maximale pour ce cluster
                y_max = expression_data[expression_data["Cluster"] == cluster][gene].max() * 1.1
                plt.text(cluster, y_max, '*', fontsize=16, color='red', ha='center')

        # Titres et labels
        plt.title(f"Violin Plot de l'expression de {gene} par cluster avec jitter et significativité", fontsize=16)
        plt.xlabel("Cluster", fontsize=14)
        plt.ylabel("Expression", fontsize=14)

        # Affichage du plot
        plt.tight_layout()
        plt.show()

"""### Annotation des Clusters"""

def annotate_clusters(df_normalized, clusters, umap_result):
    """
    Annoter les clusters en fonction de l'expression des marqueurs immunitaires pour identifier les sous-types cellulaires.

    Paramètres :
    - df_normalized : DataFrame d'expression normalisée (gènes x cellules)
    - clusters : Tableau numpy ou liste des clusters assignés aux cellules.
    - umap_result : Résultat de la projection UMAP (pour visualisation).

    Retour :
    - DataFrame avec l'annotation des clusters
    - Series associant chaque cellule à un type cellulaire
    """
    import numpy as np
    import pandas as pd

    marker_genes = [
        "IL2RA", "IL2RB", "IL2RG",
        "CD3E", "CD4", "CD8A",
        "CD19", "MS4A1",
        "CD14", "CD16",
        "NCAM1", "FCGR3A", "IL7R",
        "FOXP3", "TNFRSF4", "CTLA4",
        "CXCR5", "ICOS",
        "CCR7", "SELL",
        "PRF1", "GZMB",
        "XBP1", "MZB1", "CD38",
        "CLEC9A", "CD209", "ITGAX",
        "CD68", "CD163", "MRC1",
        "CD45"
    ]

    HIGH = 1.0
    MODERATE = 0.5

    cluster_annotation = {}
    unique_clusters = np.unique(clusters)

    for cluster in unique_clusters:
        print(f"\n🔎 Cluster {cluster} :")
        cluster_cells = df_normalized.columns[clusters == cluster]

        means = {}
        for gene in marker_genes:
            if gene in df_normalized.index:
                mean_expr = df_normalized.loc[gene, cluster_cells].mean()
                means[gene] = mean_expr
                print(f"{gene:12} : {mean_expr:.2f}")
            else:
                means[gene] = None
                print(f"{gene:12} : ❌ Non trouvé")

        # Annotation logique
        if means.get("CD45", 1) is not None and means["CD45"] <= 0.1:
            annotation = "Non leucocytaire"
        elif means["CD3E"] and means["CD3E"] > HIGH:
            if means["CD8A"] and means["CD8A"] > HIGH:
                if (means["GZMB"] and means["GZMB"] > HIGH) or (means["PRF1"] and means["PRF1"] > HIGH):
                    annotation = "Effector CD8+ T"
                elif means["CCR7"] and means["CCR7"] > MODERATE:
                    annotation = "Naive CD8+ T"
                else:
                    annotation = "CD8+ T"
            elif means["CD4"] and means["CD4"] > HIGH:
                if means["FOXP3"] and means["FOXP3"] > MODERATE and means["IL2RA"] and means["IL2RA"] > MODERATE:
                    if (means["TNFRSF4"] and means["TNFRSF4"] > MODERATE) or (means["CTLA4"] and means["CTLA4"] > MODERATE):
                        annotation = "Activated Tregs"
                    else:
                        annotation = "Tregs"
                elif means["CXCR5"] and means["CXCR5"] > MODERATE and means["ICOS"] and means["ICOS"] > MODERATE:
                    annotation = "T follicular helper (Tfh)"
                elif means["CCR7"] and means["CCR7"] > HIGH and means["SELL"] and means["SELL"] > MODERATE:
                    annotation = "Naive CD4+ T"
                elif (means["IL2RA"] and means["IL2RA"] > HIGH) or (means["TNFRSF4"] and means["TNFRSF4"] > HIGH):
                    annotation = "Activated CD4+ T"
                else:
                    annotation = "Memory CD4+ T"
            else:
                annotation = "T cells"
        elif (means["CD19"] and means["CD19"] > HIGH) or (means["MS4A1"] and means["MS4A1"] > HIGH):
            if (means["XBP1"] and means["XBP1"] > HIGH) or (means["MZB1"] and means["MZB1"] > HIGH):
                annotation = "Plasma cells"
            elif means["CD38"] and means["CD38"] > MODERATE:
                annotation = "Plasmablasts"
            else:
                annotation = "B cells"
        elif means["CD14"] and means["CD14"] > HIGH and (means["CD16"] is None or means["CD16"] < MODERATE):
            annotation = "Monocytes classiques"
        elif means["CD16"] and means["CD16"] > HIGH and (means["CD14"] is None or means["CD14"] < MODERATE):
            annotation = "Monocytes non classiques"
        elif means["CD14"] and means["CD14"] > MODERATE and means["CD16"] and means["CD16"] > MODERATE:
            annotation = "Monocytes intermédiaires"
        elif means["NCAM1"] and means["NCAM1"] > HIGH:
            if means["FCGR3A"] and means["FCGR3A"] > HIGH:
                annotation = "NK CD16+"
            elif means["IL7R"] and means["IL7R"] > MODERATE:
                annotation = "NK CD56bright"
            else:
                annotation = "NK cells"
        elif (means["CLEC9A"] and means["CLEC9A"] > HIGH) or \
             (means["ITGAX"] and means["ITGAX"] > HIGH) or \
             (means["CD209"] and means["CD209"] > HIGH):
            annotation = "Dendritic cells"
        elif (means["CD68"] and means["CD68"] > HIGH) or \
             (means["CD163"] and means["CD163"] > HIGH) or \
             (means["MRC1"] and means["MRC1"] > HIGH):
            annotation = "Macrophages"
        else:
            annotation = "Autres/Indeterminé"

        cluster_annotation[cluster] = annotation

    annotation_df = pd.DataFrame.from_dict(cluster_annotation, orient="index", columns=["Cell Type"])
    annotation_df.index.name = "Cluster"
    annotation_df.reset_index(inplace=True)

    display(annotation_df)

    cell_type_per_cell = pd.Series(clusters, index=df_normalized.columns).map(cluster_annotation)
    cell_type_per_cell.name = "Cell Type"

    display(cell_type_per_cell.value_counts())

    plot_umap_with_annotations(umap_result, cell_type_per_cell)

    return annotation_df, cell_type_per_cell, cluster_annotation

"""### Autres visualisations"""

def plot_expression_heatmap(data, genes, clusters):
    """
    Affiche une heatmap de l'expression moyenne des gènes d’intérêt dans chaque cluster.

    Paramètres :
    - data : DataFrame d'expression génique (gènes x cellules)
    - genes : liste des gènes à visualiser
    - clusters : tableau ou liste des labels de cluster par cellule

    Retour :
    - Aucun (affiche la heatmap)
    """

    cluster_labels = np.unique(clusters)
    mean_expression = []

    for cluster in cluster_labels:
        mean_values = [data.loc[gene, clusters == cluster].mean() for gene in genes]
        mean_expression.append(mean_values)

    mean_expression_df = pd.DataFrame(mean_expression, columns=genes, index=[f"Cluster {c}" for c in cluster_labels])

    plt.figure(figsize=(8, 6))
    sns.heatmap(mean_expression_df, annot=True, cmap="coolwarm", fmt=".2f")
    plt.title("Expression moyenne des récepteurs IL-2 par cluster")
    plt.xlabel("Gènes")
    plt.ylabel("Cluster")
    plt.show()

"""## Analyse Statistique

T-test
"""

def t_test_gene_expression_with_pvalues(df_expr, df_clusters, genes_of_interest):
    """
    Effectue un t-test pour chaque gène d'intérêt entre un cluster donné et les autres,
    puis génère un tableau des p-values pour chaque gène et cluster.

    Args:
        df_expr (pd.DataFrame): Expression des gènes (lignes = gènes, colonnes = cellules).
        df_clusters (pd.Series): Série indiquant le cluster de chaque cellule.
        genes_of_interest (list): Liste des gènes à tester.

    Returns:
        pd.DataFrame: Tableau des p-values pour chaque gène et cluster.
    """
    results = []

    clusters = df_clusters.unique()

    for gene in genes_of_interest:
        for cluster in clusters:
            cells_in_cluster = df_clusters[df_clusters == cluster].index
            expr_in_cluster = df_expr.loc[gene, cells_in_cluster]

            cells_outside_cluster = df_clusters[df_clusters != cluster].index
            expr_outside_cluster = df_expr.loc[gene, cells_outside_cluster]

            if len(expr_in_cluster) > 1 and len(expr_outside_cluster) > 1:
                t_stat, p_value = stats.ttest_ind(expr_in_cluster, expr_outside_cluster, equal_var=False)
            else:
                p_value = None

            results.append({"Gène": gene, "Cluster": cluster, "p-value": p_value})

    df_pval = pd.DataFrame(results)

    # Création du tableau des p-values avec les gènes en colonnes et les clusters en lignes
    pval_table = df_pval.pivot(index='Cluster', columns='Gène', values='p-value')

    return pval_table

"""Mann Whitney"""

def mann_whitney_gene_expression(df_expr, df_clusters, genes_of_interest):
    """
    Effectue un test de Mann-Whitney U entre l'expression d'un gène dans un cluster donné
    par rapport à tous les autres clusters.

    Args:
        df_expr (pd.DataFrame): DataFrame contenant l'expression des gènes
                                (lignes = gènes, colonnes = cellules).
        df_clusters (pd.Series): Série contenant l'affectation des cellules à un cluster.
        genes_of_interest (list): Liste des gènes à tester.

    Returns:
        pd.DataFrame: Tableau des p-values du test de Mann-Whitney U pour chaque gène et chaque cluster.
    """
    results = []

    # Liste des clusters uniques
    clusters = df_clusters.unique()

    for gene in genes_of_interest:
        for cluster in clusters:
            # Sélection des cellules appartenant au cluster cible
            cells_in_cluster = df_clusters[df_clusters == cluster].index
            expr_in_cluster = df_expr.loc[gene, cells_in_cluster]

            # Sélection des cellules appartenant aux autres clusters
            cells_outside_cluster = df_clusters[df_clusters != cluster].index
            expr_outside_cluster = df_expr.loc[gene, cells_outside_cluster]

            # Vérifier qu'il y a bien des valeurs dans les deux groupes
            if len(expr_in_cluster) > 1 and len(expr_outside_cluster) > 1:
                # Effectuer un test de Mann-Whitney U
                u_stat, p_value = stats.mannwhitneyu(expr_in_cluster, expr_outside_cluster, alternative="two-sided")
            else:
                p_value = None  # Pas assez de données pour tester

            results.append({"Gène": gene, "Cluster": cluster, "p-value": p_value})

    # Création du tableau des p-values
    p_value_df = pd.DataFrame(results)

    # Pivot pour mettre les clusters en lignes et les gènes en colonnes
    p_value_table = p_value_df.pivot(index='Cluster', columns='Gène', values='p-value')
    return p_value_table