from collections import Counter
from itertools import product

import hdbscan
import pandas as pd
import umap
from Bio import SeqIO


class HDBSCANService:
    def __init__(self, k=4):
        self.k = k
        bases = ["A", "T", "C", "G"]
        self.kmers = ["".join(p) for p in product(bases, repeat=k)]
        self.embeddings = []
        self.seq_ids = []
        self.df_emb_norm = None

    def _kmer_count(self, seq):
        seq = str(seq).upper()
        counts = Counter([seq[i : i + self.k] for i in range(len(seq) - self.k + 1)])
        return [counts.get(kmer, 0) for kmer in self.kmers]

    def embed_fasta(self, fasta_file):
        self.embeddings = []
        self.seq_ids = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            self.seq_ids.append(record.id)
            self.embeddings.append(self._kmer_count(record.seq))

        df_emb = pd.DataFrame(self.embeddings, index=self.seq_ids, columns=self.kmers)
        self.df_emb_norm = df_emb.div(df_emb.sum(axis=1), axis=0)

        return self.df_emb_norm

    def cluster_sequences(self, min_cluster_size=3, random_state=42):
        reducer = umap.UMAP(random_state=random_state)
        embedding_2d = reducer.fit_transform(self.df_emb_norm)

        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
        labels = clusterer.fit_predict(embedding_2d)

        df_clusters = pd.DataFrame({"Query_ID": self.seq_ids, "Cluster": labels})

        self.embedding_2d = embedding_2d
        self.labels = labels
        self.df_clusters = df_clusters

        return df_clusters

    def _get_kmers(self, seq, k=6):
        return [seq[i : i + k] for i in range(len(seq) - k + 1)]

    def annotate_by_kmer(self, query_seq, reference_fasta, k=6):
        references = list(SeqIO.parse(reference_fasta, "fasta"))
        query_kmers = Counter(self._get_kmers(query_seq, k))

        best_match = "No match"
        best_score = 0

        for ref in references:
            ref_kmers = Counter(self._get_kmers(str(ref.seq), k))
            shared = sum((query_kmers & ref_kmers).values())
            if shared > best_score:
                best_score = shared
                best_match = ref.id

        return best_match, best_score

    def _parse_header(self, header: str):
        hdr_words = header.split()
        species, genus = None, None

        for i in range(len(hdr_words) - 1):
            w1, w2 = hdr_words[i], hdr_words[i + 1]
            if w1[0].isupper() and w2[0].islower():
                species = f"{w1} {w2}"
                genus = w1
                break

        if not species and len(hdr_words) >= 2:
            species = f"{hdr_words[0]} {hdr_words[1]}"
            genus = hdr_words[0]

        return species, genus

    def build_reference_mapping(self, reference_fasta: str):
        ref_seqs = list(SeqIO.parse(reference_fasta, "fasta"))
        ref_info = {rec.id: rec.description for rec in ref_seqs}

        mapping = []
        for ref_id, header in ref_info.items():
            species, genus = self._parse_header(header)
            mapping.append({"ref_id": ref_id, "species": species, "genus": genus})

        df_mapping = pd.DataFrame(mapping)
        self.df_mapping = df_mapping
        return df_mapping

    def annotate_sequences(self, edna_fasta: str, reference_fasta: str, k=6):
        self.references = list(SeqIO.parse(reference_fasta, "fasta"))

        seq_annotation = {}
        for record in SeqIO.parse(edna_fasta, "fasta"):
            best_ref, score = self.annotate_by_kmer(
                str(record.seq), reference_fasta, k=k
            )
            seq_annotation[record.id] = best_ref

        self.seq_annotation = seq_annotation
        return seq_annotation

    def summarize_clusters(
        self, labels, seq_ids, mapping_df, output_csv="cluster_annotation_summary.csv"
    ):
        ref_to_species = {
            row["ref_id"]: row["species"] for _, row in mapping_df.iterrows()
        }
        ref_to_genus = {row["ref_id"]: row["genus"] for _, row in mapping_df.iterrows()}

        cluster_summary = []

        for cluster_id in set(labels):
            cluster_seq_ids = [
                seq_id for seq_id, lbl in zip(seq_ids, labels) if lbl == cluster_id
            ]
            num_sequences = len(cluster_seq_ids)

            if cluster_id == -1:
                cluster_summary.append(
                    {
                        "Cluster": cluster_id,
                        "Assigned_Reference": "Unknown",
                        "Assigned_Species": "Unknown",
                        "Assigned_Genus": "Unknown",
                        "Num_Sequences": num_sequences,
                        "Status": "Novel/Noise",
                    }
                )
                continue

            ref_counts = Counter(
                [self.seq_annotation.get(seq, "No match") for seq in cluster_seq_ids]
            )

            if not ref_counts:
                cluster_summary.append(
                    {
                        "Cluster": cluster_id,
                        "Assigned_Reference": "Unknown",
                        "Assigned_Species": "Unknown",
                        "Assigned_Genus": "Unknown",
                        "Num_Sequences": num_sequences,
                        "Status": "Novel",
                    }
                )
                continue

            most_common_ref, count = ref_counts.most_common(1)[0]

            if most_common_ref == "No match":
                cluster_summary.append(
                    {
                        "Cluster": cluster_id,
                        "Assigned_Reference": "Unknown",
                        "Assigned_Species": "Unknown",
                        "Assigned_Genus": "Unknown",
                        "Num_Sequences": num_sequences,
                        "Status": "Novel",
                    }
                )
            else:
                species = ref_to_species.get(most_common_ref, "Unknown")
                genus = ref_to_genus.get(most_common_ref, "Unknown")
                cluster_summary.append(
                    {
                        "Cluster": cluster_id,
                        "Assigned_Reference": most_common_ref,
                        "Assigned_Species": species,
                        "Assigned_Genus": genus,
                        "Num_Sequences": num_sequences,
                        "Status": "Known",
                    }
                )

        df_cluster_summary = pd.DataFrame(cluster_summary)
        df_cluster_summary.to_csv(output_csv, index=False)

        self.df_cluster_summary = df_cluster_summary
        return df_cluster_summary

    def build_sequence_taxonomy(self, df_mapping, df_clusters, seq_annotation):
        ref_taxonomy = {}
        for _, row in df_mapping.iterrows():
            ref_taxonomy[row["ref_id"]] = {
                "species": row.get("species", "Unknown"),
                "genus": row.get("genus", "Unknown"),
                "family": row.get("family", "Unknown"),
            }

        sequence_taxonomy = []
        for seq_id, ref_id in seq_annotation.items():
            if ref_id in ref_taxonomy:
                species = ref_taxonomy[ref_id]["species"]
                genus = ref_taxonomy[ref_id]["genus"]
                family = ref_taxonomy[ref_id]["family"]
            else:
                species = genus = family = "Unknown"

            sequence_taxonomy.append(
                {
                    "Sequence_ID": seq_id,
                    "Reference_ID": ref_id,
                    "Species": species,
                    "Genus": genus,
                    "Family": family,
                    "Cluster": df_clusters.loc[
                        df_clusters["Query_ID"] == seq_id, "Cluster"
                    ].values[0],
                }
            )

        df_seq_taxonomy = pd.DataFrame(sequence_taxonomy)
        self.df_seq_taxonomy = df_seq_taxonomy
        return df_seq_taxonomy

    def compute_abundance(self, df_seq_taxonomy=pd.DataFrame()):
        if not df_seq_taxonomy:
            df_seq_taxonomy = self.df_seq_taxonomy
        species_counts = df_seq_taxonomy["Species"].value_counts().reset_index()
        species_counts.columns = ["Species", "Abundance"]

        genus_counts = df_seq_taxonomy["Genus"].value_counts().reset_index()
        genus_counts.columns = ["Genus", "Abundance"]

        cluster_counts = df_seq_taxonomy["Cluster"].value_counts().reset_index()
        cluster_counts.columns = ["Cluster", "Num_Sequences"]

        self.species_counts = species_counts
        self.genus_counts = genus_counts
        self.cluster_counts = cluster_counts

        return species_counts, genus_counts, cluster_counts

    def export_results(self, output_dir="."):
        """
        Save annotation and abundance tables as CSVs.
        """
        if not hasattr(self, "df_seq_taxonomy"):
            raise ValueError("Run annotate_sequences_with_clusters first.")
        if not hasattr(self, "species_counts"):
            raise ValueError("Run compute_abundance first.")

        self.df_seq_taxonomy.to_csv(
            f"{output_dir}/sequence_taxonomy_annotation.csv", index=False
        )
        self.species_counts.to_csv(f"{output_dir}/species_abundance.csv", index=False)
        self.genus_counts.to_csv(f"{output_dir}/genus_abundance.csv", index=False)
        self.cluster_counts.to_csv(f"{output_dir}/cluster_summary.csv", index=False)

        print(f"Results saved in {output_dir}")

    def evaluate_accuracy(self, ground_truth_csv: str):
        """
        Compare predicted taxonomy with ground-truth CSV and return species and genus accuracy.
        """
        if not hasattr(self, "df_seq_taxonomy"):
            raise ValueError("Run build_sequence_taxonomy first.")

        ground_truth = pd.read_csv(ground_truth_csv)

        df_compare = pd.merge(
            self.df_seq_taxonomy,
            ground_truth,
            left_on="Sequence_ID",
            right_on="Sequence_ID",
            suffixes=("_pred", "_true"),
        )

        species_correct = (
            df_compare["Species_pred"] == df_compare["Species_true"]
        ).sum()
        species_total = len(df_compare)
        species_accuracy = species_correct / species_total * 100

        genus_correct = (df_compare["Genus_pred"] == df_compare["Genus_true"]).sum()
        genus_total = len(df_compare)
        genus_accuracy = genus_correct / genus_total * 100

        self.species_accuracy = species_accuracy
        self.genus_accuracy = genus_accuracy

        print(f"Species-level accuracy: {species_accuracy:.2f}%")
        print(f"Genus-level accuracy: {genus_accuracy:.2f}%")

        return species_accuracy, genus_accuracy
