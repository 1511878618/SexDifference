#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description: SimilarityMatcher with batch encoding and accelerated similarity computation
@Date       : 2025/01/08
@Author     : Tingfeng Xu
@Version    : 2.0
"""
import argparse
import os
import pickle
import torch
from transformers import BertTokenizer, BertModel
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm


class SimilarityMatcher:
    def __init__(self, model_name="bert-base-uncased", model_dir=None, device=None):
        """
        Initialize BERT tokenizer and model.
        Args:
            model_name (str): Pretrained model name from HuggingFace.
            model_dir (str): Directory path for local model weights.
            device (str): Device to run the model on ('cuda' or 'cpu').
        """
        if model_dir:
            self.tokenizer = BertTokenizer.from_pretrained(model_dir)
            self.model = BertModel.from_pretrained(model_dir)
        else:
            self.tokenizer = BertTokenizer.from_pretrained(model_name)
            self.model = BertModel.from_pretrained(model_name)

        self.device = (
            device if device else ("cuda" if torch.cuda.is_available() else "cpu")
        )
        self.model.to(self.device)
        self.model.eval()

    def _encode_batch(self, phrases):
        """
        Encode a batch of phrases into vectors using BERT.
        Args:
            phrases (list): List of input phrases to encode.
        Returns:
            numpy.ndarray: Encoded vector representations.
        """
        inputs = self.tokenizer(
            phrases, return_tensors="pt", truncation=True, padding=True, max_length=32
        )
        inputs = {key: val.to(self.device) for key, val in inputs.items()}
        with torch.no_grad():
            outputs = self.model(**inputs)
            cls_embeddings = outputs.last_hidden_state[
                :, 0, :
            ]  # Extract [CLS] token embeddings
        return cls_embeddings.cpu().numpy()

    def encode_and_save_database(self, database_file, save_path, batch_size=32):
        """
        Encode phrases from a database file in batches and save them as a .pkl file.
        Args:
            database_file (str): Path to the input database file (CSV with one column).
            save_path (str): Path to save the encoded database as a .pkl file.
            batch_size (int): Number of phrases to encode in a single batch.
        """
        df = pd.read_csv(database_file, header=None)
        phrases = df[0].tolist()

        database_vectors = {}
        print(f"Encoding database in batches of size {batch_size}...")
        for i in tqdm(range(0, len(phrases), batch_size), desc="Encoding"):
            batch_phrases = phrases[i : i + batch_size]
            batch_vectors = self._encode_batch(batch_phrases)
            for phrase, vector in zip(batch_phrases, batch_vectors):
                database_vectors[phrase] = vector

        with open(save_path, "wb") as file:
            pickle.dump(database_vectors, file)
        print(f"Database vectors saved to {save_path}")

    def load_encoded_database(self, file_path):
        """
        Load encoded database vectors from a .pkl file.
        Args:
            file_path (str): Path to the .pkl file.
        Returns:
            dict: Loaded dictionary of phrases and their encoded vectors.
        """
        with open(file_path, "rb") as file:
            database_vectors = pickle.load(file)
        return database_vectors

    def compute_batch_similarity(self, input_vectors, database_vectors, top_k=5):
        """
        Compute cosine similarity between input vectors and database vectors in batch.
        Args:
            input_vectors (numpy.ndarray): Encoded vectors for the input phrases.
            database_vectors (dict): Dictionary of phrases and their encoded vectors.
        Returns:
            list: List of similarity results for each input phrase.
        """
        db_phrases = list(database_vectors.keys())
        db_vectors = np.array(list(database_vectors.values()))

        results = []
        print("Computing similarities in batch...")
        for input_vector in tqdm(input_vectors, desc="Similarity Computation"):
            similarities = cosine_similarity([input_vector], db_vectors)[
                0
            ]  # Compute cosine similarities
            top_indices = similarities.argsort()[::-1][:top_k]  # Get top K indices
            top_matches = [(db_phrases[i], similarities[i]) for i in top_indices]
            results.append(top_matches)
        return results


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="SimilarityMatcher: Match input phrases with the most similar phrases in the database."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Path to a file containing input phrases (one phrase per line).",
    )
    parser.add_argument(
        "--str", help="Input a single phrase directly for similarity matching."
    )
    parser.add_argument(
        "-d",
        "--database",
        required=True,
        nargs="+",
        help="Database files (.pkl or single-column CSV).",
    )
    parser.add_argument(
        "--out", required=False, help="Path to save the output results as a CSV file."
    )
    parser.add_argument(
        "-k",
        "--top_k",
        type=int,
        default=5,
        help="Number of top similar phrases to return (default: 5).",
    )
    parser.add_argument(
        "-m",
        "--model",
        default="bert-base-uncased",
        help="Pretrained model name to use (default: bert-base-uncased).",
    )
    parser.add_argument(
        "--model-dir", help="Path to local model weights, only used with -m."
    )
    parser.add_argument(
        "--preprocess",
        action="store_true",
        help="Preprocess the database file and save encoded .pkl.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for encoding (default: 32).",
    )
    parser.add_argument(
        "--gpu", action="store_true", help="Force the model to use GPU for computation."
    )

    args = parser.parse_args()

    # Determine device (GPU or CPU)
    device = "cuda" if args.gpu and torch.cuda.is_available() else "cpu"

    # Check input arguments
    if not args.input and not args.str:
        raise ValueError("Must provide either -i/--input or --str argument.")
    if args.input and args.str:
        raise ValueError("Cannot provide both -i/--input and --str arguments.")

    # input_phrases = []
    if args.input:
        input_phrases = pd.read_csv(args.input, header=None)[0].tolist()

    elif args.str:
        input_phrases = [args.str]

    if not args.out:
        raise ValueError("Must provide --out argument to save the output results.")

    # Initialize the matcher
    matcher = SimilarityMatcher(
        model_name=args.model, model_dir=args.model_dir, device=device
    )

    # Load and process the database
    database_vectors = {}
    for db_file in args.database:
        if db_file.endswith(".pkl"):
            database_vectors.update(matcher.load_encoded_database(db_file))
        else:
            pkl_file = Path(db_file).with_suffix(".pkl")
            if not os.path.exists(pkl_file):
                print(f"Encoding database file {db_file} into vectors...")
                matcher.encode_and_save_database(
                    db_file, pkl_file, batch_size=args.batch_size
                )
            database_vectors.update(matcher.load_encoded_database(pkl_file))

    if args.preprocess:
        print("Preprocessing complete. Please rerun the script without --preprocess.")
        return

    # Encode input phrases in batch
    print("Encoding input phrases...")
    input_vectors = matcher._encode_batch(input_phrases)

    # Compute similarities
    similarity_results = matcher.compute_batch_similarity(
        input_vectors, database_vectors, top_k=args.top_k
    )

    # Save results to CSV
    print("Saving results to CSV...")
    results = []
    for phrase, matches in zip(input_phrases, similarity_results):
        row = [phrase]
        for match, score in matches:
            row.append(match)
            row.append(score)
        results.append(row)

    output_df = pd.DataFrame(results)
    columns = ["InputPhrase"]
    for i in range(1, args.top_k + 1):
        columns.append(f"Match{i}")
        columns.append(f"Similarity{i}")

    output_df.columns = columns
    output_df.to_csv(args.out, index=False, encoding="utf-8")
    print(f"Results saved to {args.out}")


if __name__ == "__main__":
    main()
