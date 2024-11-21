import argparse
from evolvepro.src.evolve import grid_search

def create_parser():
    parser = argparse.ArgumentParser(description="Run experiments with different combinations of grid search variables.")
    parser.add_argument("--dataset_name", type=str, help="Name of the dataset")
    parser.add_argument("--experiment_name", type=str, help="Name of the experiment")
    parser.add_argument("--model_name", type=str, help="Name of the model used for embeddings")
    parser.add_argument("--embeddings_path", type=str, help="Path to the embeddings file")
    parser.add_argument("--labels_path", type=str, help="Path to the labels file")
    parser.add_argument("--num_simulations", type=int, help="Number of simulations for each parameter combination. Example: 3, 10")
    parser.add_argument("--num_iterations", type=int, nargs="+", help="List of number of iterations. Example: 3 5 10 (must be greater than 1)")
    parser.add_argument("--measured_var", type=str, nargs="+", help="Fitness type to train on. Options: activity activity_scaled")
    parser.add_argument("--learning_strategies", type=str, nargs="+", help="Type of learning strategy. Options: random top5bottom5 top10 dist")
    parser.add_argument("--num_mutants_per_round", type=int, nargs="+", help="Number of mutants per round. Example: 8 10 16 32 128")
    parser.add_argument("--num_final_round_mutants", type=int, help="Number of mutants in final round. Example: 16")
    parser.add_argument("--first_round_strategies", type=str, nargs="+", help="Type of first round strategy. Options: random diverse_medoids representative_hie")
    parser.add_argument("--embedding_types", type=str, nargs="+", help="Types of embeddings to train on. Options: embeddings embeddings_pca")
    parser.add_argument("--pca_components", type=int, nargs="*", default=None, help="Number of PCA components to use")
    parser.add_argument("--regression_types", type=str, nargs="+", help="Regression types. Options: ridge lasso elasticnet linear neuralnet randomforest gradientboosting knn gp")
    parser.add_argument("--embeddings_file_type", type=str, help="Type of embeddings file to read. Options: csv pt")
    parser.add_argument("--output_dir", type=str, help="output directory for grid search results")
    parser.add_argument("--embeddings_type_pt", type=str, default=None, help="Type of embeddings to use (for .pt files). Options: average mutated both")
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    grid_search(
        dataset_name=args.dataset_name,
        experiment_name=args.experiment_name,
        model_name=args.model_name,
        embeddings_path=args.embeddings_path,
        labels_path=args.labels_path,
        num_simulations=args.num_simulations,
        num_iterations=args.num_iterations,
        measured_var=args.measured_var,
        learning_strategies=args.learning_strategies,
        num_mutants_per_round=args.num_mutants_per_round,
        num_final_round_mutants=args.num_final_round_mutants,
        first_round_strategies=args.first_round_strategies,
        embedding_types=args.embedding_types,
        pca_components=args.pca_components,
        regression_types=args.regression_types,
        embeddings_file_type=args.embeddings_file_type,
        output_dir=args.output_dir,
        embeddings_type_pt=args.embeddings_type_pt
    )

if __name__ == "__main__":
    main()