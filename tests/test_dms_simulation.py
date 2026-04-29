from __future__ import annotations

from evolvepro.src.evolve import directed_evolution_simulation


def test_directed_evolution_simulation_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    variants = [f"variant_{i}" for i in range(80)] + ["WT"]
    embeddings = synthetic_embeddings(n_variants=80, n_dim=12, seed=7)
    labels = synthetic_labels(variants, seed=7)

    result = directed_evolution_simulation(
        labels=labels,
        embeddings=embeddings,
        num_simulations=1,
        num_iterations=2,
        num_mutants_per_round=8,
        measured_var="activity",
        regression_type="ridge",
        learning_strategy="topn",
        final_round=8,
        first_round_strategy="random",
        embedding_type=None,
    )

    data_regression.check(
        {
            "shape": list(result.shape),
            "columns": result.columns.tolist(),
            "records": result.round(6).reset_index(drop=True).to_dict(orient="records"),
        }
    )
