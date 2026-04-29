from __future__ import annotations

import math

import pandas as pd

from evolvepro.src.model import first_round, top_layer


def _sorted_non_wt(iteration_zero: pd.DataFrame) -> list[str]:
    return sorted(
        iteration_zero.loc[iteration_zero["variant"] != "WT", "variant"].tolist()
    )


def test_first_round_random_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    variants = [f"variant_{i}" for i in range(100)] + ["WT"]
    labels = synthetic_labels(variants, seed=42)
    embeddings = synthetic_embeddings(n_variants=100, n_dim=8, seed=42)

    labels_zero, iteration_zero, this_round_variants = first_round(
        labels=labels,
        embeddings=embeddings,
        num_mutants_per_round=10,
        first_round_strategy="random",
        random_seed=42,
    )

    data_regression.check(
        {
            "labels_shape": list(labels_zero.shape),
            "iteration_zero": iteration_zero.to_dict(orient="records"),
            "selected_sorted": _sorted_non_wt(iteration_zero),
            "this_round_variants": this_round_variants.tolist(),
        }
    )


def test_first_round_diverse_medoids_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    variants = [f"variant_{i}" for i in range(50)] + ["WT"]
    labels = synthetic_labels(variants, seed=42)
    embeddings = synthetic_embeddings(n_variants=50, n_dim=16, seed=42)

    labels_zero, iteration_zero, this_round_variants = first_round(
        labels=labels,
        embeddings=embeddings,
        num_mutants_per_round=8,
        first_round_strategy="diverse_medoids",
        random_seed=42,
    )

    data_regression.check(
        {
            "labels_shape": list(labels_zero.shape),
            "iteration_zero": iteration_zero.to_dict(orient="records"),
            "selected_sorted": _sorted_non_wt(iteration_zero),
            "this_round_variants": this_round_variants.tolist(),
        }
    )


def _top_layer_case(synthetic_embeddings, synthetic_labels):
    variants = [f"variant_{i}" for i in range(40)]
    embeddings = synthetic_embeddings(n_variants=40, n_dim=8, seed=0, include_wt=False)
    labels = synthetic_labels(variants, seed=0)
    labels["iteration"] = [1.0] * 20 + [math.nan] * 20
    return embeddings, labels


def _top_layer_payload(result, decimals: int) -> dict:
    (
        train_error,
        test_error,
        train_r_squared,
        test_r_squared,
        alpha,
        median_activity_scaled,
        top_activity_scaled,
        top_variant,
        top_final_round_variants,
        activity_binary_percentage,
        spearman_corr,
        df_test,
        this_round_variants,
    ) = result
    return {
        "metrics": {
            "train_error": round(float(train_error), decimals),
            "test_error": round(float(test_error), decimals),
            "train_r_squared": round(float(train_r_squared), decimals),
            "test_r_squared": round(float(test_r_squared), decimals),
            "alpha": round(float(alpha), decimals),
            "median_activity_scaled": round(float(median_activity_scaled), decimals),
            "top_activity_scaled": round(float(top_activity_scaled), decimals),
            "top_variant": top_variant,
            "top_final_round_variants": top_final_round_variants,
            "activity_binary_percentage": round(
                float(activity_binary_percentage), decimals
            ),
            "spearman_corr": round(float(spearman_corr), decimals),
        },
        "df_test": df_test.round(decimals).to_dict(orient="records"),
        "this_round_variants": this_round_variants.tolist(),
    }


def _run_top_layer(synthetic_embeddings, synthetic_labels, regression_type: str):
    embeddings, labels = _top_layer_case(synthetic_embeddings, synthetic_labels)
    return top_layer(
        iter_train=[1],
        iter_test=None,
        embeddings_pd=embeddings,
        labels_pd=labels,
        measured_var="activity",
        regression_type=regression_type,
        final_round=10,
    )


def test_top_layer_ridge_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    result = _run_top_layer(synthetic_embeddings, synthetic_labels, "ridge")
    data_regression.check(_top_layer_payload(result, decimals=6))


def test_top_layer_randomforest_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    result = _run_top_layer(synthetic_embeddings, synthetic_labels, "randomforest")
    data_regression.check(_top_layer_payload(result, decimals=4))


def test_top_layer_knn_synthetic(
    synthetic_embeddings, synthetic_labels, data_regression
):
    result = _run_top_layer(synthetic_embeddings, synthetic_labels, "knn")
    data_regression.check(_top_layer_payload(result, decimals=6))
