#!/usr/bin/env python3
r"""
Apply QC filters to a single sample and generate QC plots.

This script:
- Loads a pickled cisTopic object
- Runs LDA models with user-defined parameters
- Evaluates models and saves evaluation plot
- Writes updated cisTopic object back to disk

Usage:
python CreateATACCountMatrix.py \
    --outdir outs/qc_output \
    --in_cistopic_obj input.pkl \
    --out_cistopic_obj output.pkl \
    --n_cpu 8 \
    --n_topics 2 5 10 15 20 \
    --n_iter 500 \
    --random_state 555 \
    --alpha 50 \
    --alpha_by_topic \
    --eta 0.1 \
    --eta_by_topic
"""

import argparse
from pycisTopic.lda_models import run_cgs_models, evaluate_models
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg") 
import pickle


def main():
    p = argparse.ArgumentParser(description="Apply QC filters and build cisTopic object.")
    p.add_argument("--outdir", type=str, required=True, help="Output directory for QC results.")
    p.add_argument("--in_cistopic_obj", required=True, help="Path to input pickled cistopic object.")
    p.add_argument("--out_cistopic_obj", required=True, help="Path to output pickled cistopic object.")
    p.add_argument("--n_cpu", type=int, default=1, help="Number of CPUs to use.")

    # run_cgs_models arguments
    p.add_argument("--n_topics", type=int, nargs="+", default=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
                   help="List of topic numbers to try (space separated).")
    p.add_argument("--n_iter", type=int, default=500, help="Number of iterations for each model.")
    p.add_argument("--random_state", type=int, default=555, help="Random seed.")
    p.add_argument("--alpha", type=float, default=50, help="Alpha hyperparameter.")
    p.add_argument("--alpha_by_topic", action="store_true", help="If set, alpha is scaled by number of topics.")
    p.add_argument("--eta", type=float, default=0.1, help="Eta hyperparameter.")
    p.add_argument("--eta_by_topic", action="store_true", help="If set, eta is scaled by number of topics.")

    # eval m

    args = p.parse_args()

    # Load input cisTopic object
    with open(args.in_cistopic_obj, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Run models
    models = run_cgs_models(
        cistopic_obj,
        n_topics=args.n_topics,
        n_cpu=args.n_cpu,
        n_iter=args.n_iter,
        random_state=args.random_state,
        alpha=args.alpha,
        alpha_by_topic=args.alpha_by_topic,
        eta=args.eta,
        eta_by_topic=args.eta_by_topic,
        save_path=f"{args.outdir}/Topics"
    )

    # Evaluate and save plot
    model = evaluate_models(
        models,
        select_model=None,
        return_model=True,
        plot=True
    )

    fig = plt.gcf()
    fig.savefig(f"{args.outdir}/Plots/model_evaluation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
    cistopic_obj.add_LDA_model(model)

    # Write updated object
    with open(args.out_cistopic_obj, "wb") as f:
        pickle.dump(cistopic_obj, f)


if __name__ == "__main__":
    main()
