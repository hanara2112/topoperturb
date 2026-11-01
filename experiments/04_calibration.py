"""Step 2c -- does each test keep its promised false-positive rate?

A test at p < 0.05 should reject about 5% of the time when nothing differs. We
check this on real cells by comparing controls with controls: for each cell line
and several group sizes we draw many control-vs-control pairs and record how often
each method rejects, with a Wilson interval, and whether its p-values look uniform
(a KS test against Uniform(0,1)). Methods whose rate sits well above 5% are not
trustworthy as evidence, however strong they look on the real knockdowns.

    python experiments/04_calibration.py              # 200 draws, sizes 75/100/150
    python experiments/04_calibration.py --draws 50
"""
import argparse
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kstest

import config
from topoperturb import data, plotting, twosample

CONTROL = data.CONTROL
SEED = config.SEED
SIZES = [75, 100, 150]
CLOUDS = {cl: data.load_clouds(config.cloud_path(cl))[0] for cl in ("k562", "rpe1")}


def run(task):
    """One control-vs-control draw at a given size: all methods' p-values."""
    cell_line, size, draw = task
    Z, labels = CLOUDS[cell_line][0], CLOUDS[cell_line][1]
    control_idx = np.where(labels == CONTROL)[0]
    rng = np.random.default_rng([SEED, 13, SIZES.index(size), draw])
    pick = rng.choice(control_idx, 2 * size, replace=False)
    row = {"cell_line": cell_line, "size": size, "draw": draw}
    for method, (p, _obs) in twosample.permutation_test(
            twosample.all_statistics, Z[pick[:size]], Z[pick[size:]], rng=rng).items():
        row[f"p_{method}"] = p
    return row


def summarise(draws):
    """False-positive rate, Wilson interval and p-value uniformity per cell/size/method."""
    rows = []
    for (cell_line, size), subset in draws.groupby(["cell_line", "size"]):
        for method in twosample.METHODS:
            p = subset[f"p_{method}"].to_numpy()
            k = int((p < 0.05).sum())
            lo, hi = twosample.wilson_interval(k, len(p))
            rows.append({"cell_line": cell_line, "size": size, "method": method,
                         "false_positive_rate": k / len(p), "ci_low": lo, "ci_high": hi,
                         "n_draws": len(p), "uniform_ks_p": float(kstest(p, "uniform").pvalue)})
    return pd.DataFrame(rows)


def main(draws_per_cell, out):
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(cl, size, d) for cl in CLOUDS for size in SIZES for d in range(draws_per_cell)]
    with Pool(8) as pool:
        draws = pd.DataFrame(pool.map(run, tasks, chunksize=1))
    draws.to_csv(out / "calibration_draws.csv", index=False)
    table = summarise(draws)
    table.to_csv(out / "calibration.csv", index=False)

    fig, axes = plotting.plt.subplots(1, len(CLOUDS), figsize=(5 * len(CLOUDS), 3.6),
                                      squeeze=False, sharey=True)
    for ax, cell_line in zip(axes[0], CLOUDS):
        for method in twosample.METHODS:
            sub = table[(table.cell_line == cell_line) & (table.method == method)].sort_values("size")
            plotting.calibration_curve(ax, sub["size"], sub.false_positive_rate,
                                       sub.ci_low, sub.ci_high, method)
        ax.axhline(0.05, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(f"{cell_line}: control-vs-control false-positive rate", fontsize=9)
    axes[0][-1].legend(fontsize=6, loc="upper right")
    plotting.save(fig, out / "calibration.png")
    print(table.pivot_table(index=["cell_line", "size"], columns="method",
                            values="false_positive_rate").round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--draws", type=int, default=200, help="control-vs-control pairs per size")
    parser.add_argument("--output", type=Path, default=config.RESULTS)
    args = parser.parse_args()
    t = time.time()
    main(args.draws, args.output)
    print(f"{time.time() - t:.0f}s")
