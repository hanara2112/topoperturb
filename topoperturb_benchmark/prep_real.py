"""Build the 30-PC clouds used by the benchmark (run once).

K562: controls + the 8 pilot knockdowns, PCA fit on those cells (as in Phases 1-2).
RPE1: controls + ZMAT2, CPSF3, POLRMT, GFM1 + positive control, own PCA (as in Phase 2.4).
Also stores K562 control-covariance summaries that parameterise the synthetic data.
Reads the cached matrices from earlier phases; falls back to the Hugging Face mirrors.
"""
import gzip, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist

HERE = Path(__file__).resolve().parent
BIO = HERE.parent
SCRATCH = Path("/private/tmp/claude-501/-Users-tanisha-Desktop-bio/a7af2154-21a5-4987-9e41-41e24be1e4b5/scratchpad")
PILOT = pd.read_csv(BIO / "topoperturb_pilot_outputs/pilot_perturbations.csv")
RPE1_SHARED = {"ENSG00000146007": "ZMAT2", "ENSG00000119203": "CPSF3",
               "ENSG00000099821": "POLRMT", "ENSG00000168827": "GFM1"}
N_PCS, N_MIN = 30, 150


def read_hf(repo, fname, keep=None):
    from huggingface_hub import hf_hub_download
    path = hf_hub_download(repo_id=repo, filename=fname, repo_type="dataset")
    with gzip.open(path, "rt") as f:
        cols = f.readline().rstrip("\n\r").split("\t")[:-1]
    dt = {g: np.float32 for g in cols}
    dt["intervention"] = "string"
    parts = []
    for ch in pd.read_csv(path, sep="\t", compression="gzip", dtype=dt, chunksize=10000):
        m = ch["intervention"] != "excluded"
        if keep is not None:
            m &= ch["intervention"].isin(keep)
        parts.append(ch.loc[m])
    d = pd.concat(parts, ignore_index=True)
    return d[cols].to_numpy(np.float32), d["intervention"].astype(str).to_numpy()


def load(cache, repo, fname, keep=None):
    if cache.exists():
        z = np.load(cache)
        return z["X"], z["lab"].astype(str)
    return read_hf(repo, fname, keep)


def pca30(X):
    return PCA(n_components=N_PCS, svd_solver="randomized", random_state=42).fit_transform(X).astype(np.float64)


def pr(C):
    return float(np.trace(C) ** 2 / np.trace(C @ C))


if __name__ == "__main__":
    Xk, lk = load(SCRATCH / "k562_pilot_cache.npz", "huggingbahl21/causalbench-k562",
                  "real/causalbench/data/causalbench_k562.tsv.gz", set(PILOT.gene_id) | {"non-targeting"})
    Zk = pca30(Xk)

    Xr, lr = load(SCRATCH / "rpe1_cache.npz", "huggingbahl21/causalbench-rpe1",
                  "real/causalbench/data/causalbench_rpe1.tsv.gz")
    cm = lr == "non-targeting"
    mu = Xr[cm].mean(0)
    counts = pd.Series(lr[~cm]).value_counts()
    shifts = {g: np.linalg.norm(Xr[lr == g].mean(0) - mu) for g, n in counts.items()
              if n >= N_MIN and g not in RPE1_SHARED}
    pos = max(shifts, key=shifts.get)
    keep = np.isin(lr, list(RPE1_SHARED) + [pos, "non-targeting"])
    Zr, lr = pca30(Xr[keep]), lr[keep]

    # K562 control geometry -> synthetic parameters
    ck = Zk[lk == "non-targeting"]
    Cc = np.cov(ck.T)
    lam = np.sort(np.linalg.eigvalsh(Cc))[::-1]
    rad = lambda P: pdist(P).mean()
    genes = []
    for g, s in zip(PILOT.gene_id, PILOT.gene_symbol):
        P = Zk[lk == g]
        genes.append({"gene": s, "centroid_shift": float(np.linalg.norm(P.mean(0) - ck.mean(0))),
                      "size_ratio": float(rad(P[:2000]) / rad(ck[:2000])), "PR": pr(np.cov(P.T))})
    meta = {"lam": lam.tolist(), "PR_ctrl": pr(Cc),
            "top2_var_frac": float(lam[:2].sum() / lam.sum()),
            "median_centroid_shift": float(np.median([g["centroid_shift"] for g in genes])),
            "k562_genes": genes,
            "k562_symbols": dict(zip(PILOT.gene_id, PILOT.gene_symbol)),
            "rpe1_symbols": {**RPE1_SHARED, pos: f"pos. control ({pos})"}}
    np.savez(HERE / "outputs/real_clouds.npz", Zk=Zk, lk=lk.astype("U32"), Zr=Zr, lr=lr.astype("U32"))
    (HERE / "outputs/real_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if k != "lam"}, indent=1))
