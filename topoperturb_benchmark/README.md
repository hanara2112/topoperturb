# TopoPerturb benchmark: geometry vs topology

Synthetic clouds with known truth plus real CausalBench knockdowns (K562, RPE1),
all run through the same battery of 10 tests (150 vs 150 cells, 30 PCs).

## Run

```bash
python prep_real.py            # builds outputs/real_clouds.npz + real_meta.json (K562/RPE1 30-PC clouds)
python synth_bench.py 40       # 12 synthetic conditions x 40 replicates  (~35 min on 8 cores)
python real_battery.py 15      # 15 real groups x 15 draws                 (~20 min on 8 cores)
python make_figures.py         # outputs/benchmark_heatmap.png
```

Stack: NumPy, SciPy, pandas, scikit-learn, ripser + persim (scikit-tda), matplotlib.

## Files

| File | What it does |
|---|---|
| `topo_battery.py` | The 10 tests: energy, covariance diff., raw/scale-normalised H1, scale-normalised H0, effective-resistance H1, diffusion H1 (199 label permutations each), and three covariance-matched Gaussian surrogate tests (99 surrogate pairs) |
| `synth_bench.py` | Synthetic conditions built on the K562 control covariance spectrum |
| `real_battery.py` | Same battery on 8 K562 + 4 RPE1 knockdowns, an RPE1 positive control and control-vs-control nulls |
| `outputs/synthetic_detection_table.csv` | Detection rate (p < 0.05) per condition x method, Wilson 95 % CI |
| `outputs/real_detection_table.csv` | Fraction of draws with p < 0.05, primary-draw p and BH q per gene x method |
| `outputs/*_replicates.csv`, `real_draws.csv` | Every p-value |

## Results (detection rate, % of 40 synthetic replicates)

| Condition | Raw H1 | Scale-norm. H1 | Diffusion H1 | Surrogate-adj. H1 |
|---|---|---|---|---|
| Null | 5 | 2 | 8 | 0 |
| Translation | 5 | 2 | 5 | 0 |
| Scaling x1.3 | **100** | 25 | 2 | 0 |
| Anisotropy (PR 16.6 to 10) | 68 | **75** | 2 | 8 |
| Heavy tails | **100** | **100** | 10 | **95** |
| Partial knockdown | 22 | 8 | 15 | 5 |
| Loop in top-2 PCs (K562-like, 26 % var.) | 2 | 10 | **68** | 5 |
| Two clusters on PC1 | 2 | 5 | 30 | 5 |
| Loop, 49 % var. | 5 | 12 | **100** | 0 |
| Loop, 32 % var. | 2 | 5 | **70** | 2 |
| Loop, 18 % var. | 8 | 8 | 5 | 8 |

1. **Raw H1 measures geometry.** It fires on every scaling and most anisotropy
   replicates and misses a covariance-matched loop at K562-like noise.
2. **Scale normalisation removes translation and most of scaling, not anisotropy.**
   This matches the POLR2B / RPL21 residuals in Phase 2.
3. **Euclidean PH cannot see a loop that holds K562's top-2-PC share of variance
   at 150 cells.** Its negative K562 result is therefore not evidence that no loop exists.
4. **Diffusion-distance H1 is the only test that ignores pure geometry (2-10 %) and
   still finds loops** (68-100 % down to 32 % variance). It is partly sensitive to
   clusters (30 %) and partial-knockdown mixtures (15 %).
5. **Covariance-matched surrogate tests remove covariance effects but have little
   power** at realistic noise. They flag heavy tails (95-100 %), and the H0 version
   rejects on 33 % of real K562 control-vs-control draws, so it is not calibrated on
   real cells.

Real data (% of 15 draws, diffusion H1): K562 RPL21 87, POLR2B / DLD / NCBP2 60,
control vs control 13. RPE1 CPSF3 100, ZMAT2 87, positive control 100, control vs control 0.
Primary-draw BH q < 0.05 for diffusion H1: K562 RPL21, RPE1 ZMAT2, RPE1 CPSF3.
The diffusion signal follows the strongest knockdowns. From the benchmark it is
structure beyond translation, scale and covariance, but it could be a loop or a
cell subpopulation (for example escaper cells). Telling these apart is the next step.

Caveats: draws reuse cells, so they are not independent replicates. There is no guide
or batch metadata in the mirrors. 40 replicates give about ±15 % (95 % CI) on a 50 % rate.
