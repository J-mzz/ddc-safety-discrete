# Data-Driven Safe Control of Discrete-Time Non-Linear Systems

MATLAB code for the paper

> J. Zheng, J. Miller, and M. Sznaier, "Data-driven safe control of discrete-time non-linear systems," *IEEE Control Systems Letters*, vol. 8, pp. 1553–1558, 2024.
> [doi:10.1109/LCSYS.2024.3412943](https://doi.org/10.1109/LCSYS.2024.3412943) · [accepted version (NSF PAR)](https://par.nsf.gov/biblio/10534780)

The method starts from noisy data of an unknown discrete-time polynomial system and solves a sum-of-squares program (Theorem 1). The program finds a Lyapunov function that certifies safe stabilization of every system consistent with the data. A safe input is then computed online at each step (Theorem 2, or Lemmas 3 and 4 for non-convex safe sets).

## Requirements

- MATLAB (tested with R2024a) with the Symbolic Math Toolbox and the Optimization Toolbox
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licenses are available); tested with MOSEK 10.0

## Usage

Run any script from MATLAB. Each one adds `utils/` to the path, generates the data, solves Theorem 1, simulates the closed loop, and plots the results.

| Script | Paper | Description |
|---|---|---|
| `example1_safe.m` | Example 1, Fig. 1 (blue) | linear system, convex safe set, 40 steps from `x0 = [2; 1]` |
| `example1_unsafe.m` | Example 1, Fig. 1 (orange) | same system, controller designed without the safety constraint |
| `example2_safe.m` | Example 2, Fig. 2 | 3D polynomial system, half-space safe set `x1 <= 1`, 10 steps |
| `example2_unsafe.m` | Example 2, Fig. 2 | same system, controller designed without the safety constraint |
| `example3_safe.m` | Example 3, Fig. 3 | polynomial system, non-convex safe set (outside a disk), input `u` in `[-2, 2]` |
| `example3_unsafe.m` | Example 3, Fig. 3 | same system, controller designed without the safety constraint |

The noise bound, decrease rate (`c3`, called `cc` in Example 2), initial state, and number of steps are set at the top of each script. The paper's figures overlay the trajectories of the `_safe` and `_unsafe` scripts of the same example. The scripts plot the trajectories and, in Examples 1 and 2, the Lyapunov function along them; they do not plot the inputs shown in the panels (b).

## Repository layout

```
example*_safe.m, example*_unsafe.m   one script per example and controller
utils/
  getSystem.m, setCons.m     true systems (used to generate data and to simulate) and problem sizes
  getDataRND.m               noisy samples of the true system
  getConsEx1/2/3.m           data matrices of the consistency set for each example
  nontrivial_constraints.m   removes redundant constraints of a polytope
  third_party/               external code, see utils/third_party/README.md
```

## Citation

If you use this code, please cite:

```bibtex
@article{zheng2024data,
  title={Data-driven safe control of discrete-time non-linear systems},
  author={Zheng, Jian and Miller, Jared and Sznaier, Mario},
  journal={IEEE Control Systems Letters},
  volume={8},
  pages={1553--1558},
  year={2024},
  publisher={IEEE},
  doi={10.1109/LCSYS.2024.3412943}
}
```

## License

The code in this repository is released under the [MIT License](LICENSE), except the files in `utils/third_party/`, which keep their original licenses.
