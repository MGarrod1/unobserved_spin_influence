# unobserved_spin_influence

This repository contains code used to perform simulations and generate plots for the paper: "Influencing dynamics on social networks without knowledge of network microstructure"

# Requirements

This repository makes use of code within the following repos:

- [ising block level influence](https://github.com/MGarrod1/ising_block_level_influence) - used for deriving and evaluating the controls.

- [spatial spin monte carlo](https://github.com/MGarrod1/spatial_spin_monte_carlo) - used for the Monte Carlo simulations.

The code was validated using Python 3.7.4 on an Anaconda distribution (conda version : 4.8.3). A full list of the packages installed can be found in the *environment.yml* file. This can be used to construct an identical anaconda environment to that used to perform the simulations. See the Anaconda [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more details.

## Data used

The simulations associated with Figures 3a,3b, S1, S2 & S3 rely on data from the social network Pokec. This can be downloaded from: [https://snap.stanford.edu/data/soc-Pokec.html](https://snap.stanford.edu/data/soc-Pokec.html). When running the simulations this data was placed inside the directory: unobserved\_spin\_influence/Figure\_3\_Pokec/Data/raw\_data.

The .csv files obtained from the simulations are included in the Data/ directory associated with each figure. Figure\_2\_threeblock/Data/MC\_phase\_data.csv is omitted due to its size.

## Summary of scripts and data associated with figures

The table below summarises the python scripts or Jupyter notebooks used to generate each of the Figures in the paper.

| Figure    | Input data                                          | Simulations in                          | Output data                                                                               | Plots made in                           | Time taken for simulations          |
|-----------|-----------------------------------------------------|-----------------------------------------|-------------------------------------------------------------------------------------------|-----------------------------------------|-------------------------------------|
| 1 a & b   | N/A                                                 | fig2-1-AB.py                            | N/A                                                                                       | fig2-1-AB.py                            | < 2 mins                            |
| 1 c,d &e  | N/A                                                 | fig2-1-CDE-simulations.py               | two_block_markup_data_spins1-0_bf_0-5.csv                                                 | fig2-1-CDE-plots.py                     | 2 H 15 mins                         |
| 2 a       | N/A                                                 | fig2-2-A.py                             | N/A                                                                                       | fig2-2-A.py                             | Seconds                             |
| 2 b       | N/A                                                 | fig2-2-B-simulations.py                 | block_level_phase_data.csv, full_MF_phase_data.csv, MC_phase_data.csv                     | three_block_phase_diagram.ipynb         |  5H in total                        |
| 2 c,d & e | N/A                                                 | fig2-2-CDE-simulations.py               | three_block_sus_data.csv                                                                  | fig2-2-CDE-plots.ipynb                  | < 1 min                             |
| 3 a       | soc-pokec-profiles.txt soc-pokec-relationships.txt  | make_bratislava_graph_and_blocks.ipynb  | bratislava_profiles.csv Bratislava_graph.graphmlBratislava_coupling.graphmlblock_info.csv | make_bratislava_graph_and_blocks.ipynb  | Minutes                             |
| 3 b       | N/A                                                 | fig3-B-simulations.ipynb                | Pokec_control_eval_data_new.csv                                                           | fig3-B-plots.ipynb                      | 8 H                                 |
| S1        | N/A                                                 | N/A                                     | N/A                                                                                       | figS1-plots.ipynb                       | Seconds                             |
| S2        | N/A                                                 | figS2-simulations.ipynb                 | Pokec_phase_diagram_data.csv                                                              | figS2-plots.ipynb                       | 2H 40 mins running this on 07/10/20 |
| S3        | N/A                                                 | figS3-simulations.ipynb                 | snapshot_as_sampfrac_data_grad_1-0.csv                                                    | figS3-plots.ipynb                       |                                     |
