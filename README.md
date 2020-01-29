# V-D and D-J junction regions generation model 
The development of a simulation model for generation of V-D and D-J regions in immunoglobulins’ nucleotide sequences

### Goals: 
- Selection of variable parameters according to real data
- Comparison with the previously established model
- Creation of a console user interface with an option to visualize the acquired results

### Methods in VD model:
- Binomial p-nucleotide length distribution
- The length of n-nucleotides depends on the nucleotide weights and the number of nucleotides sufficient for fusion
- Exonuclease activity is presented as a Markov process with two states: exonuclease removes, does not remove nucleotides

### System requirement: 
- Python 3.x (guaranteed for version 3.7)
- Installed Python packages: `biopython`, `scipy`, `numpy`, `mathplotlib`

### Steps to run:
- Clone repository 
```
git clone https://github.com/natalia0000/immunoglobulin.git
```
- For using VD model in your project needed import `main_fun` from `VD.py` and run with four parameters (reference VD documentation)

##### For sample run:
- Put `model_data.csv` with your data set in same folder
- Change parameter `index` in `optimization.py` and run

##### For visualization:
- Change parameter `index` and put parameters that returned from `optimization.py` in `plots.py` and run

### Result example:

![alt text](https://github.com/natalia0000/immunoglobulin/blob/master/sample_plot.jpg)

Parameters for that plot:
```
p = 0.30839815
p_compl = 0.10480559
exo_start_delete = 0.23485349
exo_stop_delete = 0.20396072
```

Blue is distribution from data set

Red is distribution from VD model

### Reference

Kenneth Murphy, Casey Weaver  "Janeway’s Immunobiology, Ninth Edition." 2016. Garland Science: New York, New York.
