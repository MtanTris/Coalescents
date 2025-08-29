# Introduction
This is the github repository for my master's internship "Assessing the universality of Bolthausen-Sznitman coalescent genealogies in models of rapid selection". It contains the code I used for running SLiM simulations and generating the graphes for several models.

Each model is associated a suffix, with respect to the following table :

| Suffix | Model |
|----------|----------|
| _ARGS | Model A1, diploid, with arguments that can be changed |
| X | Model A1, with hermaphrodism |
| H | Model A, haploid |
| E | Exponential model, diploid |
| EH | Exponential model, haploid |
| M | Mutations model, diploid |
| MH | Mutations model, haploid |
| A2 | Model A2, diploid |
| A2M | Model A2 except that only the 10 males with the highest fitness are chosen during reproduction |
| A2MX | Model A2 with hermaphrodism, except that only the 10 fittest individuals can reproduce |

# Run simulations
To run SLiM simulations, run `simulate3.py` with a given suffix and population size (and eventually arguments). This file also creates the tree files and extract their data into `.json` files.  

To use `simulate3.py`, three folder paths will need to be precised in the \texttt{main} function.
- `slim\_path` is the folder where the SLiM files will be stored.
- `cwd` is the folder where the tree files are stored before being analyzed. This folder is separated from the rest in case an error arises during the analysis.
- `T2\_dest` is the folder where the `.json` files will be stored after the trees are analyzed. The `.json` files are sorted by model and by population size.

For model A, some arguments can be changed before starting the simulations. They can be used to change the dominance coefficient $s$, the number $k$ of children per female (for sexual diploids) or per individual. This implementation with arguments is only available with model A as I did not want to change the other models without properly testing them, but don't hesistate do add arguments to the other models if necessary.

# Plot results
Once enough simulations are run, you can analyze data by using plot3.py. Computing the averages sadly also takes a lot of time.

# Report
For more insight about the models, don't hesitate to check my internship report (inside the repository).
