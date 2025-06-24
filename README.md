# Introduction
This is the github repository for my master's internship "Assessing the universality of Bolthausen-Sznitman coalescent genealogies in models of rapid selection". It contains the code I used for running SLiM simulations and generating the graphes for several models.

Each model is associated a suffix, with respect to the following table :

| Suffix | Model |
|----------|----------|
| _ | Model A, diploid |
| H | Model A, haploid |
| E | Exponential model, diploid |
| EH | Exponential model, haploid |
| M | Mutations model, diploid |
| MH | Mutations model, haploid |
| A2 | Model A with fitness based on chromosomes, diploid |
| A2MAX | ... |
| Row 9    | Data     |
| Row 10   | Data     |
| Row 11   | Data     |
| Row 12   | Data     |
| Row 13   | Data     |

# Run simulations
To run SLiM simulations, run simulate3.py with a given suffix and population size. Don't forget to change the paths.

# Plot results
Once enough simulations are run, you can analyze data by using plot3.py.
