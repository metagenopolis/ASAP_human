The [CroCoDeEL](https://github.com/metagenopolis/CroCoDeEL) commands used to process the data :

##### Search for contamination
```
crocodeel search_conta -s merged_gut_oral_species_abundance.tsv -c contamination_events.tsv
```

##### Visualisation of plots and manual inspection of cross-contamination events
```
crocodeel plot_conta -s merged_gut_oral_species_abundance.tsv -c contamination_events.tsv -r contamination_events.pdf
```
As indicated by crocodeel authors: 
For non-related samples, CroCoDeEL may occasionally generate false positives that can be filtered out by a human-expert. Thus, they strongly recommend inspecting scatterplots of each contamination event to discard potential false positives. 
Thus, we manually inspected the scatterplots and confirmed the cross-contamination events and removed them from the analysis. 