The [meteor](https://github.com/metagenopolis/meteor) commands used to process the data :

##### Download the reference catalogues (gut & oral)

```         
meteor download -i hs_10_4_gut -c -o /ASAP_humain/catalogue
```

```         
meteor download -i hs_8_4_oral -c -o /ASAP_humain/catalogue
```

##### Mapping the reads against the reference catalogues

```         
meteor mapping -i ~/ASAP_human/fecal_down18 -r hs_10_4_gut -o ~/ASAP_human/down18/mapping/gut/
```

```         
meteor mapping -i ~/ASAP_human/fecal_down18 -r hs_8_4_oral -o ~/ASAP_human/down18/mapping/oral/
```

##### Compute taxonomical and functional abundances

```         
meteor profile -i ~/ASAP_human/down18/mapping/gut -o ~/ASAP_human/down18/profile/gut/ -r hs_10_4_gut -n fpkm
```

```         
meteor profile -i ~/ASAP_human/down18/mapping/oral -o ~/ASAP_human/down18/profile/oral/ -r hs_8_4_oral -n fpkm
```

##### Merging output from different samples in a single table

```         
meteor merge -i ~/ASAP_human/down18/profile/gut/ -r hs_10_4_gut -o ~/ASAP_human/down18/merging/gut/
```

```         
meteor merge -i ~/ASAP_human/down18/profile/oral/ -r hs_8_4_oral -o ~/ASAP_human/down18/merging/oral/
```
