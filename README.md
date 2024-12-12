# TidyLocalSynteny
A Tidyverse-based workflow for local synteny visualization 

![Example output](https://github.com/cxli233/TidyLocalSynteny/blob/main/Results/example1_2.png)

* Author: Chenxin Li, Ph.D., Assistant Research Scientist at Department of Crop & Soil Sciences and Center for Applied Genetic Technologies, University of Georgia. 
* Contact: [Chenxin.Li@uga.edu](Chenxin.Li@uga.edu) | [@ChenxinLi2](https://twitter.com/ChenxinLi2) | [@chenxinli2.bsky.social](https://bsky.app/profile/chenxinli2.bsky.social)

# Introduction
This is a workflow to generate visualization for local synteny using Tidyverse in R. 
This tutorial comes with a hypothetical example dataset and a real-world example dataset (Wiersma et al., 2024, [https://doi.org/10.1002/tpg2.20523](https://doi.org/10.1002/tpg2.20523)). 

## Dependencies
```{r}
library(tidyverse)
library(readxl)
library(RColorBrewer)
```
Technically `tidyverse` is the only dependency. `readxl` is only required when the input data is excel format. 
`RColorBrewer` is only required if you are using color palettes in the [RColorBrewer package](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html). 

# Data
Required input:

1. Region tables: Genomic coordinates of genes with 4 columns: Chr, start, end, gene_ID.
2. Each region should have its own region table. Each region table can be a different species, genotype, or a different region of the same genome.   
3. Attribute table with the following required columns: gene_ID, group, strand. gene_IDs must be unique across genomes/species/accessions, and _match exactly_ with the gene_IDs in the region tables. `group` here refers to homologous groups. One can obtain homology information using BLASTp or software such as [GENESPACE](https://github.com/jtlovell/GENESPACE). This is the table where you include any other information useful for visualization, such as gene symbol. 

__Limitations__:  For the group column in attribute table, each species/accession/genome should have 1 gene per group, i.e., 1:1:1 relationships across species, genotypes, or regions. When 1:n (n > 1) cases occurs, one of the paralogs must be declared as _the_ homolog (same group as homolog in the other species), whereas the others declared as tandem arrays (a different group assignment from _the_ homolog).  

Here is an example using hypothetical data. 
You can find them in the `Data/` folder of this repository. 

```{r}
region1 <- read_excel("../Data/ExampleRegion1.xlsx")
region2 <- read_excel("../Data/ExampleRegion2.xlsx")
region3 <- read_excel("../Data/ExampleRegion3.xlsx")

attributes <- read_excel("../Data/ExampleAttributes1.xlsx")
```

# Define range and midpoint
The first step is to align the regions of interest to their mid point using the following function: 
```{r}
shift_midpoint <- function(region){
  
  colnames(region) <- c("Chr", "start", "end", "gene_ID")
  
  region %>%
    mutate(min = min(start)) %>% 
    mutate(max = max(end)) %>% 
    mutate(mid = min+(max - min)/2) %>% 
    mutate(start_new = start - mid) %>% 
    mutate(end_new = end - mid)
}

```

Then we create a data frame for visualization. 
```{r}
my_regions <- rbind(
  region1 %>% 
  shift_midpoint() %>% 
  mutate(species = "species1"),
  
  region2 %>% 
  shift_midpoint() %>% 
  mutate(species = "species2"),
  
  region3 %>% 
  shift_midpoint() %>% 
  mutate(species = "species3")
) %>% 
  left_join(attributes, by = "gene_ID") %>% 
  mutate(gene_mid = start_new + (end_new-start_new)/2) %>% 
  mutate(
    strand.text = case_when(
      strand == "+" ~ ">",
      strand == "-" ~ "<",
      T ~ " "
    )
  ) %>% 
  mutate(group = case_when(
    is.na(group) ~ "Other",
    T ~ group
  )) %>% 
  mutate(start_strand = case_when(
    strand == "+" ~ start_new,
    strand == "-" ~ end_new,
    T ~ start_new
  )) %>% 
  mutate(end_strand = case_when(
    strand == "+" ~ end_new,
    strand == "-" ~ start_new,
    T ~ end_new
  )) 

my_regions
```
What happened here: 

1. We added species information to each region and calculated new coordinates relative to the mid points of the respective regions.
2. We bound these regions together as rows. 
3. We calculated mid point for each gene.
5. Based on the strandedness of the gene, we added ">" sign for + strand and "<" sign for - strand.
6. Flip start and end positions based on strandedness of genes. If a gene is on the - strand, start becomes end, and end becomes start. 

# Visualization 
Now we can graph it. We need species/genotype (as a numeric variable) on y axis, and genomic coordinates relative to the region midpoint on x-axis. We also need to specify which order the genotypes will appear. 

We can specify the order of y axis using `factor(levels = c(...))`. In side `c(...)`,
the first in the list will show up at the bottom, and the last in the list will show up at the top of the ggplot. 

```{r}
my_regions %>% 
  mutate(species.f = factor(species, levels = c(
    "species1", "species2", "species3"
  ))) %>% 
  mutate(species.n = as.numeric(species.f)) %>% 
  ggplot(aes(x = gene_mid, y = species.n)) +
  geom_hline(aes(yintercept = species.n)) +
  geom_ribbon(data = . %>%
                filter(group != "Other"),
              aes(group = group, fill = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3) +
  geom_rect(aes(xmin = start_strand, xmax = end_strand, 
                ymin = species.n - 0.05, ymax = species.n + 0.05,
                fill = group)
            ) +
  geom_text(aes(label = strand.text, x = gene_mid)) +
  scale_fill_manual(values = c(brewer.pal(8, "Set2")[1:4], "grey80")) +
  scale_y_continuous(breaks = c(1, 2, 3),
                     labels = c("Species1", "Species2", "Species3")) +
  labs(fill = "Homologs") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"))

ggsave("../Results/example1.svg", height = 2, width = 4)
ggsave("../Results/example1.png", height = 2, width = 4)
```
![Example output1](https://github.com/cxli233/TidyLocalSynteny/blob/main/Results/example1.png)

A few elements are visualized in this code chunk.

1. `geom_hline(aes(yintercept = species.n))` draws a horizontal black line for each species. 
2. `geom_ribbon()` draws ribbons connecting genes between species. If genes are in the same group, a ribbon will be drawn connecting them. This is why the `group` column is required in the attribute table. 
3. Within `geom_ribbon()`, we have the option to highlight only certain groups by filtering the data: `data = . %>% filter(...)` will allow you to do that. In this example, we filtered for groups that are not "Others". 
4. Within `geom_ribbon()`, we have the option to color the ribbons by any columns in the attribute table, such as `group`, using `aes(fill = ...)`.
5. `geom_rect()` draws rectangles for genes. 
6. Within `geom_rect()`, we have the option to color the rectangles by any columns in the attribute table, using `aes(fill = ...)`. 
7. We have option to label the strandedness of genes, using `geom_text()`, where ">" symbols are used to label + stranded genes, and "<" symbols are used to label - stranded genes. The ">" or "<" symbol will appear at the midpoint of the gene. 
8. We have the option to specify the colors using `scale_fill_manual(...)`. 
9. In `scale_y_continuous(...)`, we have to label species on the y axis. To do so, first we specify y axis should break at 1, 2, and 3 in this example. In you have 5 species, the y axis should break at 1, 2, 3, 4, 5. Then we use `labels = c(...)` to fill in the names of species. The order inside `c(...)` should match how we specified the order of y axis. 
10. The rest of the code turns off axis titles, axis lines, axis ticks, x axis text, and adjust y axis text color and alignment. 

# Flipping a species (if need). 
We can see that in this example, the orientation of species 3 is flipped relative to the rest of the 2 genotypes. We can flip that if we want. 
```{r}
flip_region <- function(region){
  region %>% 
    mutate(start = max(end)-start) %>% 
    mutate(end = max(end)-end) %>% 
    mutate(start1 = end) %>% 
    mutate(end1 = start) %>% 
    mutate(start = start1) %>% 
    mutate(end = end1) %>% 
    select(Chr, start, end, gene_ID)
}
```

To flip a region, we simply take the end point of a region and subtract each coordinate of that region to the end point.
Then we flip start to end, and end to start.

## With flipped species3
Then we simply run it with `flip_region()` run before `shift_midpoint()` for the regions that need to be flipped.
Since we have flipped the region for species 3, we also need to flip the strand (change + into - and vice versa) for species  3. 

```{r}
my_regions2 <- rbind(
  region1 %>% 
  shift_midpoint() %>% 
  mutate(species = "species1"),
  
  region2 %>% 
  shift_midpoint() %>% 
  mutate(species = "species2"),
  
  region3 %>% 
  flip_region() %>% 
  shift_midpoint() %>% 
  mutate(species = "species3")
) %>% 
  left_join(attributes, by = "gene_ID") %>% 
  mutate(strand = case_when( 
    species == "species3" & strand == "+" ~ "-",
    species == "species3" & strand == "-" ~ "+",
    T ~ strand
  )) %>% 
  mutate(gene_mid = start_new + (end_new-start_new)/2) %>% 
  mutate(
    strand.text = case_when(
      strand == "+" ~ ">",
      strand == "-" ~ "<",
      T ~ " "
    )
  ) %>% 
  mutate(group = case_when(
    is.na(group) ~ "Other",
    T ~ group
  )) %>% 
  mutate(start_strand = case_when(
    strand == "+" ~ start_new,
    strand == "-" ~ end_new,
    T ~ start_new
  )) %>% 
  mutate(end_strand = case_when(
    strand == "+" ~ end_new,
    strand == "-" ~ start_new,
    T ~ end_new
  )) 

my_regions2
```
```{r}
my_regions2 %>% 
  mutate(species.f = factor(species, levels = c(
    "species1", "species2", "species3"
  ))) %>% 
  mutate(species.n = as.numeric(species.f)) %>% 
  ggplot(aes(x = gene_mid, y = species.n)) +
  geom_hline(aes(yintercept = species.n)) +
  geom_ribbon(data = . %>%
                filter(group != "Other"),
              aes(group = group, fill = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3) +
  geom_rect(aes(xmin = start_strand, xmax = end_strand, 
                ymin = species.n - 0.05, ymax = species.n + 0.05,
                fill = group)
            ) +
  geom_text(aes(label = strand.text, x = gene_mid)) +
  scale_fill_manual(values = c(brewer.pal(8, "Set2")[1:4], "grey80")) +
  scale_y_continuous(breaks = c(1, 2, 3),
                     labels = c("Species1", "Species2", "Species3")) +
  labs(fill = "Homologs") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"))

ggsave("../Results/example1_2.svg", height = 2, width = 4)
ggsave("../Results/example1_2.png", height = 2, width = 4)
```
![Example output2](https://github.com/cxli233/TidyLocalSynteny/blob/main/Results/example1_2.png)

Done!
