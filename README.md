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

# Use case example using real-world data
In this example, we will be using real data from this study (Wiersma et al., 2024, [https://doi.org/10.1002/tpg2.20523](https://doi.org/10.1002/tpg2.20523)). 
In this paper, the authors used genome wide association to study anthracnose and bean common moscaic virus resistance in common beans.
They identified a locus with many disease resistance genes associated with both anthracnose and bean common moscaic virus resistance. 
In this example, let's visualize this locus across 2 resistant genotypes and 3 susceptible genotypes.

## Regions
There are 5 genotypes, and thus 5 regions. They are:

1. G19833, the reference genome and susceptible. 
2. Sacramento, susceptible
3. TARS-HT2, susceptible
4. TARS-HT1, resistant
5. Red Hawk, resistant 

```{r}
bean_region1 <- read_excel("../Data/BeanRegion1.xlsx")
bean_region2 <- read_excel("../Data/BeanRegion2.xlsx")
bean_region3 <- read_excel("../Data/BeanRegion3.xlsx")
bean_region4 <- read_excel("../Data/BeanRegion4.xlsx")
bean_region5 <- read_excel("../Data/BeanRegion5.xlsx")
```


## Attributes
```{r}
bean_attributes <- read_excel("../Data/BeanAttributes.xlsx")
head(bean_attributes)
```
The attribute table has gene_ID, functional annotation, accession (short-hand code for genotypes), and line names (genotype). 
There is no group information in this attribute table, so we have to pull that information from somewhere. 

## Load orthogroup/syntelog information 
There are multiple ways to pull orthogroup/syntelog information. 
In this example, we are parsing [GENESPACE](https://github.com/jtlovell/GENESPACE) pangenes.txt output file. 
To run this example, please unzip the `pvul_v2_pangenes.txt.gz` file in the Data/ folder. 

```{r}
bean_pangenes <- read_delim("../Data/pvul_v2_pangenes.txt", delim = "\t")

bean_pangenes_roi <- bean_pangenes %>% 
  filter(str_detect(flag, "NSOrtho") == F) %>% 
  filter(id %in% str_remove_all(bean_attributes$gene_ID, "\\.\\d+$")) %>%  # Filter genes in region of interest 
  select(id, og)
```

To parse the pangenes.txt output from GENESPACE, we filter for genes in the attribute table.
I noticed that the gene IDs in the attribute table contains isoform information e.g., ".1" at the end of the gene ID, while the `id` column in the pangenes.txt output does not.  
So to do the filtering, I trimmed the gene_ID from the attribute table using `str_remove_all(attributes$gene_ID, "\\.\\d+$")`. 
In regular expression, `\\.` matches ".", and `\\d+` matches any numbers, and `$` matches the end of line. 
What we need from the pangenes.txt output is the `og` (orthogroup) column.
This indicates if two genes from different genomes are homologous. 
I removed non-syntenic ortholog, which is reflected by the `flag` column. 
The value `array` in the `flag` column refers to tandem duplicates. 

## Add information to attribute table 
```{r}
bean_attributes_with_og <- bean_attributes %>% 
  mutate(id = str_remove_all(gene_ID, "\\.\\d+$")) %>% 
  left_join(bean_pangenes_roi, by = "id") %>% 
  mutate(tag = case_when(
    str_detect(functional_annotation, "TIR-NBS-LRR") ~ "TIR-NBS-LRR",
    T ~ "Other"
  )) %>% 
  mutate(strand = ".")

bean_attributes_with_og
```

In this study, the authors were interested in disease resistance protein in the TIR-NBS-LRR class in these regions. 
A new attribute called `tag` is added. When the functional annotation contains the words "TIR-NBS-LRR", the tag will be "TIR-NBS-LRR". 
As you can see, there are quite a number of these TIR-NBS-LRR genes in these regions. 
Lastly, since strandedness is not provided with this dataset, we will in the 'strand' column with a space holder ".". This does not affect the workflow.  

# Define range and midpoint
Next, we can align these regions to their mid point for visualization. 
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
Now we can create a data frame for visualization. 

```{r}
my_bean_regions <- rbind(
  bean_region1 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "G19833"),
  
  bean_region2 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "Sacramento"),
  
  bean_region3 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "TARS-HT2"),
  
  bean_region4 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "TARS-HT1"),
  
  bean_region5 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "Red Hawk")
  
) %>% 
  left_join(bean_attributes_with_og, by = "gene_ID") %>% 
  mutate(gene_mid = start_new + (end_new-start_new)/2) %>% 
  mutate(
    strand.text = case_when(
      strand == "+" ~ ">",
      strand == "-" ~ "<",
      T ~ " "
    )
  ) %>% 
  mutate(group = as.character(og)
  ) %>% 
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

head(my_bean_regions)
```

What happend here? 

1. We added genotype information to each region and shifted them relative to their mid points.
2. We bound these regions together as rows. 
3. Since `group` information was not provided in the attribute, we used the orthogroup (`og`) column from GENESPACE output. This dictates how genes from different genotypes should be grouped and connected by ribbons in the visualization.  
4. We calculated mid point for each gene.
5. Based on the strandedness of the gene, we added ">" or "<" symbols. Since this information is not provided, it does not apply to this example. 
6. Flip start and end position based on strandedness of genes, again not applicable to this example. 

## Other decorations to the data frame for visualization 
I noticed some of the genes are in the same orthogroup as TIR proteins, but they are not annotated as such.
I want to add the "TIR" tag to those genes too. 

```{r}
TIR_ogs <- my_bean_regions %>% 
  filter(tag == "TIR-NBS-LRR")

my_bean_regions <- my_bean_regions %>% 
  mutate(tag2 = case_when(
    og %in% TIR_ogs$og ~ "TIR-NBS-LRR",
    T ~ "Other"
  ))
```

# Plot
Now we can graph it. We need species/genotype (as a numeric variable) on y axis, and genomic coordinates relative to the region midpoint on x-axis. We also need to specify which order the genotypes will appear. 

We can specify the order of y axis using `factor(levels = c(...))`. In side `c(...)`,
the first in the list will show up at the bottom, and the last in the list will show up at the top of the ggplot. 
Let's say we want the two resistant lines (Red Hawk and TARS-HT1) at the bottom, and the three susceptible lines on top. 
For no particular region, let's say we want the order to be `"Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"`. 

```{r}
my_bean_regions %>% 
  mutate(line_name.f = factor(line_name, levels = c(
    "Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"
  ))) %>% 
  mutate(line_name.n = as.numeric(line_name.f)) %>% 
  ggplot(aes(x = gene_mid, y = line_name.n)) +
  geom_hline(aes(yintercept = line_name.n)) +
  geom_ribbon(aes(group = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3, fill = "grey80") +
  geom_ribbon(data = . %>% 
                 filter(str_detect(tag2, "TIR")),
     aes(group = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3, fill = brewer.pal(8, "Accent")[6]) +
  geom_rect(aes(xmin = start_strand, xmax = end_strand, 
                ymin = line_name.n - 0.05, ymax = line_name.n + 0.05,
                fill = tag2)
            ) +
  geom_text(aes(label = strand.text, x = gene_mid)) +
  scale_fill_manual(values = c(brewer.pal(8, "Accent")[6], "grey70"),
                    limits = c("TIR-NBS-LRR", "Other")) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c(
    "Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"
  )) +
  labs(fill = "Gene family") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", hjust = 0.5),
        legend.position = "bottom")

ggsave("../Results/bean_TIL_NBS_LRR.svg", height = 3, width = 4)
ggsave("../Results/bean_TIL_NBS_LRR.png", height = 3, width = 4)
```
![bean example 1](https://github.com/cxli233/TidyLocalSynteny/blob/main/Results/bean_TIL_NBS_LRR.png)

A few elements are visualized in this code chunk.

1. `geom_hline(aes(yintercept = line_name.n))` draws a horizontal black line for each genotype.  
2. `geom_ribbon()` draws ribbons connecting genes between genotypes. If genes are in the same group, a ribbon will be drawn connecting them. This is why the `group` column is required in the attribute table. The default is drawing grey ribbon between genes of the same group across species. 
3. Within `geom_ribbon()`, we have the option to highlight only certain groups by filtering the data: `data = . %>% filter(...)` will allow you to do that.
4. Within `geom_ribbon()`, we have the option to color the ribbons by any columns in the attribute table, such as `group`, using `aes(fill = ...)`. In this  example, we drew pink ribbon for genes in the TIR-NBS-LRR disease resistant gene family. 
5. `geom_rect()` draws rectangles for genes. 
6. Within `geom_rect()`, we have the option to color the rectangles by any columns in the attribute table, using `aes(fill = ...)`. In this  example, we drew pink boxes for genes in the TIR-NBS-LRR disease resistant gene family, and grey boxes for other genes.  
7. We have option to label the strandedness of genes, using `geom_text()`, where ">" symbols are used to label + stranded genes, and "<" symbols are used to label - stranded genes. The ">" or "<" symbol will appear at the midpoint of the gene. This is not available for this dataset, as the attribute table does not contain strand information. 
8. We have the option to specify the colors using `scale_fill_manual(...)`. In this example, I specified pink (the 6th color of "Accent" palette from the [RColorBrewer package](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html)) and grey.  
9. In `scale_y_continuous(...)`, we have to label genotypes on the y axis. To do so, first we specify y axis should break at 1, 2, 3, 4, 5 in this example, because we have 5 genotypes. Then we use `labels = c(...)` to fill in the names of species. The order inside `c(...)` should match how we specified the order of y axis, in this case it should be `"Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"`. 
10. The rest of the code turns off axis titles, axis lines, axis ticks, x axis text, and adjust y axis text color and alignment. 

# Flipping an accession, if needed 
We can see that in this example, the orientation of TARS-HT2 genotype is flipped relative to the rest of the 4 genotypes.
We can flip that. 

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

## With flipped TARS-HT2
Then we simply run it with `flip_region()` run before `shift_midpoint()` for the regions that need to be flipped. 

```{r}
my_bean_regions2 <- rbind(
  bean_region1 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "G19833"),
  
  bean_region2 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "Sacramento"),
  
  bean_region3 %>% 
  flip_region() %>% 
  shift_midpoint() %>% 
  mutate(genotype = "TARS-HT2"),
  
  bean_region4 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "TARS-HT1"),
  
  bean_region5 %>% 
  shift_midpoint() %>% 
  mutate(genotype = "Red Hawk")
  
) %>% 
  left_join(bean_attributes_with_og, by = "gene_ID") %>% 
  mutate(gene_mid = start_new + (end_new-start_new)/2) %>% 
  mutate(strand = case_when( 
    genotype == "TARS-HT2" & strand == "+" ~ "-",
    genotype == "TARS-HT2" & strand == "-" ~ "+",
    T ~ strand
  )) %>% 
  mutate(
    strand.text = case_when(
      strand == "+" ~ ">",
      strand == "-" ~ "<",
      T ~ " "
    )
  ) %>% 
  mutate(group = as.character(og)) %>% 
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

head(my_bean_regions2)
```

```{r}
my_bean_regions2 <- my_bean_regions2 %>% 
  mutate(tag2 = case_when(
    og %in% TIR_ogs$og ~ "TIR-NBS-LRR",
    T ~ "Other"
  ))
```

```{r}
my_bean_regions2 %>% 
  mutate(line_name.f = factor(line_name, levels = c(
    "Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"
  ))) %>% 
  mutate(line_name.n = as.numeric(line_name.f)) %>% 
  ggplot(aes(x = gene_mid, y = line_name.n)) +
  geom_hline(aes(yintercept = line_name.n)) +
  geom_ribbon(aes(group = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3, fill = "grey80") +
  geom_ribbon(data = . %>% 
                 filter(str_detect(tag2, "TIR")),
     aes(group = group,
                  xmin = start_strand, xmax = end_strand),
              alpha = 0.3, fill = brewer.pal(8, "Accent")[6]) +
  geom_rect(aes(xmin = start_strand, xmax = end_strand, 
                ymin = line_name.n - 0.05, ymax = line_name.n + 0.05,
                fill = tag2)
            ) +
  geom_text(aes(label = strand.text, x = gene_mid)) +
  scale_fill_manual(values = c(brewer.pal(8, "Accent")[6], "grey70"),
                    limits = c("TIR-NBS-LRR", "Other")) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5),
                     labels = c(
    "Red Hawk", "TARS-HT1", "TARS-HT2", "G19833", "Sacramento"
  )) +
  labs(fill = "Gene family") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", hjust = 0.5),
        legend.position = "bottom")

ggsave("../Results/bean_TIL_NBS_LRR_2.svg", height = 3, width = 4)
ggsave("../Results/bean_TIL_NBS_LRR_2.png", height = 3, width = 4)
```
![bean example 2](https://github.com/cxli233/TidyLocalSynteny/blob/main/Results/bean_TIL_NBS_LRR_2.png)

Done! 
