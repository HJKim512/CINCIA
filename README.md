# CINCIA

CINCIA, which represents 'Computational Inference of Cell-Cell Interactions from scRNA-seq Data', is an R package for inferring potential hetero-typic cell-cell interactions using scRNA-seq data. In this package, putative doublets are predicted by three different doublet callers, subjected to a meta-caller inferring a posterior probability of cell-cell interaction.  



### Contents

 CINCIA contains four R scripts where each set of functions implements different parts.  
 
* **Processing.R**  
  In this set of functions, the input expression matrix is processed through Seurat v4 package (Hao, Y., et al., 2021) resulting in the reference set which contains putative singlets with minimized number of putative multiplets.  
  
  
* **Detect_Doublets.R**
  Those functions detect putative doublets through three widely-used doublet detection tools including Hybrid (Bais and Kostka, 2019), scDblFinder (Germain, et al., 2021) and Scrublet (Wolock, et al., 2019).  


* **Define_Doublets.R**  
  In this script, putative doublets are classified by the Gaussian Naive Bayes based doublet meta-caller, which uses the doublet score of each droplet calculated previously (Detect_Doublets.R), and each droplet in the input data is assigned a specific class label, either singlet or doublet.


* **Infer_interactions.R**  
  Finally, doublets are assigned the two most likely cell types, and putative heterotypic cell-cell interactions are inferred from the deconvoluted putative doublets.  
  
<br>

### Installation

Install from github
```
devtools::install_github("HJKim512/CINCIA")
library(CINCIA)
```


<br>

### Analysis Example  


#### 1. Processing


```
# UMI/gene counts plots are saved at OutPath. Put the input data into the function in the form of 
# either an expression matrix (param 'expMtx') or a Seurat object(param 'rawObj')
# If you pass the gene expression matrix, the raw Seurat object ('rawObj') would be saved in the OutPath as 'outPath/1.Processing/tmp/rawObj.Rds'.
Lower_plot_Preprocess(expMtx = matrix, OutPath = "out/dir/path/", project.id = "proj_name")

# First, low quality cells with lower gene/UMI count under the thresholds (user input parameters, 
# 'thres.gene' and 'thres.umi') are filtered out. Upper gene counts are filtered by user-input threshold 
# (thres.UpGene) or by EM algorithm for Gaussian Mixture models (thres.prob).
obj.flt <- Upper_Preprocess( obj, thres.gene=200, thres.umi=500, OutPath="out/dir/path/", thres.prob=0.7 )

# Percent of mitochondrial RNA is calculated for each droplet, and user can manually determine the 
# threshold based on the plots saved in OutPath. 
obj.flt <- mtRNA_plot_Preprocess(obj.flt, prefix = 'MT-', OutPath="out/dir/path/")
obj.flt.mtcut <- subset(obj.flt, subset = percent.mt < 20)  # put the proper threshold
saveRDS(obj.flt.mtcut, "out/dir/path/1.Processing/tmp/rawObj.mtcut.Rds")

# This function will normalize the data, find highly variable genes, scale the data and perform PCA.
# Check the JackStrawPlot generated in the OutPath, and determine the optimal number of PCs (nPC) by comparing their p-values.
obj.norm <- NormToPCA_Process(Obj.flt=obj.flt.mtcut, OutPath="out/dir/path", num.dim=50)

# Generate a vector of marker genes for generating feature plots.
markers <-  c("Dclk1", "Chga", "Muc2", "Alpi", "Slc26a3", "Mki67", "Lgr5", "Lyz1")
# UMAP, clustering and DE analysis are sequentially performed on the input data.
# You can use any value you want for the random seed number ('seed') and clustering resolution ('res').
obj.norm <- UMAPtoDEG_Process(Obj.norm=obj.norm, nPC=nPC, genes=markers, seed=77777, res=1.0, "out/dir/path")

# Manually determine the cell type of each droplet and save as a column named 'celltype' in the meta.data.
# of the object. Make sure you set the cell type without character '-'. 
#    example code :  
obj.norm$celltype <- mapvalues(obj.norm.flt$seurat_clusters, 
                               from = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                               to = c("Paneth", "stem", "Transit amplifying", "Progenitor early", 
                                      "Progenitor late", "Enterocyte", "Goblet", "Enteroendocrine", "Tuft")) 
saveRDS(obj.norm, 'out/dir/path/1.Processing/tmp/rawObj.Annotated.Rds')

# Running DoubletFinder() to remove putative doublets in the reference set
# Put nPC, estimated multiplet rate of the input data (est.rates), random seed number (rseed). 
obj.sgl <- Run_DoubletFinder(obj.norm, npc = 25, nLib = 1, nRun = 1, rseed = 4044, est.rates = 0.1, OutPath = "out/dir/path" )
saveRDS(obj.sgl, 'out/dir/path/1.Processing/tmp/Reference.Set.Rds')

# Dimension plot of the Reference set
DimPlot(obj.sgl, reduction = "umap", group.by = 'celltype', pt.size = .05, label = T, label.size = 1.5)+
  theme(axis.title = element_text(size=7), axis.text = element_text(size=7), plot.title = element_text(size=7),
      legend.text = element_text(size=5), legend.key.height = unit(0.5,'cm'), legend.key.width = unit(0.5,'cm'))
ggsave('out/dir/path/1.Processing/stats/2.Dimplot.Celltype.Singlets.pdf', units = 'cm', width = 13, height = 10)

```



#### 2. Doublet detection through meta-caller

```

# Here, doublet simulation using the reference set is performed. For nDbl, use estimated doublet numbers according to  
Simulate_Hetero_Doublets( obj.sgl, nDbl=500, col.ref="celltype", OutPath="out/path/for/detection")

# In this step, doublet score of each droplet of the input data (type='raw') and simulated data (type='simul) is calculated 
# by three doublet detection tools. In 'PATH' directory, 'tmp/Reference.Set.Rds', input expression matrix ("raw/~.ExpMtx.nonZero.txt") 
# and label table ("raw/~.LabelTable.txt"), which are previously generated in the processing step, should be exist. For the parameter 
# 'tool', one among 'Hyb' (Hybrid), 'Scr' (Scrublet) and 'scDF' (scDblFinder) should be given, where 'PyPath' should be provided 
# when running Scrublet (tool='Scr'). It can take quite long for simulated datasets (type='simul').
#     example :
DoubletScoring( PATH='out/path/for/detection/', type='simul', tool='Scr', est.rate=0.1, nSet=100, nRun=10, PyPath="path/to/python")
DoubletScoring( PATH='out/path/for/detection/', type='simul', tool='Hyb', est.rate=0.1, nSet=100, nRun=10)
DoubletScoring( PATH='out/path/for/detection/', type='raw', tool='scDF', est.rate=0.1, nSet=100, nRun=10)


# Putative doublets included in the data are finally defined.
dblobj <- DefineDoublets(PATH="out/path/for/detection/", PATH_out = "out/dir/path/", nRun=10, nSet_s = 100, nSet_a = 1)


```

#### 3. Cell-type assignment and inference

```

# For each putative doublet, top two cell types that are most likely to compose itself are assigned.
# In this function, the Seurat object of the input raw data (rawObj) which is saved in the output directory, selected number of PCs, and threshold for cell type prediction score filtering (default 0.98) are used.

Deconv.obj <- Assign_Celltypes( refObj = obj.sgl, rawObj = rawObj, NBres = dblobj, col.label = "celltype",
                                OutPath = "out/dir/path/", nPC = PCnum, thres = 0.98)
                                
# Finally, putative heterotypic cell-cell interactions are inferred. Cell type fraction of the input data is calculated first to estimate the expected number of doublets. For each cell type pair observed in the previous step, the number of corresponding kind of doublets is compared to the expected number. The result object ('Interaction.obj') contains significant heterotypic cell-cell interactions inferred from the whole processes of CINCIA.

Interaction.obj <- FindInteraction(DeconvObj = Deconv.obj, refObj = ReferenceSet, col.label = "celltype", 
                                   lower.frac = 0.1, OutPath = "out/dir/path/")
                                       
```


<br>

#### Doublet detection tools used

[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
McGinnis, C.S. et al., Cell Systems 2019
<br>[Hybrid(scds)](https://github.com/kostkalab/scds)
Bais, A.S. and Kostka, D., Bioinformatics 2019
<br>[Scrublet](https://github.com/swolock/scrublet)
Wolock, S.L., Lopez, R. and Klein, A.M. Cell Syst 2019
<br>[scDblFinder](https://github.com/plger/scDblFinder) 
Germain, P., et al., F1000Research 2021
