Several downstream processes, like characterizing heterogeniety in an exploratory analysis, we need `genes` which are differentially expressed between any two cells. So, the process of selecting such `genes` which may tell about different existing "cell groups" within a population is called Feature selection and the group of `genes` which help us doing this are called `HVGs` (Highly variable genes). We will proceed with a 10X PBMC data set.

```r
raw.path <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.tar.gz")
out.path <- file.path(tempdir(), "pbmc4k")
untar(raw.path, exdir=out.path)
fname <- file.path(out.path, "raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

# gene annotation #
rownames(sce.pbmc) <- uniquifyFeatureNames(
    rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)
    
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")
    
 ```
 
 **explanation**
 
When we look at the `rawData(sce.pbmc)` we will find that in "ID" column all the IDs are unique for all the transcript, but it may so happens, and it happens quite often, that for any two unique IDs there are same names. So, to "uniquify" the row names  of a `singleCellExperiment` object (here sce.pbmc) `uniquifyfeaturenames(ID,name)` will append "ID" to the "name" to make a unique combination for each of `rownames()`. This function will attempt to use "name" if it is unique. If not, it will append the ID to any non-unique value of names. Missing names will be replaced entirely by ID.
The output is guaranteed to be unique, assuming that ID is also unique. This combination can be directly used as the `rownames()` of a `SingleCellExperiment` object. 


```r
# cell detection #
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

# Quality control #
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

# Normalization # 

set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)
```

**explanation**

With the droplet based techniques it may happen that some of the droplets remains empty. `emptyDrops()` function identifies droplets that contain cells. `e.out` is a dataMatrix where rows are cells and 5 columns are some values associated with each cell. One of the column is FDR. 




