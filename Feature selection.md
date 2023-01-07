Several downstream processes, like characterizing heterogeniety in an exploratory analysis, we need `genes` which are differentially expressed between any two cells. So, the process of selecting such `genes` which may tell about different existing "cell groups" within a population is called Feature selection and the group of `genes` which help us doing this are called `HVGs` (Highly variable genes). We identify `HVGs` to focus on the genes that are driving heterogeneity across the population of cells. This requires estimation of the variance in expression for each gene, followed by decomposition of the variance into biological and technical components. `HVGs` are then identified as those genes with the largest biological components We will proceed with a 10X PBMC data set.

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

With the droplet based techniques it may happen that some of the droplets remains empty. `emptyDrops()` function identifies droplets that contain cells. `e.out` is a dataMatrix where rows are cells and 5 columns are some values associated with each cell. One of the column is FDR. `which(location=="MT")` will check "location" for `MTs`, and will give out all the "location" which are stamped as `MT`.  

We need to take another data set "416b"

```r
> sce.416b <- LunSpikeInData(which = "416b")
> class(sce.416b$block)
[1] "integer"
# class(sce.416b$block) is just and integer , wnd we need to convert it in "factor" to really divide all cells in to  "block"#
sce.416b$block <- factor(sce.416b$block)

# Gene annotation #
ens.mm.v97 <- AnnotationHub()[["AH73905"]] # getting annotation object from AnnotationHub(), here a mouse database #

rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)

rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
    keytype="GENEID", column="SYMBOL")
    
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
    keytype="GENEID", column="SEQNAME")

rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
    rowData(sce.416b)$SYMBOL)

# Quality control #
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

# Normalization #
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)


 
```

### Quantifying per-gene variation 

The simplest approach to quantifying per-gene variation is to compute the variance of the log-normalized expression values (i.e., “log-counts” ) for each gene across all cells (A. T. L. Lun, McCarthy, and Marioni 2016). The advantage of this approach is that the feature selection is based on the same log-values that are used for later downstream steps. In particular, genes with the largest variances in log-values will contribute most to the Euclidean distances between cells during procedures like clustering and dimensionality reduction. By using log-values here, we ensure that our quantitative definition of heterogeneity is consistent throughout the entire analysis. Calculation of the per-gene variance is simple but feature selection requires modelling of the mean-variance relationship. The log-transformation is not a variance stabilizing transformation in most cases, which means that the total variance of a gene is driven more by its abundance than its underlying biological heterogeneity. To account for this effect, we use the `modelGeneVar()` function to fit a trend to the variance with respect to abundance across all genes. 

```r
dec.pbmc <- modelGeneVar(sce.pbmc)
# Visualizing the fit #
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
    ylab="Variance of log-expression")
    
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2) 
```

![image5](https://user-images.githubusercontent.com/85447250/210852749-9748987f-8123-4c57-92f9-134a8c0895e8.png)
 
 Fig: Variance in the PBMC data set as a function of the mean. Each point represents a gene while the blue line represents the trend fitted to all genes.
 
 ### explanation 
 
 `modelGeneVar()` calculates variation of genes (features or rows of `sce`) across all cells. This will give a `data.frame` where rows corresponds to the genes and there will be certain numbers of columns which indicates total variance, technical component of variance, biological component of variance, p.value, FDR, and mean.
 
 `metadata()` will give a `list` object which includes _mean_, _var_ , _trend_ and, _std.dev_. Among these _trend_ is a function which will define a curve that fits to the mean-variance plot. 
 
 ```r
 > fit.pbmc$trend
function (x) 
{
    output <- FUN(x) * med
    names(output) <- names(x)
    output
}
```
### end explanation
 
 
At any given abundance, we assume that the variation in expression for most genes is driven by uninteresting technical processes. Under this assumption, the fitted value of the trend at any given gene’s abundance represents an estimate of its uninteresting variation, which we call the technical component. We then define the biological component for each gene as the difference between its total variance and the technical component. This biological component represents the “interesting” variation for each gene and can be used as the metric for HVG selection. It may be noted here that if we have access to fitting fitting a mean-dependent trend to the variance of the spike-in transcript, that would constitute a better strategy to figure out the technical component of the variation.  

```r
# Ordering by most interesting genes for inspection
dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),]
```
It is important to note that, the interpretation of the fitted trend as the technical component assumes that the expression profiles of most genes are dominated by random technical noise. 

### Quantifying technical noises

The assumption in previous section may be problematic in rare scenarios where many genes at a particular abundance are affected by a biological process. For example, strong upregulation of cell type-specific genes may result in an enrichment of HVGs at high abundances. This would inflate the fitted trend in that abundance interval and compromise the detection of the relevant genes. We can avoid this problem by fitting a mean-dependent trend to the variance of the spike-in transcript, if they are available (Ideally, the technical component would be estimated by fitting a mean-variance trend to the spike-in transcriptsWe may recall that the same set of spike-ins was added in the same quantity to each cell. This means that the spike-in transcripts should exhibit no biological variability, i.e., any variance in their counts should be technical in origin). The premise here is that spike-ins should not be affected by biological variation, so the fitted value of the spike-in trend should represent a better estimate of the technical component for each gene.

```r
dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
dec.spike.416b[order(dec.spike.416b$bio, decreasing=TRUE),]

DataFrame with 46604 rows and 6 columns
              mean     total      tech       bio      p.value          FDR
         <numeric> <numeric> <numeric> <numeric>    <numeric>    <numeric>
Lyz2       6.61097   13.8497   1.57131   12.2784 1.48993e-186 1.54156e-183
Ccl9       6.67846   13.1869   1.50035   11.6866 2.21855e-185 2.19979e-182
Top2a      5.81024   14.1787   2.54776   11.6310  3.80015e-65  1.13040e-62
Cd200r3    4.83180   15.5613   4.22984   11.3314  9.46221e-24  6.08574e-22
Ccnb2      5.97776   13.1393   2.30177   10.8375  3.68706e-69  1.20193e-66
...            ...       ...       ...       ...          ...          ...
Rpl5-ps2   3.60625  0.612623   6.32853  -5.71590     0.999616     0.999726
Gm11942    3.38768  0.798570   6.51473  -5.71616     0.999459     0.999726
Gm12816    2.91276  0.838670   6.57364  -5.73497     0.999422     0.999726
Gm13623    2.72844  0.708071   6.45448  -5.74641     0.999544     0.999726
Rps12l1    3.15420  0.746615   6.59332  -5.84670     0.999522     0.999726

```

![image6](https://user-images.githubusercontent.com/85447250/211096635-1c1e5366-a220-4b45-bc58-f559a8f6beff.png)










