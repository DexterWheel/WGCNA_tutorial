library("WGCNA")

#Qs
#line 212: what is approximate scale-free topology
#line 712: -abs?

#####
#Tutorial 1: Analysis of a single empirical gene expression data set
#####

#The data are gene expression measurements from livers of female mouse of a specific F2 intercross

#Paper: Ghazalpour et al (2006), Integrating Genetics and Network Analysis to Characterize Genes Related to Mouse Weight

#data set contains 3600 measured expression profiles, filtered from the original over 20,000 profiles keeping only the most variant and most connected probes

     #==== 
     #1. Data input and cleaning
     #====
  
          #====
          #1.a Loading expression data
          #====
          
          ### "expression data is contained in LiverFemale3600.csv"
          
          #Load package
          library(WGCNA)
          options(stringsAsFactors = FALSE) ### factors only store one name one time, if there are repeats you will lose them unless they are strings ###
          
          #load data
          femData = read.csv("data/LiverFemale3600.csv") 
          #each row corresponds to a gene and column to a sample or auxiliary information
          
          #dimensions and columns
          #dim(femData); names(femData)
          
          ### "the data files contain extra information about the surveyed probes we do not need"
          
          #remove auxiliary data and transpose the expression data
          datExpr0 = as.data.frame( t(femData[, -c(1:8)]) )
          
          #Replacing col names with that of the data from the substance BXH column
          names(datExpr0) = femData$substanceBXH
          
          #Replacing rownames with the column names starting with the 9th column
          rownames(datExpr0) = names(femData)[-c(1:8)]
          rownames(datExpr0)
          
          #====
          #1.b Checking data for excessive missing values and identification of outlier microarray samples
          #====
          
          ### we need to remove genes and samples with too many missing values
          
          #search for samples with good samples genes: https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/goodSamplesGenes
          #it checks data for missing entries, 
          #entries with weights below a threshold 
          #and zero-variance genes.
          #It returns a list of samples and genes that pass criteria on maximum number of missing or low weight values. 
          #If necessary, the filtering is iterated.
          gsg = goodSamplesGenes(datExpr0, verbose = 3)
          gsg$allOK
          
          ### If the last statement returns TRUE, all genes have passed the cuts
          
          ###otherwise we remove those that didn't with
          #if (!gsg$allOK)
          #{
            # Optionally, print the gene and sample names that were removed:
           # if (sum(!gsg$goodGenes)>0)
              #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
            #if (sum(!gsg$goodSamples)>0)
             # printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
            # Remove the offending genes and samples from the data:
            #datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
          #}
          
          
          ### Next we cluster the samples to see if there are any obvious outliers. 
          ###in contrast to clustering genes that will come later)
          
          #hclust package: Hierarchical Clustering
          sampleTree = hclust(dist(datExpr0), method = "average")
          
          ###dist= distance meaning distance between samples
          ###euclidian method often used
          
          # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
          # The user should change the dimensions if the window is too large or too small.
          sizeGrWindow(12,9)
          #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
          #cex font size
          par(cex = 0.6);
          #mar is margin in a unit of measuremen
          #c(down, left, up, right)
          par(mar = c(0,4,2,0))
          
          plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
               cex.axis = 1.5, cex.main = 2);
          
          ### It appears there is one outlier sample F2_221
          ### One can remove it by hand, or use an automatic approach
          ###Choose a height cut that will remove the offending sample, say 15 (the red line in the plot), and use a branch cut at that height.
          
          # Plot a line to show the cut
          abline(h = 15, col = "red")
          # Determine cluster under the line
          clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
          table(clust)
          # clust 1 contains the samples we want to keep.
          keepSamples = (clust==1)
          datExpr = datExpr0[keepSamples, ]
          nGenes = ncol(datExpr)
          nSamples = nrow(datExpr)
          
          ###The variable datExpr now contains the expression data ready for network analysis.

          
          #====
          #1.c Loading clinical trait data
          #====
          
          
          ### We need to match the trait samples to the same samples in the expression data.
          
          ### Read in the trait data
          traitData = read.csv("data/ClinicalTraits.csv")
          
          dim(traitData); names(traitData)
          
          # remove columns that hold information we do not need.
          allTraits = traitData[, -c(31, 16)]
          allTraits = allTraits[, c(2, 11:36) ]
          dim(allTraits)
          names(allTraits)
          
          ### Form a data frame analogous to expression data that will hold the clinical traits.
          
          #taking the rows of the finished expression dataset
          femaleSamples = rownames(datExpr)
          
          #mice is the row that contains the sample names
          #this is using the match function to take only the samples that are present in both datasets
          traitRows = match(femaleSamples, allTraits$Mice)
          
          #subsetting alltraits to just contain the matched samples
          datTraits = allTraits[traitRows, -1]
          
          #renaming the rows 
          rownames(datTraits) = allTraits[traitRows, 1]
          
          #collects temporary files
          collectGarbage()
          
          ###We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable
          ###datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits
          ###relate to the sample dendrogram.
          
          # Re-cluster samples
          sampleTree2 = hclust(dist(datExpr), method = "average")
          
          # Convert traits to a color representation: white means low, red means high, grey means missing entry
          traitColors = numbers2colors(datTraits, signed = FALSE);
          
          # Plot the sample dendrogram and the colors underneath.
          plotDendroAndColors(sampleTree2, traitColors,
                              groupLabels = names(datTraits),
                              main = "Sample dendrogram and trait heatmap")
          
          ###In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.
          
          ###save expression and trait data for use in the next steps
          
          save(datExpr, datTraits, file = "data-processed/FemaleLiver-01-dataInput.RData")
          
          
          
     #==== 
     #2. Network construction and module detection
     #====
          
          # Load the data saved in the first part
          lnames = load(file = "data-processed/FemaleLiver-01-dataInput.RData");
          #The variable lnames contains the names of loaded variables.
          lnames
          
          
          #This step is the bedrock of all network analyses using the WGCNA methodology/ 
          #We present three different ways of constructing a network and identifying modules:
          
          #a. Using a convenient 1-step network construction and module detection function, suitable for users wishing to arrive at the result with minimum effort;
          
          #b. Step-by-step network construction and module detection for users who would like to experiment with customized/alternate methods;
          
          #c. An automatic block-wise network construction and module detection method for users who wish to analyze data
          #sets too large to be analyzed all in one.
          
          ###In this tutorial section, we illustrate the 1-step, automatic network construction and module detection.###
          
          
          
          #====
          #2.a Automatic, one-step network construction and module detection:
          #====
          
          #====
          #2.a.1: Choosing the soft-thresholding power: analysis of network topology
          #====
               
               ###Weighted gene network construction requires a soft thresholding power β to which co-expression similarity is raised to calculate adjacency. 
                
               ###It is suggested to choose the soft thresholding power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details  zhang and horvath 2005; 
               
               
               ###the function pickSoftThreshold performs analysis of network topology and aids the user in choosing a proper soft-thresholding power. 
               
               ###The user chooses a set of candidate powers (the function provides suitable default values), 
               ###and the function returns a set of network indices that should be inspected, for example as follows:
               
               # Choose a set of soft-thresholding powers
               powers = c(c(1:10), seq(from = 12, to=20, by=2))
               
               powers
               
               # Call the network topology analysis function
               sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
               
               # Plot the results:
               sizeGrWindow(9, 5)
               par(mfrow = c(1, 2));
               cex1 = 0.9;
               
               # Scale-free topology fit index as a function of the soft-thresholding power
               plot(sft$fitIndices[,1], 
                    -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                    xlab = "Soft Threshold (power)", 
                    ylab = "Scale Free Topology Model Fit,signed R^2",type = "n",
                    main = paste("Scale independence"));
               
               text(sft$fitIndices[,1],
                    -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                    labels=powers,cex=cex1,col="red");
               
               # this line corresponds to using an R^2 cut-off of h
               abline(h=0.90,col="red")
               
               # Mean connectivity as a function of the soft-thresholding power
               plot(sft$fitIndices[,1], 
                    sft$fitIndices[,5],
                    xlab = "Soft Threshold (power)", 
                    ylab = "Mean Connectivity", type = "n",
                    main = paste("Mean connectivity"))
               
               text(sft$fitIndices[,1], 
                    sft$fitIndices[,5], 
                    labels=powers, 
                    cex=cex1,col="red")
               
               #The result is shown in Fig. 1. We choose the power 6, which is the lowest power for which the scale-free topology fit
               #index curve flattens out upon reaching a high value (in this case, roughly 0.90).
               
               
          #====
          #2.a.2: One-step network construction and module detection
          #====
               
               #Constructing the gene network and identifying modules is now a simple function call:
                 
               net = blockwiseModules(datExpr, power = 6,
                                      TOMType = "unsigned", minModuleSize = 30,
                                      reassignThreshold = 0, mergeCutHeight = 0.25,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      saveTOMs = TRUE,
                                      saveTOMFileBase = "data-processed/femaleMouseTOM",
                                      verbose = 3)
               
               
               
               
               
##               #We have chosen the soft thresholding power 6, a relatively large minimum module size of 30, and a medium sensitivity (deepSplit=2) to cluster splitting. 
               
##               ###The parameter mergeCutHeight is the threshold for merging of modules. 
               
##               ###We have also instructed the function to return numeric, rather than color, labels for modules, and to save the Topological Overlap Matrix. 
               
##               ###The output of the function may seem somewhat cryptic, but it is easy to use. 
##               ###For example, net$colors contains the module assignment, and net$MEs contains the module eigengenes of the modules.
               
#Important#    ###A word of caution for the readers who would like to adapt this code for their own data. The function blockwiseModules has many parameters, and in this example most of them are left at their default value.
               
#Important#    ###We have attempted to provide reasonable default values, but they may not be appropriate for the particular data set the reader wishes to analyze.
               
##               ###We encourage the user to read the help file provided within the package in the R environment and experiment with tweaking the network construction and module detection parameters. The potential reward is, of course, better (biologically more relevant) results of the analysis.
               
##               ###A second word of caution concerning block size. 
##               ###In particular, the parameter maxBlockSize tells the function how large the largest block can be that the reader’s computer can handle.
               
##               ###The default value is 5000 which is appropriate for most modern desktops. Note that if this code were to be used to analyze a data set with more than 5000 probes, the function blockwiseModules will split the data set into several blocks. This will break some of the plotting code below, that is executing the code will lead to errors. Readers wishing to analyze larger data sets need
               #to do one of the following:
              
##               #• If the reader has access to a large workstation with more than 4 GB of memory, the parameter maxBlockSize
##               #can be increased. A 16GB workstation should handle up to 20000 probes; a 32GB workstation should handle
##               #perhaps 30000. A 4GB standard desktop or a laptop may handle up to 8000-10000 probes, depending on
##               #operating system and how much memory is in use by other running programs.
               
               ###• If a computer with large-enough memory is not available, the reader should follow Section 2.c, Dealing withlarge datasets, and adapt the code presented there for their needs. In general it is preferable to analyze a data set in one block if possible, although in Section 2.c we present a comparison of block-wise and single-block analysis that indicates that the results are very similar.
               
               
               
               
               
               #We now return to the network analysis. 
               #To see how many modules were identified and what the module sizes are, one can use table(net$colors).
               
               table(net$colors)
               
               #This indicates that there are 18 modules, 
               #labeled 1 through 18 in order of descending size, 
               #with sizes ranging from 609 to 34 genes. 
               #The label 0 is reserved for genes outside of all modules.
               
               #The hierarchical clustering dendrogram (tree) used for the module identification is returned in net$dendrograms[[1]]; #$.
               
               #The dendrogram can be displayed together with the color assignment using the following code:
               
               
               # open a graphics window
               sizeGrWindow(12, 9)
               # Convert labels to colors for plotting
               mergedColors = labels2colors(net$colors)
               # Plot the dendrogram and the module colors underneath
               plotDendroAndColors(net$dendrograms[[1]], 
                                   mergedColors[net$blockGenes[[1]]],
                                   "Module colors",
                                   dendroLabels = FALSE, hang = 0.03,
                                   addGuide = TRUE, guideHang = 0.05)
               
               
               
#Important     #We note that if the user would like to change some of the tree cut, module membership, and module merging criteria, the package provides the function recutBlockwiseTrees that can apply modified criteria without having to recompute the network and the clustering dendrogram. This may save a substantial amount of time.
               
               
               #We now save the module assignment and module eigengene information necessary for subsequent analysis.
               
               moduleLabels = net$colors
               moduleColors = labels2colors(net$colors)
               MEs = net$MEs;
               geneTree = net$dendrograms[[1]];
               save(MEs, moduleLabels, 
                    moduleColors, 
                    geneTree,
                    file = "data-processed/FemaleLiver-02-networkConstruction-auto.RData")
               
          #====
          #2.b  Step-by-step network construction and module detection:     
          #====
               
          #====
          #2.b.1   Choosing the soft-thresholding power: analysis of network topology     
          #====     
               
               ###Code for this section is the same as 2.a.1
               
          #====
          #2.b.2   Co-expression similarity and adjacency     
          #====      
               
               ###We now calculate the adjacencies, using the soft thresholding power 6:
               
               softPower = 6;
               adjacency = adjacency(datExpr, power = softPower)
               
          #====
          #2.b.3   Topological Overlap Matrix (TOM)
          #==== 
               
               ###To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
               
               # Turn adjacency into topological overlap
               
               TOM = TOMsimilarity(adjacency)
               dissTOM = 1-TOM
               
          #====
          #2.b.4   Clustering using TOM
          #====
               
               
               ###We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. 
               
               ###Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.
               
               # Call the hierarchical clustering function
               geneTree = hclust(as.dist(dissTOM), method = "average")
               
               # Plot the resulting clustering tree (dendrogram)
               sizeGrWindow(12,9)
               plot(geneTree, xlab="", sub="", 
                    main = "Gene clustering on TOM-based dissimilarity",
                    labels = FALSE, hang = 0.04)
               
               
               ###The clustering dendrogram plotted by the last command is shown in Figure 2. 
               
               ###In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene 
               ###Branches of the dendrogram group together densely interconnected, highly co-expressed genes.
               
               ###Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”).
               
               ###There are several methods for branch cutting; our standard method is:
               
               ###Dynamic Tree Cut from the package dynamicTreeCut. 
               
               ###The next snippet of code illustrates its use.
               
               
               # We like large modules, so we set the minimum module size relatively high:
               minModuleSize = 30;
               # Module identification using dynamic tree cut:
               dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                           deepSplit = 2, pamRespectsDendro = FALSE,
                                           minClusterSize = minModuleSize);
               table(dynamicMods)
               
               
               #The function returned 22 modules labeled 1–22 largest to smallest. Label 0 is reserved for unassigned genes.
               
               #The above command lists the sizes of the modules. We now plot the module assignment under the gene dendrogram:
                 
               
               # Convert numeric lables into colors
               dynamicColors = labels2colors(dynamicMods)
               table(dynamicColors)
               # Plot the dendrogram and colors underneath
               sizeGrWindow(8,6)
               plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                                   dendroLabels = FALSE, hang = 0.03,
                                   addGuide = TRUE, guideHang = 0.05,
                                   main = "Gene dendrogram and module colors")
               
               
               
               
               
               
               
               
               
               
               
          #====
          #2.b.5  Merging of modules whose expression profiles are very similar
          #====
               
               #The Dynamic Tree Cut may identify modules whose expression profiles are very similar.
               
               #It may be prudent to merge such modules since their genes are highly co-expressed.
               
               #To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:
                 
               # Calculate eigengenes
               MEList = moduleEigengenes(datExpr, colors = dynamicColors)
               MEs = MEList$eigengenes
               
               # Calculate dissimilarity of module eigengenes
               MEDiss = 1-cor(MEs);
               
               # Cluster module eigengenes
               METree = hclust(as.dist(MEDiss), method = "average");
               
               # Plot the result
               sizeGrWindow(7, 6)
               plot(METree, main = "Clustering of module eigengenes",
                    xlab = "", sub = "")
               
               #We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge (see Fig. 4):
               
               MEDissThres = 0.25
               # Plot the cut line into the dendrogram
               abline(h=MEDissThres, col = "red")
               # Call an automatic merging function
               merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
               # The merged module colors
               mergedColors = merge$colors;
               # Eigengenes of the new merged modules:
               mergedMEs = merge$newMEs;
               
               #To see what the merging did to our module colors, we plot the gene dendrogram again, 
               #with the original and merged module colors underneath (Figure 5).
               
               sizeGrWindow(12, 9)
               #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
               plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                                   c("Dynamic Tree Cut", "Merged dynamic"),
                                   dendroLabels = FALSE, hang = 0.03,
                                   addGuide = TRUE, guideHang = 0.05)
               #dev.off()
               
               #In the subsequent analysis, we will use the merged module colors in mergedColors. We save the relevant variables for use in subsequent parts of the tutorial:
               
               # Rename to moduleColors
               moduleColors = mergedColors
               
               # Construct numerical labels corresponding to the colors
               colorOrder = c("grey", standardColors(50));
               moduleLabels = match(moduleColors, colorOrder)-1;
               MEs = mergedMEs;
               # Save module colors and labels for use in subsequent parts
               save(MEs, moduleLabels, moduleColors, geneTree, file = "data-processed/FemaleLiver-02-networkConstruction-stepByStep.RData")
               
               
     #==== 
     #3.  Relating modules to external information and identifying important genes
     #====
               
           # Load the WGCNA package
           library(WGCNA)
           # The following setting is important, do not omit.
           options(stringsAsFactors = FALSE);
           # Load the expression and trait data saved in the first part
           lnames = load(file = "data-processed/FemaleLiver-01-dataInput.RData");
           #The variable lnames contains the names of loaded variables.
           lnames
           # Load network data saved in the second part.
           lnames = load(file = "data-processed/FemaleLiver-02-networkConstruction-auto.RData");
           lnames
               
          #====
          #3.a  Quantifying module–trait associations
          #====
                
               #In this analysis we would like to identify modules that are significantly associated with the measured clinical traits.
                
               #Since we already have a summary profile (eigengene) for each module, 
               #we simply correlate eigengenes with external traits and look for the most significant associations:
               
               
               # Define numbers of genes and samples
               nGenes = ncol(datExpr)
               nSamples = nrow(datExpr)
               
               # Recalculate MEs with color labels
               MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
               MEs = orderMEs(MEs0)
               moduleTraitCor = cor(MEs, datTraits, use = "p")
               moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
               
               #Since we have a moderately large number of modules and traits, 
               #a suitable graphical representation will help in reading the table. 
               #We color code each association by the correlation value:
               
               sizeGrWindow(10,6)
               # Will display correlations and their p-values
               textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                                  signif(moduleTraitPvalue, 1), ")", sep = "")
               
               dim(textMatrix) = dim(moduleTraitCor)
               
               par(mar = c(6, 8.5, 3, 3))
               
               # Display the correlation values within a heatmap plot
               labeledHeatmap(Matrix = moduleTraitCor,
                              xLabels = names(datTraits),
                              yLabels = names(MEs),
                              ySymbols = names(MEs),
                              colorLabels = FALSE,
                              colors = greenWhiteRed(50),
                              textMatrix = textMatrix,
                              setStdMargins = FALSE,
                              cex.text = 0.5,
                              zlim = c(-1,1),
                              main = paste("Module-trait relationships"))
                 
               
               #The analysis identifies the several significant module–trait associations. 
               #We will concentrate on weight as the trait of interest.
               
               
               
          #====
          #3.b  Gene relationship to trait and important modules: Gene Significance and Module Membership
          #====    
                 
               ###We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. 
                 
               ###For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile.
                 
               ###This allows us to quantify the similarity of all genes on the array to every module.
                 
                 
               # Define variable weight containing the weight column of datTrait
               weight = as.data.frame(datTraits$weight_g);
               #rename column
               names(weight) = "weight"
                 
               # names (colors) of the modules
               modNames = substring(names(MEs), 3)
              
               geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
               MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
                 
               names(geneModuleMembership) = paste("MM", modNames, sep="")
               names(MMPvalue) = paste("p.MM", modNames, sep="")
                 
               geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
               GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
                 
               names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
               names(GSPvalue) = paste("p.GS.", names(weight), sep="")
              
          #====
          #3.c   Intramodular analysis: identifying genes with high GS and MM
          #====      
              
               ###Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.
               
               ###As an example, we look at the brown module that has the highest association with weight.
               
               ###We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
               
               
               module = "brown"
               column = match(module, modNames);
               
               moduleGenes = moduleColors==module;
               
               sizeGrWindow(7, 7);
               
               par(mfrow = c(1,1));
               
               verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                                  abs(geneTraitSignificance[moduleGenes, 1]),
                                  xlab = paste("Module Membership in", module, "module"),
                                  ylab = "Gene significance for body weight",
                                  main = paste("Module membership vs. gene significance\n"),
                                  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
               
               #The plot is shown in Fig. 2. Clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait.
               
               #The reader is encouraged to try this code with other significance trait/module correlation (for example, the magenta, midnightblue, and red modules with weight).
               
               
          #====
          #3.d    Summary output of network analysis results
          #====      
               
               ### We have found modules with high association with our trait of interest, 
               #and have identified their central players by the Module Membership measure.
               
               ###We now merge this statistical information with gene annotation and write out a file that
               #summarizes the most important results and can be inspected in standard spreadsheet software 
               #such as MS Excel or Open Office Calc.
               
               ###Our expression data are only annotated by probe ID names: the command
               names(datExpr)
               ###will return all probe IDs included in the analysis. Similarly,
               
               names(datExpr)[moduleColors=="brown"]
               ###will return probe IDs belonging to the brown module. 
               
               ###To facilitate interpretation of the results, 
               #we use a probe annotation file provided by the manufacturer of the expression arrays
               #to connect probe IDs to gene names and universally recognized identification numbers (Entrez codes).
               
               
               annot = read.csv(file = "data/GeneAnnotation.csv")
               dim(annot)
               names(annot)
               probes = names(datExpr)
               probes2annot = match(probes, annot$substanceBXH)
               
               # The following is the number or probes without annotation:
               sum(is.na(probes2annot))
               # Should return 0.
               
               ###We now create a data frame holding the following information for all probes: 
               #probe ID, gene symbol, Locus Link ID (Entrez code), module color, gene significance for weight, and module membership and p-values in all modules. 
               
               ###The modules will be ordered by their significance for weight, 
               #with the most significant ones to the left.
               
               # Create the starting data frame
               geneInfo0 = data.frame(substanceBXH = probes,
                                      geneSymbol = annot$gene_symbol[probes2annot],
                                      LocusLinkID = annot$LocusLinkID[probes2annot],
                                      moduleColor = moduleColors,
                                      geneTraitSignificance,
                                      GSPvalue)
               
               # Order modules by their significance for weight
               modOrder = order(-abs(cor(MEs, weight, use = "p")))
               
               ###-abs?#
               ##
               
               # Add module membership information in the chosen order
               for (mod in 1:ncol(geneModuleMembership)){
                 oldNames = names(geneInfo0)
                 geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                                        MMPvalue[, modOrder[mod]])
                 
                 names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                      paste("p.MM.", modNames[modOrder[mod]], sep=""))
                 }
               # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
               geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
               geneInfo = geneInfo0[geneOrder, ]
               
               #This data frame can be written into a text-format spreadsheet, for example by
               write.csv(geneInfo, file = "C:/Users/Dexter weighell/Documents/Biology/Yr 4/Project/R/tutorial/WGCNA_tutorial/data-processed/geneInfo.csv")
               
               #The reader is encouraged to open and view the file in a spreadsheet software, or inspect it directly within R using the command 
               fix(geneInfo)
               
     #==== 
     #4.  Interfacing network analysis with other data such as functional annotation and gene ontology
     #====      
               library(WGCNA)
               # The following setting is important, do not omit.
               options(stringsAsFactors = FALSE);
               # Load the expression and trait data saved in the first part
               lnames = load(file = "data-processed/FemaleLiver-01-dataInput.RData");
               #The variable lnames contains the names of loaded variables.
               lnames
               # Load network data saved in the second part.
               lnames = load(file = "data-processed/FemaleLiver-02-networkConstruction-auto.RData");
               lnames
               
               
          ###Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight.
               
          ###To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, whether they are significantly enriched in certain functional categories etc
               
          #==== 
          #4.a   Output gene lists for use with online software and services
          #====     
               
               ###One option is to simply export a list of gene identifiers 
               #that can be used as input for several popular gene ontology
               #and functional enrichment analysis suites such as David or AmiGO. 
               
               ###For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:
               getwd()
               # Read in the probe annotation
               annot = read.csv(file = "data/GeneAnnotation.csv")
               
               # Match probes in the data set to the probe IDs in the annotation file
               probes = names(datExpr)
               probes2annot = match(probes, annot$substanceBXH)
               
               # Get the corresponding Locuis Link IDs
               allLLIDs = annot$LocusLinkID[probes2annot]
               
               # $ Choose interesting modules
               intModules = c("brown", "red", "salmon")
               
               for (module in intModules){
                 # Select module probes
                 modGenes = (moduleColors==module)
                 # Get their entrez ID codes
                 modLLIDs = allLLIDs[modGenes]
                 # Write them into a file
                 fileName = paste("C:/Users/Dexter weighell/Documents/Biology/Yr 4/Project/R/tutorial/WGCNA_tutorial/LocusLinkIDs-", module, ".txt", sep="")
                 
                 write.table(as.data.frame(modLLIDs), file = fileName,
                             row.names = FALSE, col.names = FALSE)
               }
               
               # As background in the enrichment analysis, we will use all probes in the analysis.
               fileName = paste("LocusLinkIDs-all.txt", sep="");
               write.table(as.data.frame(allLLIDs), file = fileName,
                           row.names = FALSE, col.names = FALSE)
               
          #==== 
          #4.b   Enrichment analysis directly within R
          #====
               
               ###The WGCNA package now contains a function to perform GO enrichment analysis using a simple, single step.
               
               ###To run the function, Biconductor packages GO.db, AnnotationDBI, and the appropriate organism-specific annotation package(s) need to be installed before running this code.
               
               ###The organism-specific packages have names of the form org.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc.
               
               ###The only exception is yeast, for which no org.Xx.eg.db package is available; instead, the package carries the name org.Sc.sgd.db.
               
               ###Please visit the Bioconductor main page at http://www.bioconductor.org to download and install the required packages.
               
               ###In our case we are studying gene expressions from mice, so this code needs the package org.Mm.eg.db.
               
               ###Calling the GO enrichment analysis function GOenrichmentAnalysis is very simple. 
               
               ###The function takes a vector of module labels, and the Entrez (a.k.a. Locus Link) codes for the genes whose labels are given.
               
               
               GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10)
               
               #The function runs for awhile and returns a long list, the most interesting component of which is:
               
               tab = GOenr$bestPTerms[[4]]$enrichment
               
               #This is an enrichment table containing the 10 best terms for each module present in moduleColors.
               
               #Names of the columns within the table can be accessed by:
               names(tab)
               
               
               
               
               
               
               
               
               
               
               
               
               