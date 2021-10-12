library("WGCNA")

#####
#Qs
#####
#Line 28
#Line 69-108

#####
#Tutorial 1: Analysis of a single empirical gene expression data set
#####

#The data are gene expression measurements from livers of female mouse of a specific F2 intercross

#Paper: Ghazalpour et al (2006), Integrating Genetics and Network Analysis to Characterize Genes Related to Mouse Weight

#data set contains 3600 measured expression profiles, filtered from the original over 20,000 profiles keeping only the most variant and most connected probes

     ##### 
     #1. Data input and cleaning
     #####
  
          #####
          #1.a Loading expression data
          #####
          
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
          
          #####
          #1.b Checking data for excessive missing values and identification of outlier microarray samples
          #####
          
          ### we need to remove genes and samples with too many missing values
          
          #search for samples with good samples genes: https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/goodSamplesGenes
          #it checks data for missing entries, 
          #entries with weights below a threshold 
          #and zero-variance genes.
          #It returns a list of samples and genes that pass criteria on maximum number of missing or low weight values. 
          #If necessary, the filtering is iterated.
          gsg = goodSamplesGenes(datExpr0, verbose = 4)
          
          gsg$allOK
          
          ### If the last statement returns TRUE, all genes have passed the cuts
          
          ###otherwise we remove those that didn't with
          if (!gsg$allOK)
          {
            # Optionally, print the gene and sample names that were removed:
            if (sum(!gsg$goodGenes)>0)
              printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
            if (sum(!gsg$goodSamples)>0)
              printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
            # Remove the offending genes and samples from the data:
            datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
          }
          
          ### Next we cluster the samples to see if there are any obvious outliers. 
          ###in contrast to clustering genes that will come later)
          
          #hclust package: Hierarchical Clustering
          sampleTree = hclust(dist(datExpr0), method = "average")
          
          ###dist= distance meaning distance between samples
          ###euclidian method often used
          
          # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
          # The user should change the dimensions if the window is too large or too small.
          sizeGrWindow(12,9)
          pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
          #cex font size
          par(cex = 0.6)
          
          #mar is margin in a unit of measuremen
          #c(down, left, up, right)
          par(mar = c(5,4,2,0)) 
          
          
          plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
               cex.axis = 1.5, cex.main = 2)
        
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

          
          #####
          #1.c Loading clinical trait data
          #####
          
          
          ### We need to match the trait samples to the same samples in the expression data.
          
          ### Read in the trait data
          traitData = read.csv("data/ClinicalTraits.csv")
          
          dim(traitData); names(traitData)
          
          # remove columns that hold information we do not need.
          allTraits = traitData[, -c(31, 16)]
          allTraits = allTraits[, c(2, 11:36) ]
          dim(allTraits)
          names(allTraits)
          
          # Form a data frame analogous to expression data that will hold the clinical traits.
          
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
          
          #We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable
          #datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits
          #relate to the sample dendrogram.
          
          # Re-cluster samples
          sampleTree2 = hclust(dist(datExpr), method = "average")
          
          # Convert traits to a color representation: white means low, red means high, grey means missing entry
          traitColors = numbers2colors(datTraits, signed = FALSE);
          
          # Plot the sample dendrogram and the colors underneath.
          plotDendroAndColors(sampleTree2, traitColors,
                              groupLabels = names(datTraits),
                              main = "Sample dendrogram and trait heatmap")
          
          #In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.
          
          #save expression and trait data for use in the next steps
          
          save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
          