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
          options(stringsAsFactors = FALSE) ### Q: Why do we do this? ###
          
          #load data
          femData = read.csv("data/LiverFemale3600.csv") 
          #each row corresponds to a gene and column to a sample or auxiliary information
          
          #dimensions and columns
          #(femData); names(femData)
          
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
              printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
            if (sum(!gsg$goodSamples)>0)
              printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
            # Remove the offending genes and samples from the data:
            datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
          }
          
          ### Next we cluster the samples to see if there are any obvious outliers. 
          ###in contrast to clustering genes that will come later)
          
          #hclust package: Hierarchical Clustering
          sampleTree = hclust(dist(datExpr0), method = "average")
          
          # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
          # The user should change the dimensions if the window is too large or too small.
          sizeGrWindow(12,9)
          #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
          par(cex = 0.6);
          par(mar = c(0,4,2,0))
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