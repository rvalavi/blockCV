## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("rvalavi/blockCV")

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
# loading the package
library(blockCV)

## ---- fig.height=5, fig.width=7.2, warning=FALSE, message=FALSE----------
# loading raster library
library(raster)
# import raster data
awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))


## ---- fig.height=4.5, fig.width=7.1--------------------------------------
# import presence-absence species data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pa_data <- SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))
# see the first few rows
head(pa_data)
# plot species data on the map
plot(awt[[1]]) # plot raster data
points(pa_data[which(pa_data$Species==1), ], col="red") # add presence points
points(pa_data[which(pa_data$Species==0), ], col="blue") # add absence points
legend(x=500000, y=8250000, legend=c("Presence","Absence"), col=c(2, 4), pch=c(1,1), bty="n")


## ------------------------------------------------------------------------
# import presence-background species data
PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pb_data <- SpatialPointsDataFrame(PB[,c("x", "y")], PB, proj4string=crs(awt))
# number of presence and background records
table(pb_data$Species)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by specified range and random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species",
                   rasterLayer = awt,
                   theRange = 68000, # size of the blocks
                   k = 5,
                   selection = 'random',
                   iteration = 250, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by rows and columns
sb3 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 6,
                    k = 5,
                    selection = 'random',
                    maskBySpecies = TRUE,
                    biomod2Format = TRUE)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by rows and random assignment
sb2 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    rasterLayer = awt,
                    rows = 6,
                    k = 3,
                    selection = 'systematic',
                    biomod2Format = TRUE)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# adding points on saptialBlock plot
sb$plots + geom_point(data = as.data.frame(coordinates(pa_data)), aes(x=x, y=y), alpha=0.7)


## ---- warning=FALSE, message=FALSE---------------------------------------
# buffering with presence-absence data
bf1 <- buffering(speciesData= pa_data,
                 species = "Species", # to count the number of presences and absences
                 theRange= 68000,
                 spDataType = "PA",
                 progress = T)

## ------------------------------------------------------------------------
# buffering with presence-background data
bf2 <- buffering(speciesData= pb_data, # presence-background data
                 species= "Species",
                 theRange= 68000,
                 spDataType = "PB",
                 addBG = TRUE, # add background data to testing folds
                 progress = T)

## ---- warning=FALSE, message=FALSE---------------------------------------
# environmental clustering
eb <- envBlock(rasterLayer = awt,
               speciesData = pa_data,
               species = "Species",
               k = 5,
               standardization = "standard", # rescale variables between 0 and 1
               rasterBlock = FALSE,
               numLimit = 50)

## ---- eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7.2----
sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        border = NULL,
                        showPlots = TRUE,
                        plotVariograms = FALSE,
                        doParallel = FALSE)


## ------------------------------------------------------------------------
# class of the output result
class(sac)

## ---- eval=TRUE----------------------------------------------------------
# summary statistics of the output
summary(sac)

## ---- fig.height=4, fig.width=7------------------------------------------
plot(sac$variograms[[1]])


## ---- eval=FALSE---------------------------------------------------------
#  # explore generated folds
#  foldExplorer(sb, awt, pa_data)

## ---- eval=FALSE---------------------------------------------------------
#  # explore the block size
#  rangeExplorer(rasterLayer = awt) # the only mandatory input
#  
#  # add species data to add them on the map
#  rangeExplorer(rasterLayer = awt,
#                speciesData = pa_data,
#                species = "Species",
#                rangeTable = NULL,
#                minRange = 30000, # limit the search domain
#                maxRange = 100000)

## ---- fig.height=4, fig.width=7------------------------------------------
# loading the libraries
library(maxnet)
library(plotROC)
# extract the raster values for the species points as a dataframe
mydata <- extract(awt, pb_data, df=TRUE)
# the extract function creates an ID column that should be excluded for the modelling
# remove the extra column
mydata <- mydata[,2:ncol(mydata)]
# create a vector of 1 (for presence) and 0 (for background)
pb <- pb_data$Species

# extract the folds in buffering object created in the previous section (with presence-background data)
folds <- bf2$folds
# create an empty vector to store the AUC of each fold
AUCs <- vector()
for(k in 1:length(folds)){
  trainSet <- unlist(folds[[k]][1]) # extract the training set indices
  testSet <- unlist(folds[[k]][2]) # extract the testing set indices
  # fitting a maxent model using linear, quadratic and hinge features
  mx <- maxnet(pb[trainSet], mydata[trainSet, ], maxnet.formula(pb[trainSet], mydata[trainSet, ], classes="lqh"))
  testTable <- pb_data@data[testSet, ] # a table for testing predictions and reference data
  testTable$pred <- predict(mx, mydata[testSet, ], type="cloglog") # predict the test set
  # calculate AUC using calc_auc function in plotROC package
  auc <- calc_auc(ggplot(testTable, aes(m=pred, d=Species)) + geom_roc(n.cuts = 0))[3]
  AUCs[k] <- as.numeric(auc)
}

# print the mean and standard deviation of AUCs
print(mean(AUCs))
print(sd(AUCs))


## ---- warning=FALSE, message=FALSE, fig.height=3.7, fig.width=7----------
# loading the libraries
library(randomForest)
library(plotROC)
# extract the raster values for the species points as a dataframe
mydata <- extract(awt, pa_data, df=TRUE)
# adding species column to the dataframe
mydata$Species <- as.factor(pa_data$Species)
# remove extra column (ID)
mydata <- mydata[,2:ncol(mydata)]

# extract the folds in BufferedBlock object created in the previous section
folds <- bf1$folds
# create a data.frame to store the prediction of each fold (record)
testTable <- pa_data@data
testTable$pred <- NA
for(k in 1:length(folds)){
  trainSet <- unlist(folds[[k]][1]) # extract the training set indices
  testSet <- unlist(folds[[k]][2]) # extract the testing set indices
  rf <- randomForest(Species~., mydata[trainSet, ], ntree=250) # model fitting on training set
  testTable[testSet,"pred"] <- predict(rf, mydata[testSet, ], type="prob")[,2] # predict the test set
}
# calculate Area Under the ROC curve and plot the result using plotROC package
ggROC <- ggplot(testTable, aes(m=pred, d=Species)) + geom_roc(n.cuts=0, color='red') + 
  coord_equal() + geom_abline(intercept = 0, slope = 1) + theme_bw() +
  labs(x="False positive rate (1 - specificity)", y="True positive rate (sensitivity)")
auc <- calc_auc(ggROC)[3]
plot(ggROC + ggtitle('', subtitle=paste("AUC for testing dataset:", signif(auc, 4))))


## ---- results="hide", warning=FALSE, message=FALSE, eval=FALSE-----------
#  # loading the library
#  library(biomod2)
#  # species occurrences
#  DataSpecies <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#  # the name of studied species
#  myRespName <- "Species"
#  # the presence/absences data for our species
#  myResp <- as.numeric(DataSpecies[,myRespName])
#  # the XY coordinates of species data
#  myRespXY <- DataSpecies[,c("x","y")]
#  # change the RasterBrick to RasterStack
#  awt <- stack(awt)
#  
#  # 1. Formatting Data
#  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
#                                       expl.var = awt, # explanatory raster data
#                                       resp.xy = myRespXY,
#                                       resp.name = myRespName,
#                                       na.rm = TRUE)
#  
#  # 2. Defining the folds for DataSplitTable
#  # note that biomodTable should be used here not folds
#  DataSplitTable <- sb$biomodTable # use generated folds from spatialBlock in previous section
#  
#  # 3. Defining Models Options using default options.
#  myBiomodOption <- BIOMOD_ModelingOptions()
#  
#  # 4. Model fitting
#  myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
#                                       models = c('GLM','MARS','GBM'),
#                                       models.options = myBiomodOption,
#                                       DataSplitTable = DataSplitTable, # blocking folds
#                                       VarImport = 0,
#                                       models.eval.meth = c('ROC'),
#                                       do.full.models=FALSE,
#                                       modeling.id="test")
#  

## ---- eval=FALSE---------------------------------------------------------
#  # 5. Model evaluation
#  # get all models evaluation
#  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
#  myBiomodModelEval["ROC","Testing.data",,,]
#  

