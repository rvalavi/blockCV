## ---- eval=FALSE---------------------------------------------------------
#  remotes::install_github("rvalavi/blockCV", dependencies = TRUE)

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
# loading the package
library(blockCV)


## ---- fig.height=5, fig.width=7.2, warning=FALSE, message=FALSE----------
# loading raster library
library(raster)
library(sf)

# import raster data
awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))


## ---- fig.height=4.5, fig.width=7.1--------------------------------------
# import presence-absence species data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pa_data <- st_as_sf(PA, coords = c("x", "y"), crs = crs(awt))
# see the first few rows
pa_data

# plot species data on the map
plot(awt[[1]]) # plot raster data
plot(pa_data[which(pa_data$Species==1), ], pch = 16, col="red", add=TRUE) # add presence points
plot(pa_data[which(pa_data$Species==0), ], pch = 16, col="blue", add=TRUE) # add absence points
legend(x=500000, y=8250000, legend=c("Presence","Absence"), col=c(2, 4), pch=c(16,16), bty="n")


## ------------------------------------------------------------------------
# import presence-background species data
PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pb_data <- st_as_sf(PB, coords = c("x", "y"), crs = crs(awt))
# number of presence and background records
table(pb_data$Species)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by specified range with random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species",
                   rasterLayer = awt,
                   theRange = 70000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by rows and columns with checkerboard assignment
sb2 <- spatialBlock(speciesData = pb_data, # presence-background data
                    species = "Species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 6,
                    k = 5,
                    selection = "systematic",
                    biomod2Format = TRUE)


## ----eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by rows with systematic assignment
sb3 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    rasterLayer = awt,
                    rows = 6,
                    selection = "checkerboard",
                    biomod2Format = TRUE)


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7-------------
# adding points on saptialBlock plot
library(ggplot2)

sb$plots + geom_sf(data = pa_data, alpha = 0.5)


## ----eval=FALSE, warning=FALSE, message=FALSE----------------------------
#  # buffering with presence-absence data
#  bf1 <- buffering(speciesData = pa_data,
#                   theRange = 70000,
#                   species = "Species", # to count the number of presences and absences/backgrounds
#                   spDataType = "PA", # presence-absence  data type
#                   progress = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  # buffering with presence-background data
#  bf2 <- buffering(speciesData = pb_data, # presence-background data
#                   theRange = 70000,
#                   species = "Species",
#                   spDataType = "PB", # presence-background data type
#                   addBG = TRUE, # add background data to testing folds
#                   progress = TRUE)
#  

## ----eval=FALSE, warning=FALSE, message=FALSE----------------------------
#  # environmental clustering
#  eb <- envBlock(rasterLayer = awt,
#                 speciesData = pa_data,
#                 species = "Species",
#                 k = 5,
#                 standardization = "standard", # rescale variables between 0 and 1
#                 rasterBlock = FALSE,
#                 numLimit = 50)

## ---- eval=TRUE, warning=FALSE, message=FALSE, fig.height=5, fig.width=7.2----
sac <- spatialAutoRange(rasterLayer = awt,
                        sampleNumber = 5000,
                        doParallel = TRUE,
                        showPlots = TRUE)


## ------------------------------------------------------------------------
# class of the output result
class(sac)

## ---- eval=TRUE----------------------------------------------------------
# summary statistics of the output
summary(sac)

## ---- eval=FALSE, fig.height=4, fig.width=7------------------------------
#  library(automap)
#  
#  plot(sac$variograms[[1]])
#  

## ---- eval=FALSE---------------------------------------------------------
#  # explore generated folds
#  foldExplorer(blocks = sb,
#               rasterLayer = awt,
#               speciesData = pa_data)

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

## ---- eval=FALSE, fig.height=4, fig.width=7------------------------------
#  # loading the libraries
#  library(maxnet)
#  library(precrec)
#  # library(ggplot2)
#  
#  # extract the raster values for the species points as a dataframe
#  mydata <- raster::extract(awt, pb_data)
#  mydata <- as.data.frame(mydata)
#  # create a vector of 1 (for presence) and 0 (for background samples)
#  pb <- pb_data$Species
#  
#  # extract the folds in spatialBlock object created
#  # in the previous section (with presence-background data)
#  # the foldID only works for spatialBlock and envBlock folds
#  folds <- sb2$foldID
#  
#  # create an empty vector to store the AUC of each fold
#  AUCs <- vector(mode = "numeric")
#  for(k in seq_len(5)){
#    # extracting the training and testing indices
#    # this way only works with foldID
#    trainSet <- which(folds != k) # training set indices
#    testSet <- which(folds == k) # testing set indices
#    # fitting a maxent model using linear, quadratic and hinge features
#    mx <- maxnet(p = pb[trainSet],
#                 data = mydata[trainSet, ],
#                 maxnet.formula(p = pb[trainSet],
#                                data = mydata[trainSet, ],
#                                classes = "default"))
#    testTable <- pb_data[testSet, ] # a table for testing predictions and reference data
#    testTable$pred <- predict(mx, mydata[testSet, ], type = "cloglog") # predict the test set
#    # calculate area under the ROC curve
#    precrec_obj <- evalmod(scores = testTable$pred, labels = testTable$Species)
#    AUCs[k] <- auc(precrec_obj)[1,4] # extract AUC-ROC
#  }
#  
#  # print the mean of AUCs
#  print(mean(AUCs))
#  

## ---- echo=FALSE---------------------------------------------------------
# The model fitting is not run to save the vignette generation time
# this AUC is based on the actual run
print(0.8664762)

## ---- eval=FALSE, fig.height=3.7, fig.width=7----------------------------
#  # loading the libraries
#  library(randomForest)
#  library(precrec)
#  
#  # extract the raster values for the species points as a dataframe
#  mydata <- raster::extract(awt, pa_data, df = TRUE)
#  # adding species column to the dataframe
#  mydata$Species <- as.factor(pa_data$Species)
#  # remove extra column (ID)
#  mydata <- mydata[,-1]
#  
#  # extract the foldIDs in SpatialBlock object
#  # created in the previous section
#  # the folds (list) works for all three blocking strategies
#  folds <- bf1$folds
#  
#  # create a data.frame to store the prediction of each fold (record)
#  testTable <- pa_data
#  testTable$pred <- NA
#  
#  for(k in seq_len(length(folds))){
#    # extracting the training and testing indices
#    # this way works with folds list (but not foldID)
#    trainSet <- unlist(folds[[k]][1]) # training set indices
#    testSet <- unlist(folds[[k]][2]) # testing set indices
#    rf <- randomForest(Species~., mydata[trainSet, ], ntree = 250) # model fitting on training set
#    testTable$pred[testSet] <- predict(rf, mydata[testSet, ], type = "prob")[,2] # predict the test set
#  }
#  
#  # calculate Area Under the ROC and PR curves and plot the result
#  precrec_obj <- evalmod(scores = testTable$pred, labels = testTable$Species)
#  
#  autoplot(precrec_obj)
#  

## ----warning=FALSE, message=FALSE, eval=FALSE----------------------------
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

