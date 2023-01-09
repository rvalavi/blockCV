## ---- eval=FALSE, fig.height=4, fig.width=7-----------------------------------
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

## ---- echo=FALSE--------------------------------------------------------------
# The model fitting is not run to save the vignette generation time
# this AUC is based on the actual run
print(0.8664762)

## ---- eval=FALSE, fig.height=3.7, fig.width=7---------------------------------
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
#  # extract the fold indices from buffering object
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

## ----warning=FALSE, message=FALSE, eval=FALSE---------------------------------
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
#  # use generated folds from spatialBlock in previous section
#  DataSplitTable <- sb$biomodTable
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

## ---- eval=FALSE--------------------------------------------------------------
#  # 5. Model evaluation
#  # get all models evaluation
#  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
#  myBiomodModelEval["ROC","Testing.data",,,]
#  

