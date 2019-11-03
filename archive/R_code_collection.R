newdata1 <- with(.data, data.frame(Gewinn, Verlust)
                 
                 
                 #### have to check if there are missing NA's in the model.pred vector
                 names.numeric <- as.numeric((names(model.pred)))
                 trials <- seq(from = 1, to = length(.data[,1]), by = 1) # the original length of trial number vector
                 
                 for (l in 1:length(.data[,1])) {test.na <- which(names.numeric==trials[l]) 
                                                 if (length(test.na) == 0){
                                                   dropped.cur.name <- as.character(trials[l])
                                                   model.pred[dropped.cur.name] <- NA}
                                                 
                 }
                 names.numeric <- as.numeric(names(model.pred))
                 
                 sortbyname <- function(x, y, ...) x[order(y, ...)]
                 model.pred <- sortbyname(model.pred, names.numeric)