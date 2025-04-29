# doubleMutualist()
# A function to visualize double mutualistic interactions in a very simple and
# comprehensive plot

# Adjusted from plotweb2() in the R-package "bipartite"

doubleMutualist <- 
  
  function (

    web, web2, # two interaction matrices, one species group (double mutualist needs to be shared between matrices)
    method="cca", method2="cca", # default plotting method for web and web2 (see plotweb())
    empty=FALSE, empty2=TRUE, 
    
    labsize=1, lab.space=1, 
    
    ybig=1, y_width=0.1, 
    arrow="no", arrow2="no", 
    
    high.abun=NULL, low.abun=NULL, 
    high.abun2=NULL, low.abun2=NULL,
    
    lablength=NULL, lablength2=NULL, 
    sequence=NULL, sequence.pred2=NULL,
     
    # Define default colors. Can be customized to highlight e.g. different species
    # Colorsettings for interactions between species (in lower and upper web)
    col.interaction1="grey80", col.interaction2="grey80",
    col.borderinteraction1="black", col.borderinteraction2="black",
    
    # Colorsettings for the doublemutualists
    col.doubleMutualist="grey10",  col.doubleMutualist2="grey10", 
    col.border.doubleMutualist="grey10", col.border.doubleMutualist2="grey10",
    
    # Colorsettings for species the doublemutualists interact with (lower and upper web)
    col.intPartner="grey90", col.intPartner2="grey85",
    col.border.intPartner="grey10", col.border.intPartner2="grey10",
    
    # Colorsettings for the labels
    col.text.low="black", col.text.high="black",
    col.text.doubleMutualist="black"
    
  ) 

# Let's define the function now

  {
    
    # We start with the first interaction matrix.
    # If empty is TRUE, we delete those rows / columns that are, well, empty...
    if (empty) 
      web <- empty(web)
    else method <- "normal"
    
    # We make sure web is formatted as a matrix object
    web <- as.matrix(web)
    meths <- c("normal", "cca")
    meths.match <- pmatch(method, meths)
    
    if (is.na(meths.match)) 
      stop("Choose plot-method: normal/cca.\n")
    
    # If method "cca" is provided, we reorder the matrix accordingly
    if (meths.match == 2) {
      ca <- cca(web)
      web <- web[order(summary(ca)$sites[, 1], decreasing = TRUE), 
                 order(summary(ca)$species[, 1], decreasing = TRUE)]
    }
    
    # If a specific sequence is provided to reorder the web, we can 
    # order the matrix accordingly
    if (!is.null(sequence)) {
      cs <- sequence$seq.pred %in% colnames(web)
      rs <- sequence$seq.prey %in% rownames(web)
      web <- web[sequence$seq.prey[rs], sequence$seq.pred[cs]]
    }
    
    # Calculate the total number of interactions in the lower web
    # I.e. the web that will be plotted at the bottom
    websum <- sum(web)
    
    # difff and diffh are vectors that calculate the difference between the frequeny
    # of species in the original web (i.e. interactions) and frequencies that could 
    # specifically be provided (e.g. independent observations) as high.abun & low.abun
    # If no specific frequency is provided, these differences are obviously 0
    difff <- diffh <- 0
    
    # Define the frequency of interactions for the "predator" in the lower web, 
    # i.e. species that will be plotted in the middle
    # Will only be executed if high.abun is provided
    # These are the double mutualists (here seen from the first web)
    # Here we will also calculate the difference between the "original" frequency of each 
    # species and their frequency in the high.abun vector
    if (!is.null(high.abun)) {
      highfreq = colSums(web)
      dummy <- highfreq
      for (i in 1:length(high.abun)) {
        ind <- which(names(high.abun)[i] == names(dummy))
        highfreq[ind] <- highfreq[ind] + high.abun[i]
      }
      diffh = (highfreq - colSums(web))/websum
    }
    
    # Define the frequency of interactions for the "prey" in the lower web, 
    # i.e. species that will be plotted at the very bottom
    # This will only be executed, if a low.abun vector is provided
    # These are the specoes that double mutualists interact with (here seen from the first web
    # Here we will also calculate the difference between the "original" frequency of each 
    # species and their frequency in the low.abun vector
    if (!is.null(low.abun)) {
      lowfreq <- rowSums(web)
      dummy <- lowfreq
      for (i in 1:length(low.abun)) {
        ind <- which(names(low.abun)[i] == names(dummy))
        lowfreq[ind] <- lowfreq[ind] + low.abun[i]
      }
      difff <- (lowfreq - rowSums(web))/websum
    }
    
    # We calculate the proportional importance of each species in the middle dimension
    # These are the double mutualists (here seen from the first web)
    # If high.abun is provided:
    if (is.null(high.abun)) 
      pred_prop <- colSums(web)/websum
    # And if not, we use what we calculated one step above
    else pred_prop <- highfreq/websum
    
    # We can also define the abundance of the double mutualist as the number of
    # individuals that we counted / captured
    # Abundance as a measure of individuals that were caught
    # Here, we use birds as double mutualists
    # birdAbun should be a vector, containing the number of counted individuals for
    # those species that appear in the network.
    # (maybe we could also include here those species where we have no interaction?)
    prop_doubleMutCounted <- birdAbun/sum(birdAbun)
    
    # We calculate the proportional importance of each species in the lower dimension
    # These are the specoes that double mutualists interact with (here seen from the first web
    if (is.null(low.abun)) 
      prey_prop <- rowSums(web)/websum
    else prey_prop <- lowfreq/websum
    
    # Number of predator (here doublemutualists) and prey (their food) species in web1
    n.pred <- length(pred_prop)
    n.prey <- length(prey_prop)
    
    # Get vectors, that contain the porportions of species in the other web
    # Vectors containing the proportions of the second web 
    # The doublemutualists
    propDoubleMutWeb2 <- rowSums(web2)/sum(web2)
    # And the species they interact with
    propIntPartWeb2 <- colSums(web2)/sum(web2)
    
    # boxes if a species is not in web1
    propMiddle <- c()
    # Now we loop through both webs to check the proportions of each respective predator
    for(pred in 1:n.pred){
      propMiddle[pred] <- ifelse(pred_prop[pred] == 0, propDoubleMutWeb2[pred], pred_prop[pred])
    }
    
    # Now we should try to standardize proportions between webs and layers to plot them nicely
    # First we define the maximum dimension of the x-axis by figuring out, the dimension with the 
    # highest species number and then calculating how far we reach on the x-axis when we add spacing 
    # in between the species-boxes
    maxXlim <- 1 + (max(nrow(web), ncol(web), ncol(web2), nrow(web2))) * 0.05
    # We need to consider the fact that we added spacing when calculating the new proportions
    # The added space would be e.g. nrow(web1) * 0.05
    # I.e. To reach the same dimensions as the longest vector we would need to consider
    # maxXlim - (length(currentVector) * 0.05)
    # Now we can standardize proportions
    
    # The interaction partner of the double mutualists in Web 1
    # Calculate the adjusted length for the vector
    adjLenIntPartWeb1 <- maxXlim - (length(prey_prop) * 0.05)
    # Standardize
    standPropIntPartWeb1 <- prey_prop / sum(prey_prop) * adjLenIntPartWeb1 
    
    # The interaction partner of the double mutualists in Web 2
    # Calculate the adjusted length for the vector
    adjLenIntPartWeb2 <- maxXlim - (length(propIntPartWeb2) * 0.05)
    # Standardize
    standPropIntPartWeb2 <- propIntPartWeb2 / sum(propIntPartWeb2) * adjLenIntPartWeb2 
    
    # The double mutualists
    # Calculate the adjusted length for the vector
    adjLenPropMiddle <- maxXlim - (length(prop_doubleMutCounted) * 0.05)
    # Standardize
    standPropMiddle <- propMiddle / sum(prop_doubleMutCounted)  * adjLenPropMiddle
    
    # Another standardized middle layer (i.e. double mutualist layer) is using
    # the proportions of counted individuals
    standPropMiddleCounted <- prop_doubleMutCounted / sum(prop_doubleMutCounted) * adjLenPropMiddle
    
    # Initialize positions and offsets for predators and prey in the plot later
    pred_x <- 0
    pred_xold <- -1
    pred_versatz <- 0
    pred_y <- 1.5
    prey_x <- 0
    prey_xold <- -1
    prey_versatz <- 0
    prey_y <- 0.5
    
    # Ensure column and row names are not NULL
    if (length(colnames(web)) == 0) 
      colnames(web) <- colnames(web, do.NULL = FALSE)
    if (length(rownames(web)) == 0) 
      rownames(web) <- rownames(web, do.NULL = FALSE)
    
    # Truncate labels if lablength is specified
    if (!is.null(lablength)) 
      colnames(web) <- substr(colnames(web), 1, lablength)
    if (!is.null(lablength)) 
      rownames(web) <- substr(rownames(web), 1, lablength)
    
    # Initiate plotting process
    par(mai = c(0.1, 0.01, 0.1, 0.01))
    
    # Define plot boundaries
    wleft <- 0
    wright <- 1 + (max(nrow(web), ncol(web), ncol(web2), nrow(web2))) * 0.05
    wup <- 2.6 + y_width + lab.space * 0.05
    wdown <- 0.4 - y_width - lab.space * 0.05
    
    # Create the plot area
    plot(0, type = "n", xlim = c(wleft, wright), 
         ylim = range(wdown / ybig, wup * ybig), axes = FALSE, 
         xlab = "", ylab = "")
    
    # Initialize positions for predator boxes
    pred_x <- 0
    hoffset <- 0
    left <- 0
    right <- 0
    height <- strheight(colnames(web)[1], cex = 0.6)
    
    # Create a vector to save the starting positions of each middle species
    # Initialize an empty vector to save pred_x values
    middle_x_values <- c()

    # Plot double mutualist species boxes from perspective of the first web
    for (i in 1:n.pred) {
      
      # Save the current value of pred_x to the middle species vector
      middle_x_values <- c(middle_x_values, pred_x)
      
      # Draw predator species box if pred_prop is greater than 0
      rect(pred_x, pred_y - y_width, 
           pred_x + standPropMiddleCounted[i],
           pred_y,
           col=ifelse(colSums(web)[i] > 0 & rowSums(web2)[i] == 0, col.doubleMutualist[i], "transparent"),
           border=ifelse(colSums(web)[i] > 0 & rowSums(web2)[i] == 0, 
                         col.border.doubleMutualist[i], "transparent"), lwd=1.2) # take out "+ y_width", to plot the predators only half as high

      width <- strwidth(colnames(web)[i], cex = 0.6 * labsize)
      left <- pred_x + standPropMiddleCounted[i] / 2 - width / 2
      
      # Adjust text offset if labels overlap
      if (left < right && i > 1) 
        hoffset <- hoffset + height
      else {
        right <- pred_x + standPropMiddleCounted[i] / 2 + width / 2
        hoffset <- 0
      }
      
      # Move to the next position for the next box
      # If a "middle" species is not present in the lower web, we should instead use
      # it's abundance in the higher web
      pred_x <- pred_x + standPropMiddleCounted[i] + 0.05
      
    }
    
    # Initialize positions for prey boxes
    prey_x <- 0
    left <- 0
    right <- 0
    height <- strheight(rownames(web)[1], cex = 0.6)
    hoffset <- height
    
    # Initialize a vector that will save the starting x coordinates of the interaction
    # partners of the doublemutualists in web1
    xPositionIntPartWeb1 <- c()
    
    # Plot prey species boxes
    for (i in 1:n.prey) {
      
      # Save the current value of pred_x to the middle species vector
      xPositionIntPartWeb1 <- c(xPositionIntPartWeb1, prey_x)
      
      # Draw prey species box
      rect(prey_x, prey_y - y_width, prey_x + standPropIntPartWeb1[i], 
           prey_y + y_width, col = col.intPartner[i], border=col.border.intPartner[i], 
           lwd=1.2)
      
      width <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
      left <- prey_x + standPropIntPartWeb1[i] / 2 - width / 2
      
      # Adjust text offset if labels overlap
      if (left < right && i > 1) 
        hoffset <- hoffset + height
      else {
        right <- prey_x + standPropIntPartWeb1[i] / 2 + width / 2
        hoffset <- height
      }
      
      # Plot prey species names
      text(prey_x + standPropIntPartWeb1[i] / 2, prey_y - y_width - hoffset, 
           rownames(web)[i], cex = 0.6 * labsize, offset = 0,
           col=col.text.low)
      
      # calculate the next x position for the next box
      prey_x <- prey_x + standPropIntPartWeb1[i] + 0.05
      
    }
    
    ############################
    ## Plot interaction lines ##
    ############################
    
    # Initialize interaction plotting variables
    pred_x <- 0
    zwischenweb <- web
    XYcoords <- matrix(ncol = 2, nrow = length(zwischenweb))
    
    # Determine the coordinates of interactions
    for (i in 1:length(zwischenweb)) {
      XYcoords[i, ] <- which(zwischenweb == max(zwischenweb), arr.ind = TRUE)[1, ]
      zwischenweb[XYcoords[i, 1], XYcoords[i, 2]] <- -1
    }
    
    # Set the y coordinates
    y1 <- pred_y - y_width
    y2 <- y1
    y3 <- prey_y + y_width
    y4 <- y3
    
    # We will need to factor in the fact, that we have standardized proportions and thus, things have slightly changed.
    # To account for the standardized scale, we will thus add the percentual change
    # In this step we need to apply the conversion to the doublemutualists in web1
    # We will calculate this case by case, as some species may not be in web1
    percentDoubleMutWeb <- numeric(length(pred_prop))
    
    for (i in seq_along(pred_prop)) {
      if (pred_prop[i] > 0) {
        percentDoubleMutWeb[i] <- (standPropMiddleCounted[i] - pred_prop[i]) / pred_prop[i]
      } else {
        percentDoubleMutWeb[i] <- 0
      }
    }
    names(percentDoubleMutWeb) <- names(pred_prop)

    # Conversion for the interaction partners of the double mutualists in web
    percentIntPartWeb1 <- (sum(standPropIntPartWeb1) - sum(prey_prop)) / sum(prey_prop)
    
    # now find the respective x coordinates
    if (sum(web > 0)) {
      
      for (p in 1:sum(web > 0)) {
        
        i <- XYcoords[p, 1]
        j <- XYcoords[p, 2]
        
        # Calculate x coordinates for interaction lines
        # First position is the position of the corresponding middle layer
        
        # Reset x1 for each new i
        x1 <- middle_x_values[j]
        
        # Adjust x1 based on the cumulative sum for each j within i
        if (i > 1) {
          x1 <- middle_x_values[j] + ((sum(web[1:(i - 1), j]) / sum(web) +
                                        ((sum(web[1:(i - 1), j]) / sum(web)) * percentDoubleMutWeb[j]))) * pollenPercent[j]
        }
        
        x2 <- x1 + ((web[i, j] / sum(web) + ((web[i, j] / sum(web)) * percentDoubleMutWeb[j]))) * pollenPercent[j]
        
        if (arrow == "up" || arrow == "both") {
          x2 <- (x1 + x2) / 2
          x1 <- x2
        }
        
        tweb <- t(web)
        
        # Reset x3 for each new i
        x3 <- xPositionIntPartWeb1[i]
        
        if (j > 1)
          x3 <- xPositionIntPartWeb1[i] + sum(tweb[1:(j - 1), i]) / sum(web) + ((sum(tweb[1:(j - 1), i]) / sum(web)) * percentIntPartWeb1)

        x4 <- x3 + tweb[j, i] / sum(web) + ((tweb[j, i] / sum(web)) * percentIntPartWeb1) 
        
        if (arrow == "down" || arrow == "both") {
          x4 <- (x3 + x4) / 2
          x3 <- x4
        }
        
        # Draw the polygon representing the interaction
        polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), 
                col = col.interaction1[j],
                border = col.borderinteraction1)
      }
      
    }
    
    #####################################################
    ## Replot box-borders on top to make it look nicer ##
    #####################################################
    
    # Define plot boundaries
    wleft <- 0
    wright <- 1 + (max(nrow(web), ncol(web), ncol(web2), nrow(web2))) * 0.05
    wup <- 2.6 + y_width + lab.space * 0.05
    wdown <- 0.4 - y_width - lab.space * 0.05
    
    # Initialize positions for predator boxes
    pred_x <- 0
    hoffset <- 0
    left <- 0
    right <- 0
    height <- strheight(colnames(web)[1], cex = 0.6)
    
    # Create a vector to save the starting positions of each middle species
    # Initialize an empty vector to save pred_x values
    middle_x_values <- c()
    
    # Remember to 
    # Plot predator species boxes
    for (i in 1:n.pred) {
      
      # Save the current value of pred_x to the middle species vector
      middle_x_values <- c(middle_x_values, pred_x)
      
      # Draw predator species box if pred_prop is greater than 0
      rect(pred_x, pred_y - y_width, 
           pred_x + standPropMiddleCounted[i],
           pred_y,
           col="transparent",
           border=ifelse(colSums(web)[i] > 0 & rowSums(web2)[i] == 0, 
                         col.border.doubleMutualist[i], "transparent"), lwd=1.2) # take out "+ y_width", to plot the predators only half as high
      
      width <- strwidth(colnames(web)[i], cex = 0.6 * labsize)
      left <- pred_x + standPropMiddleCounted[i] / 2 - width / 2
      
      # Adjust text offset if labels overlap
      if (left < right && i > 1) 
        hoffset <- hoffset + height
      else {
        right <- pred_x + standPropMiddleCounted[i] / 2 + width / 2
        hoffset <- 0
      }
      
      # Move to the next position for the next box
      # If a "middle" species is not present in the lower web, we should instead use
      # it's abundance in the higher web
      pred_x <- pred_x + standPropMiddleCounted[i] + 
        0.05
      
    }
    
    # Initialize positions for prey boxes
    prey_x <- 0
    left <- 0
    right <- 0
    height <- strheight(rownames(web)[1], cex = 0.6)
    hoffset <- height
    
    # Initialize a vector that will save the starting x coordinates of the interaction
    # partners of the doublemutualists in web1
    xPositionIntPartWeb1 <- c()
    
    # Plot prey species boxes
    for (i in 1:n.prey) {
      
      # Save the current value of pred_x to the middle species vector
      xPositionIntPartWeb1 <- c(xPositionIntPartWeb1, prey_x)
      
      # Draw prey species box
      rect(prey_x, prey_y - y_width, prey_x + standPropIntPartWeb1[i], 
           prey_y + y_width, col="transparent", border=col.border.intPartner[i], 
           lwd=1.2)
      
      width <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
      left <- prey_x + standPropIntPartWeb1[i] / 2 - width / 2
      
      # Adjust text offset if labels overlap
      if (left < right && i > 1) 
        hoffset <- hoffset + height
      else {
        right <- prey_x + standPropIntPartWeb1[i] / 2 + width / 2
        hoffset <- height
      }
      
      # Plot prey species names
      text(prey_x + standPropIntPartWeb1[i] / 2, prey_y - y_width - hoffset, 
           rownames(web)[i], cex = 0.6 * labsize, offset = 0,
           col=col.text.low)
      
      # calculate the next x position for the next box
      prey_x <- prey_x + standPropIntPartWeb1[i] + 0.05
      
    }
    
    ######################################
    #### The second interaction matrix ###
    ######################################
    web2 <- as.matrix(web2)
    
    for (i in 1:dim(web)[2]) {
      dn <- dimnames(web)[[2]][i]
      if (is.na(match(dn, dimnames(web2)[[1]]))) {
        dummy <- matrix(rep(0, dim(web2)[2]), nrow = 1)
        rownames(dummy) <- dn
        web2 <- rbind(web2, dummy)
      }
    }
    
    web2 <- web2[order(dimnames(web2)[[1]]), ]
    web2 <- web2[rank(dimnames(web)[[2]]), ]
    difff <- diffh <- 0
    dummy <- colSums(web)
    lowfreq = rowSums(web2)
    for (i in 1:length(dummy)) {
      dummy[i] <- dummy[i] - lowfreq[which(names(lowfreq) == 
                                             names(dummy[i]))]
    }
    
    # Use the lower matrix, to calculate the positions for the middle species in the 
    # upper matrix
    low.abun2 <- dummy
    lowfreq = lowfreq + low.abun2 # 
    difff = low.abun2/websum
    
    if (!is.null(high.abun2)) {
      highfreq = colSums(web2)
      dummy <- highfreq
      for (i in 1:length(high.abun2)) {
        ind <- which(names(high.abun2)[i] == names(dummy))
        highfreq[ind] <- highfreq[ind] + high.abun2[i]
      }
      diffh = (highfreq - colSums(web2))/websum
    }
    
    if (is.null(high.abun2)) 
      pred_prop <- colSums(web2)/websum
    else pred_prop <- highfreq/websum
    
    if (is.null(low.abun2)) 
      prey_prop <- rowSums(web2)/websum
    else prey_prop <- lowfreq/websum
    
    n.pred <- length(pred_prop)
    n.prey <- length(prey_prop)
    pred_x <- 0
    pred_xold <- -1
    pred_versatz <- 0
    pred_y <- 2.5
    prey_x <- 0
    prey_xold <- -1
    prey_versatz <- 0
    prey_y <- 1.5
    
    if (length(colnames(web2)) == 0) 
      colnames(web2) <- colnames(web2, do.NULL = FALSE)
    if (length(rownames(web2)) == 0) 
      rownames(web2) <- rownames(web2, do.NULL = FALSE)
    if (!is.null(lablength2)) 
      colnames(web2) <- substr(colnames(web2), 1, lablength2)
    if (!is.null(lablength2)) 
      rownames(web2) <- substr(rownames(web2), 1, lablength2)
    
    pred_x = 0
    hoffset <- 0
    left <- 0
    right <- 0
    height <- strheight(colnames(web2)[1], cex = 0.6)
    
    # Initialize a vector that stores the new x positions of the interaction partners
    # of the doublemutualists in the second web
    xPositionIntPartWeb2 <- c()
    
    # We plot the boxes for the topmost species set
    for (i in 1:n.pred) {
      
      # store the new x positions of the interaction partners now
      xPositionIntPartWeb2 <- c(xPositionIntPartWeb2, pred_x)
      
      # This colors the uppermost boxes
      rect(pred_x, pred_y - y_width, pred_x + standPropIntPartWeb2[i], 
           pred_y + y_width, col=col.intPartner2[i], border=col.border.intPartner2[i],
           lwd=1.2)
      
      width <- strwidth(colnames(web2)[i], cex = 0.6 * labsize)
      left <- pred_x + standPropIntPartWeb2[i]/2 - width/2
      
      if (left < right && i > 1) 
        hoffset = hoffset + height
      
      else {
        right <- pred_x + standPropIntPartWeb2[i]/2 + width/2
        hoffset <- 0
      }
      
      text(pred_x + standPropIntPartWeb2[i]/2, pred_y + y_width + height + 
             hoffset, colnames(web2)[i], cex = 0.6 * labsize, 
           offset = 0, col=col.text.high)
      
      pred_x <- pred_x + standPropIntPartWeb2[i] + 0.05
      
    }
    
    # And then we continue with the doublemutualists from the higher web perspective
    prey_x <- 0
    left <- 0
    right <- 0
    height <- strheight(rownames(web2)[1], cex = 0.6)
    hoffset <- height
    
    # Now we plot the middle layer of the tripartite network from the perspective of the
    # second web. E.g. double mutualists
    for (i in 1:n.prey) {
      
      # This part colors the middle boxes from the upper networks perspective
      if (rowSums(web2)[i] > 0) { # Here I use rowSums(web2) to test if a middle species is present in the second web
        if (prey_prop[i] > 0) {
          rect(middle_x_values[i], prey_y - y_width, 
               middle_x_values[i] + standPropMiddleCounted[i], 
               prey_y + y_width, col=col.doubleMutualist2[i], 
               border=col.border.doubleMutualist2[i], 
               lwd=1.2)
        }
      } 
      
      if (!is.null(low.abun2)) {
        if (prey_prop[i] == 0 ) {
          rect(middle_x_values[i] + standPropMiddleCounted[i], prey_y, 
               middle_x_values[i], prey_y + y_width,
               col=col.doubleMutualist2[i], border=col.border.doubleMutualist2[i], 
               lwd=1.2)
        }
      }
      
      width <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
      left <- prey_x + prey_prop[i]/2 - width/2
      
      if (left < right && i > 1) 
        hoffset = hoffset + height
      
      else {
        right <- prey_x + standPropMiddleCounted[i]/2 + width/2
        hoffset <- height
      }
      
      text(middle_x_values[i] + standPropMiddleCounted[i]/2, prey_y, 
           rownames(web2)[i], cex = 0.6 * labsize, offset = 0,
           col=col.text.doubleMutualist)
      
    }
    
    ##########################
    ## Add the interactions ##
    ##########################
    pred_x <- 0
    zwischenweb <- web2
    XYcoords <- matrix(ncol = 2, nrow = length(zwischenweb))
    for (i in 1:length(zwischenweb)) {
      XYcoords[i, ] <- which(zwischenweb == max(zwischenweb), 
                             arr.ind = TRUE)[1, ]
      zwischenweb[XYcoords[i, 1], XYcoords[i, 2]] <- -1
    }
    
    y1 <- pred_y - y_width
    y2 <- y1
    y3 <- prey_y + y_width
    y4 <- y3
    
    # We will need to factor in the fact, that we have standardized proportions and thus, things have slightly changed.
    # To account for the standardized scale, we will thus add the percentual change
    # Conversion for the double mutualists in web2. Again, we calculate case by case, as not all species are present in web2
    # Now what we actually plot in web 1
    percentDoubleMutWeb2 <- numeric(length(propDoubleMutWeb2))
    
    for (i in seq_along(propDoubleMutWeb2)) {
      if (propDoubleMutWeb2[i] > 0) {
        percentDoubleMutWeb2[i] <- (standPropMiddleCounted[i] - propDoubleMutWeb2[i]) / propDoubleMutWeb2[i]
      } else {
        percentDoubleMutWeb2[i] <- 0
      }
    }
    names(percentDoubleMutWeb2) <- names(propDoubleMutWeb2)
    
    # Conversion for the interaction partners of the double mutualists in web2
    percentIntPartWeb2 <- (sum(standPropIntPartWeb2) - sum(propIntPartWeb2)) / sum(propIntPartWeb2)
    
    # Let's plot the interactions between species now
    if (sum(web2 > 0)) {
      
      for (p in 1:sum(web2 > 0)) {
        i <- XYcoords[p, 1]
        j <- XYcoords[p, 2]
        
        # restart x1 after each iteration
        x1 <- xPositionIntPartWeb2[j]
        
        # Adjust x1 based on the cumulative sum for each j within i
        if (i > 1) {
          x1 <- xPositionIntPartWeb2[j] + sum(web2[1:(i - 1), j]) / sum(web2) + 
            ((sum(web2[1:(i - 1), j]) / sum(web2)) * percentIntPartWeb2)
        }
        
        x2 <- x1 + web2[i, j] / sum(web2) + ((web2[i, j] / sum(web2)) * percentIntPartWeb2) 
        
        if (arrow == "up" || arrow == "both") {
          x2 <- (x1 + x2) / 2
          x1 <- x2
        }
        
        tweb <- t(web2)
        
        # Reset x3 for each new i
        x3 <- middle_x_values[i]
        
        # Adjust x3 based on the cumulative sum for each j within i
        if (j > 1) {
          x3 <- middle_x_values[i] + ((sum(tweb[1:(j - 1), i]) / sum(tweb) +
                                         ((sum(tweb[1:(j - 1), i]) / sum(tweb)) * percentDoubleMutWeb2[i])))  * seedPercent[i]
        }
        
        x4 <- x3 + ((tweb[j, i] / sum(tweb) + 
                       ((tweb[j, i] / sum(tweb)) * percentDoubleMutWeb2[i]))) * seedPercent[i] 
        
        if (arrow2 == "down" || arrow2 == "both") {
          x4 <- (x3 + x4) / 2
          x3 <- x4
        }
        
        # Draw the polygons
        polygon(c(x1, x2, x4, x3), c(y1, y2, y4, y3), 
                col = col.interaction2[i],
                border = col.borderinteraction2)
      }
      
    }
    
    #############################
    ## replot the boxes in top ##
    #############################
    
    pred_x = 0
    hoffset <- 0
    left <- 0
    right <- 0
    height <- strheight(colnames(web2)[1], cex = 0.6)
    
    # Initialize a vector that stores the new x positions of the interaction partners
    # of the doublemutualists in the second web
    xPositionIntPartWeb2 <- c()
    
    # We plot the boxes for the topmost species set
    for (i in 1:n.pred) {
      
      # store the new x positions of the interaction partners now
      xPositionIntPartWeb2 <- c(xPositionIntPartWeb2, pred_x)
      
      # This colors the uppermost boxes
      rect(pred_x, pred_y - y_width, pred_x + standPropIntPartWeb2[i], 
           pred_y + y_width, col="transparent", border=col.border.intPartner2[i], 
           lwd=1.2)
      
      width <- strwidth(colnames(web2)[i], cex = 0.6 * labsize)
      left <- pred_x + standPropIntPartWeb2[i]/2 - width/2
      
      if (left < right && i > 1) 
        hoffset = hoffset + height
      
      else {
        right <- pred_x + standPropIntPartWeb2[i]/2 + width/2
        hoffset <- 0
      }
      
      text(pred_x + standPropIntPartWeb2[i]/2, pred_y + y_width + height + 
             hoffset, colnames(web2)[i], cex = 0.6 * labsize, 
           offset = 0, col=col.text.high)
      
      pred_x <- pred_x + standPropIntPartWeb2[i] + 0.05
      
    }
    
    # And then we continue with the doublemutualists from the higher web perspective
    prey_x <- 0
    left <- 0
    right <- 0
    height <- strheight(rownames(web2)[1], cex = 0.6)
    hoffset <- height
    
    # Now plot the middle layer of the tripartite network from the perspective of the
    # second web. E.g. double mutualists
    for (i in 1:n.prey) {
      
      # This part colors the middle boxes from the upper networks perspective
      if (rowSums(web2)[i] > 0) { # Here I use rowSums(web2) to test if a middle species is present in the second web
        if (prey_prop[i] > 0) {
          rect(middle_x_values[i], prey_y - y_width, 
               middle_x_values[i] + standPropMiddleCounted[i], 
               prey_y + y_width, 
               col="transparent", border=col.border.doubleMutualist2[i], 
               lwd=1.2)
        }
      } 
      
      if (!is.null(low.abun2)) {
        if (prey_prop[i] == 0 ) {
          rect(middle_x_values[i] + standPropMiddleCounted[i], prey_y, 
               middle_x_values[i], prey_y + y_width,
               col="transparent", border=col.border.doubleMutualist2[i], 
               lwd=1.2)
        }
      }
      
      width <- strwidth(rownames(web)[i], cex = 0.6 * labsize)
      left <- prey_x + prey_prop[i]/2 - width/2
      
      if (left < right && i > 1) 
        hoffset = hoffset + height
      
      else {
        right <- prey_x + standPropMiddleCounted[i]/2 + width/2
        hoffset <- height
      }
      
      text(middle_x_values[i] + standPropMiddleCounted[i]/2, prey_y, 
           rownames(web2)[i], cex = 0.6 * labsize, offset = 0,
           col=col.text.doubleMutualist)
      
    }
    
  }