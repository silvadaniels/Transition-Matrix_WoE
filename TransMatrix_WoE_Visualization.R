# LCM visualization based on the script 'TransMatrix_WoE_function.R'
  library(ggplot2)
  library("markovchain")
  library(diagram)

# plot map of layers
  plot(DEM, col=topo.colors(6), breaks=c(100, 150, 200, 250, 300, 350), main="Digital Elevation Model")
  plot(lots, add=T)
  plot(roads, add=T, col="red")
  plot(hydro, add=T, col="blue")
  #plot(terrain, col=topo.colors(6))

# 1. Markov chain transition probabilities ####

  # 1.1. Plot transition prob scheme ####
  markov2 <- new("markovchain",
                 transitionMatrix = tp,
                 states = c('veg','agric','urban', 'water'))
  
    row.names(tp) <- c("Vegetation","Agriculture","Urban","Water"); colnames(tp) <- c("Vegetation","Agriculture","Urban","Water")
    plotmat(round(tp,2),pos = c(2,2), lwd = 1, box.lwd = 0, cex.txt = 0.8, shadow.size = 0,
          box.size = 0.1, box.type = "circle", box.prop = 0.5, box.col = "#ffd966",
          arr.length=.1, arr.width=.1, self.cex = .6, self.shifty = -.01, self.shiftx = .14,
          main = "Markov Chain")

  # 1.2. Plot markov chain probabilities for several steps ####
    # Initiate probability vectors
      initState <- c(1,0,0,0) # The initial state - describe the probability of the chain begining at the state I (i.e. from veg to any other LU)
      
      Veg <- c()
      Agr <- c()
      Urban <- c()
      Water = c()
  
    # Calculate probabilities for 10 steps.
      for(k in 1:10){
        nsteps <- initState*markov2^k
        Veg[k] <- nsteps[1,1]
        Agr[k] <- nsteps[1,2]
        Urban[k] <- nsteps[1,3]
        Water[k] = nsteps[1,4]
        }
  
    # Make dataframes and merge them
      Veg <- as.data.frame(Veg)
      Veg$Group <- 'Veg'
      Veg$Iter <- 1:10
      names(Veg)[1] <- 'Value'
      
      Agr <- as.data.frame(Agr)
      Agr$Group <- 'Agr'
      Agr$Iter <- 1:10
      names(Agr)[1] <- 'Value'
      
      Urban <- as.data.frame(Urban)
      Urban$Group <- 'Urban'
      Urban$Iter <- 1:10
      names(Urban)[1] <- 'Value'
      
      Water <- as.data.frame(Water)
      Water$Group <- 'Water'
      Water$Iter <- 1:10
      names(Water)[1] <- 'Value'
      
      steps <- rbind(Veg,Agr,Urban,Water)
      
    # Plot the projected probabilities (by default this chart considers veg as initial state; to change it, go back to initial prob vector)
      ggplot(steps, aes(x = Iter, y = Value, col = Group))+
        geom_line()+ xlab('Chain Step') +ylab('Probability') 
        +ggtitle('10 Step Chain Probability Prediction') +theme(plot.title = element_text(hjust = 0.5))

# 2. WoE and land change (LC) probabilities ####

  # WoE and posterior prob
    plot_infotables(IV, IV$Summary$Variable[c(1,3,4,6)], same_scale=FALSE)

  # Map (raster) of posterior probabilities
    plot(Dp, col=heat.colors(4, rev = T))
