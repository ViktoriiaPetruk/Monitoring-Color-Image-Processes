library(foreach)
library(doParallel)


#remade from Ivan's code MRL_simul_180_300_Main_stat_U_n_out_of_control
ARL_Main_stat_U_n_out_of_control_both <- function(CL, CL_gray){    
  
  
  width<-300
  height<-300
  ROIs_size_2<-20
  load(file = "Mean_im_RGB_300.Rdata")
  load(file = "chol_decomp_of_cov_of_means_2_300.Rdata")
  load(file = "sq_trace_cov_for_20_300.Rdata")
  
  load(file = "sq_trace_cov_for_20_gray.Rdata")
  load(file = "make_grayscale_function.Rdata")
  
  
  #load(file = "Mean_im_RGB_changed_hue(_2)_stripe.Rdata")
  #load(file = "Mean_im_RGB_changed_red.Rdata")
  #load(file = "Mean_im_RGB_changed_blue.Rdata")
  load(file = "Mean_im_RGB_changed_hue.Rdata")
  
  
  
  Mean_im_gray<-make_grayscale(Mean_im_RGB)
  
  counts_sim <- foreach(n=1:350, .combine="rbind", .inorder = FALSE) %dopar% {
    
    r <- 1;
    r_gray <- 1 
    
    steps<-2000
    xr <- matrix(0, 3*width*height/(ROIs_size_2^2), (steps+20))
    
    xr_gray <- matrix(0, width*height/(ROIs_size_2^2), (steps+20))
    
    

    # 20 runs in-control
    for(j in 1:(20)){

      ROI <- chol_decomp_of_cov_of_means_2%*%rnorm( 3*width*height/(ROIs_size_2^2) );
      dim(ROI) <- NULL;
      xr[,j] <- ROI + Mean_im_RGB;
      xr_gray[,j] <- make_grayscale(ROI + Mean_im_RGB)

      # some_help <- numeric(0);
      # 
      # for (tau in 1:(j-1)) {
      # 
      #   mat_h <- t( xr[,j:tau,drop=FALSE]-Mean_im_RGB )%*%( xr[,j:tau,drop=FALSE]-Mean_im_RGB);
      #   some_help[tau] <- 1/(sqrt(2*(j-tau+1)*(j-tau)))*( sum( mat_h ) - sum( diag(mat_h) ) )/sq_trace_cov_for_20*sqrt(2);
      # }
      # Result <- max(some_help)


    }
    
    color<-TRUE
    gray<-TRUE
    
    for(j in 21:(steps+20)){
      
      ROI <- chol_decomp_of_cov_of_means_2%*%rnorm( 3*width*height/(ROIs_size_2^2) );
      dim(ROI) <- NULL;
      xr[,j] <- ROI + Mean_im_RGB_changed;  #Mean_im_RGB for in-control test 
      xr_gray[,j] <- make_grayscale(ROI + Mean_im_RGB_changed) #Mean_im_RGB_changed for out-of-control test
      
      
      if(color){
        some_help <- numeric(0);
        
        for (tau in 1:(j-1)) {
          
          mat_h <- t( xr[,j:tau,drop=FALSE]-Mean_im_RGB )%*%( xr[,j:tau,drop=FALSE]-Mean_im_RGB);
          some_help[tau] <- 1/(sqrt(2*(j-tau+1)*(j-tau)))*( sum( mat_h ) - sum( diag(mat_h) ) )/sq_trace_cov_for_20*sqrt(2);
        }
        Result <- max(some_help)
        if(Result < CL){r <- r+1} else {color<-FALSE}
      }
      
      if(gray){
        some_help <- numeric(0);
        
        for (tau in 1:(j-1)) {
          
          mat_h <- t( xr_gray[,j:tau,drop=FALSE]-Mean_im_gray )%*%( xr_gray[,j:tau,drop=FALSE]-Mean_im_gray);
          some_help[tau]<-(1/(sqrt(2*(j-tau+1)*(j-tau))))*(sum(mat_h)-sum(diag(mat_h)))/sq_trace_cov_for_20_gray*sqrt(2) #why *sqrt(2) ?
          
        }
        Result <- max(some_help)
        if(Result < CL_gray){r_gray <- r_gray+1} else {gray<-FALSE}
      }
      
      if((!color)&(!gray)){break}
      
      
    }
    res<-c(r,r_gray)
    return(res)
    }
  
  
  
  #return( list=c( median(counts_sim[,1]), mean(counts_sim[,1]), sd(counts_sim[,1]) ) )
  return(list=c(mean(counts_sim[,1]),mean(counts_sim[,2])) )
  
}

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

sim_for_ARL_with_change_both<-ARL_Main_stat_U_n_out_of_control_both(3.1, 3.129)
sim_for_ARL_with_change_both

sim_for_ARL_with_change_both<-ARL_Main_stat_U_n_out_of_control_both(3.1, 3.129)
sim_for_ARL_with_change_both

stopCluster(cl)






