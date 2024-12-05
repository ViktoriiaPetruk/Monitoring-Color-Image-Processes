library(doParallel)
library(foreach)

library(imager)

library(fastmatrix)

###image---

#open image
colored_flag_image <- load.image("deutsch flag color 300x300.jpg")
data_frame_color<-as.data.frame(colored_flag_image)

#extract RGB colors as vectors
color_image_vector<-data_frame_color$value
color_len<-length(color_image_vector)/3
R_image_vector<-color_image_vector[1:color_len]
G_image_vector<-color_image_vector[(color_len+1):(2*color_len)]
B_image_vector<-color_image_vector[(2*color_len+1):(3*color_len)]

#calc. cov between colors in picture
data <- data.frame(R = R_image_vector,
                   G = G_image_vector,
                   B = B_image_vector)
cov_between_colors<-cov(data)

#make R, G and B matrices
imagematrix_R<-matrix(R_image_vector,nrow=sqrt(color_len),byrow=TRUE)
imagematrix_G<-matrix(G_image_vector,nrow=sqrt(color_len),byrow=TRUE)
imagematrix_B<-matrix(B_image_vector,nrow=sqrt(color_len),byrow=TRUE)



###model----------------

##funcions--

#matrices for defining location of ROIs in image 
#(from Ivan's code)
Matrices_for_ROIS_size_1 <- function(width, height, ROIs_size){           # only for squared ROIs... so-called e_{ij}, sizes 10 and 20 non-overlap
  
  x <- list();
  
  for(j in 1:(height/ROIs_size)){
    
    for(i in 1:(width/ROIs_size)){
      
      ROIS_Matrix <- matrix(0,width,height);
      
      ROIS_Matrix[((i-1)*ROIs_size+1):(i*ROIs_size),((j-1)*ROIs_size+1):(j*ROIs_size)] <- 1;
      
      dim(ROIS_Matrix) <- NULL;                    # will be "vech'ed, columnwise
      
      x[[(j-1)*(width/ROIs_size)+i]] <- ROIS_Matrix;
      
    }
  }
  
  return(x);              
  
};


#covariance matrix of ROIs for one color (originaly - grayscale) 
#(from Ivan's code)
cov_of_means <- function(Matrices_for_ROIS_size_qua_collection, width, height, ROIs_size){       #for non-overlapping
  
  help_matrix <- matrix( c( rep(1:width,height), rep(1:height,rep(width,height)) ),width*height,2);
  #  help_matrix[278:363,]
  #  dim(help_matrix)
  
  tt <- (width*height/(ROIs_size^2));
  
  covar_mat <- matrix(0,tt,tt);
  
  sum_help_1 <- 0;
  
  var_of_pix <- 0.03^2;  #could be changed
  
  lambda <- 0.9;         #could be changed
  
  for(i in which(unlist(Matrices_for_ROIS_size_qua_collection[1])!=0)){
    for(j in which(unlist(Matrices_for_ROIS_size_qua_collection[1])!=0)){
      
      sum_help_1 <- sum_help_1+var_of_pix*lambda^sqrt( (help_matrix[i,][1]-help_matrix[j,][1])^2 +
                                                         (help_matrix[i,][2]-help_matrix[j,][2])^2);
      #      cat("i=",i,"j=",j,"",help_matrix[i,],help_matrix[j,],"\n");
    }
  }
  
  cl <- makeCluster(detectCores()-2);
  registerDoParallel(cl);
  
  covar_saver <- foreach(l=1:(tt-1), .combine="rbind", .inorder = FALSE) %dopar% {
    
    var_of_pix <- 0.03^2;  #could be changed
    
    lambda <- 0.9;         #could be changed
    
    for(m in (l+1):(tt)){
      
      sum_help_2 <- 0;
      
      for(i in which(unlist(Matrices_for_ROIS_size_qua_collection[l])!=0)){
        for(j in which(unlist(Matrices_for_ROIS_size_qua_collection[m])!=0)){
          
          sum_help_2<-sum_help_2+var_of_pix*lambda^sqrt((help_matrix[i,][1]-help_matrix[j,][1])^2 +
                                                          (help_matrix[i,][2]-help_matrix[j,][2])^2);
        }
        
      }
      
      covar_mat[l,m] <- 1/((ROIs_size^2)^2)*sum_help_2;
      
    }
    
    covar_mat[l,]
    
  }
  
  stopCluster(cl);
  
  dimnames(covar_saver) <- NULL;
  
  cat(dim(covar_saver));
  
  covar_mat <- rbind(covar_saver,rep(0,tt))+t(rbind(covar_saver,rep(0,tt)));
  
  diag(covar_mat) <- 1/((ROIs_size^2)^2)*sum_help_1;
  
  return(covar_mat);
  
}

#covariance matrix of ROIs in RGB case where order of ROIs is: 
#(r g b  r g b  r g b  r g b...)x(r g b  r g b  r g b  r g b...)
#instead of one element there is this element multiplied by 3x3 matrix cov_of_colors
cov_of_ROIs_color <- function(width, height, matrix_for_one_dim, cov_matrix_of_RGB, ROIs_size){
  
  help_matrix<- matrix(0, (width*height/(ROIs_size^2))*3, (width*height/(ROIs_size^2))*3)
  
  for(i in 1:(width*height/(ROIs_size^2))){
    for (j in 1:(width*height/(ROIs_size^2))){
      help_matrix[(1+3*(i-1)):(3*i), (1+3*(j-1)):(3*j)]<-matrix_for_one_dim[i,j]*cov_matrix_of_RGB
    }
  }
  return(help_matrix)
}


#vector of values of ROIs (calculate for each color separately) from imagematrix_A=Book_cover 
#(from Ivan's code)
Mean_for_ROIs <- function(Matrices_for_ROIS_size_qua_collection, Book_cover, ROIs_size){
  
  Mean_im <- numeric();
  
  for(pp in 1:length(Matrices_for_ROIS_size_qua_collection)){
    
    Mean_im[pp] <- sum(unlist(Matrices_for_ROIS_size_qua_collection[pp])*Book_cover)/(ROIs_size^2);
    
  }
  
  return(Mean_im);
  
}

#organize vectors of ROIs of each color into one vector: [r g b  r g b  r g b  r g b...]
Mean_for_ROIs_color<- function(height, width, means_R, means_G, means_B){
  
  help_vector<-numeric();
  for(i in 1:length(means_R)){
    help_vector[(1+3*(i-1))]<-means_R[i]
    help_vector[(2+3*(i-1))]<-means_G[i]
    help_vector[(3*i)]<-means_B[i]
    
  }
  
  return(help_vector)
  
}


#covariance matrix for grayscale picture
cov_from_color_to_gray<-function(cov_of_ROIs_RGB){
  #> dim(cov_of_ROIs_RGB) [1] 675 675
  dimention_of_new_m<-dim(cov_of_ROIs_RGB)[1]/3
  cov_of_grays<-matrix(0,nrow = dimention_of_new_m, ncol = dimention_of_new_m)
  
  for(i in 1:dimention_of_new_m){
    for (j in 1:dimention_of_new_m){
      cov_of_grays[i,j]<-mean(cov_of_ROIs_RGB[((i-1)*3+1):((i-1)*3+3),((j-1)*3+1):((j-1)*3+3)])
    }
  }
  
  return(cov_of_grays)
  
}



#grayscale function for a RGB vector: [r g b  r g b  r g b  r g b...]
make_grayscale<-function(color_vector){
  gray_vector<-c()
  for(i in 1:(length(color_vector)/3)){
    gray_vector[i]<-mean(c(color_vector[(i-1)*3+1],color_vector[(i-1)*3+2],color_vector[(i-1)*3+3]))
  }
  return(gray_vector)
}
#save(make_grayscale, file = "make_grayscale_function.Rdata")


##calculations--

ROIs_size_2 <- 20

width <- sqrt(color_len)
height <- sqrt(color_len)

cov_matrix_of_RGB<-cov_between_colors

Matrices_for_ROIS_size_20x20 <- Matrices_for_ROIS_size_1(width, height, 
                                                         ROIs_size=ROIs_size_2)

#takes time
cov_of_means_ROI<- cov_of_means(Matrices_for_ROIS_size_20x20,
                               width,height,ROIs_size=ROIs_size_2)
#save(cov_of_means_ROI, file = "cov_of_means_ROI_for_300x300.Rdata")
#load("cov_of_means_ROI_for_300x300.Rdata")

cov_of_ROIs_RGB<-cov_of_ROIs_color(width, height, cov_of_means_ROI, 
                                   cov_matrix_of_RGB, ROIs_size=ROIs_size_2)

chol_decomp_of_cov_of_means_2 <- t(chol(cov_of_ROIs_RGB))
trace_cov_for_20 <- sum(diag(cov_of_ROIs_RGB));
sq_trace_cov_for_20 <- sqrt(2*sum(diag(cov_of_ROIs_RGB%*%cov_of_ROIs_RGB)))

inv_cov_matrix_RGB<-solve(cov_of_ROIs_RGB)

#save(inv_cov_matrix_RGB, file = "inv_cov_matrix_RGB.Rdata")
#save(chol_decomp_of_cov_of_means_2, file = "chol_decomp_of_cov_of_means_2_300.Rdata")
#save(sq_trace_cov_for_20, file = "sq_trace_cov_for_20_300.Rdata")
#save(trace_cov_for_20, file = "trace_cov_for_20_300.Rdata")




cov_of_means_ROI_gray<-cov_from_color_to_gray(cov_of_ROIs_RGB)

chol_decomp_of_cov_of_means_gray <- t( chol(cov_of_means_ROI_gray) );
trace_cov_for_20_gray <- sum( diag( cov_of_means_ROI_gray ) );
sq_trace_cov_for_20_gray <- sqrt(2*sum(diag(cov_of_means_ROI_gray%*%cov_of_means_ROI_gray)));

inv_cov_matrix_gray<-solve(cov_of_means_ROI_gray)





##alternative for cov of grayscale?

I_3<-t(t(rep(1, 3)))
#I_r<-t(t(rep(1, r_gray)))
I_r<-diag(1, r_gray, r_gray )

kronec<-kronecker.prod(I_r, t(I_3))

alternative_cov_of_gray<-(kronec%*%cov_of_ROIs_RGB%*%kronecker.prod(I_r, I_3))/9 #works!


alternative_mean_im_gray<-(kronec%*%Mean_im_RGB)/3 #the same =)

alt_inv_cov_matrix_gray<-solve(alternative_cov_of_gray)


alt_trace_cov_for_20_gray <- sum( diag( alternative_cov_of_gray ) )
alt_sq_trace_cov_for_20_gray <- sqrt(2*sum(diag(alternative_cov_of_gray%*%alternative_cov_of_gray)))


#save(alternative_cov_of_gray, file = "alternative_cov_of_gray.Rdata") 
#save(alt_inv_cov_matrix_gray, file = "alt_inv_cov_matrix_gray.Rdata")
#save(alt_trace_cov_for_20_gray, file = "alt_trace_cov_for_20_gray.Rdata") 
#save(alt_sq_trace_cov_for_20_gray, file = "alt_sq_trace_cov_for_20_gray.Rdata")



#save(inv_cov_matrix_gray, file = "inv_cov_matrix_gray.Rdata") 
#save(trace_cov_for_20_gray, file = "trace_cov_for_20_gray.Rdata")
#save(sq_trace_cov_for_20_gray, file = "sq_trace_cov_for_20_gray.Rdata")

Mean_im_R <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_R, ROIs_size_2);
Mean_im_G <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_G, ROIs_size_2)
Mean_im_B <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_B, ROIs_size_2)

Mean_im_RGB<-Mean_for_ROIs_color(height, width, Mean_im_R, Mean_im_G, Mean_im_B)
#save(Mean_im_RGB, file = "Mean_im_RGB_300.Rdata")







###changes (not all changes yet)----------------

#change in whole picture
flag_as_hsi<- RGBtoHSI(colored_flag_image)
flag_as_hsi_changed_hue<-imchange(flag_as_hsi,~ c==1,~ .+(1.5)) #for hue change
#flag_as_hsi_changed<-imchange(flag_as_hsi,~ c==2,~ ./(1.01) ) #for saturation change
colored_flag_changed<-HSItoRGB(flag_as_hsi_changed_hue)

#further calculations for changed picture
data_frame_color_changed<-as.data.frame(colored_flag_changed)

color_image_vector_changed<-data_frame_color_changed$value

R_image_vector_changed<-color_image_vector_changed[1:color_len]
G_image_vector_changed<-color_image_vector_changed[(color_len+1):(2*color_len)]
B_image_vector_changed<-color_image_vector_changed[(2*color_len+1):(3*color_len)]


imagematrix_R_changed<-matrix(R_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)
imagematrix_G_changed<-matrix(G_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)
imagematrix_B_changed<-matrix(B_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)

Mean_im_R_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_R_changed, ROIs_size_2)
Mean_im_G_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_G_changed, ROIs_size_2)
Mean_im_B_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_B_changed, ROIs_size_2)

Mean_im_RGB_changed<-Mean_for_ROIs_color(height, width, Mean_im_R_changed, Mean_im_G_changed, Mean_im_B_changed)
#save(Mean_im_RGB_changed, file = "Mean_im_RGB_changed_hue.Rdata")

#or--------
#making change in intensity of one of the colors in the rectangular area
change_in_mean<-function(matrix_to_change, from_i, from_j, size_i, size_j, 
                         value_of_change ){
  
  matrix_to_change[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]<-
    matrix_to_change[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]+value_of_change
  
  return(matrix_to_change)
}

imagematrix_R_changed<-change_in_mean(imagematrix_R, 140, 60, 20, 180, 
                                      0.003 )

#further calculations for imagematrix_R_changed
Mean_im_R_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_R_changed, ROIs_size_2)

Mean_im_RGB_changed<-Mean_for_ROIs_color(height, width, Mean_im_R_changed, Mean_im_G, Mean_im_B)
#save(Mean_im_RGB_changed, file = "Mean_im_RGB_changed_red.Rdata")


#green
imagematrix_G_changed<-change_in_mean(imagematrix_G, 140, 60, 20, 180, 
                                      0.003 )

#further calculations for imagematrix_R_changed
Mean_im_G_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_G_changed, ROIs_size_2)

Mean_im_RGB_changed<-Mean_for_ROIs_color(height, width, Mean_im_R, Mean_im_G_changed, Mean_im_B)
#save(Mean_im_RGB_changed, file = "Mean_im_RGB_changed_green.Rdata")

#blue
imagematrix_B_changed<-change_in_mean(imagematrix_B, 140, 60, 20, 180, 
                                      0.003 )

#further calculations for imagematrix_R_changed
Mean_im_B_changed <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_B_changed, ROIs_size_2)

Mean_im_RGB_changed<-Mean_for_ROIs_color(height, width, Mean_im_R, Mean_im_G, Mean_im_B_changed)
#save(Mean_im_RGB_changed, file = "Mean_im_RGB_changed_blue.Rdata")


#or------
#it can be easier to first change the whole picture and then make a picture with only its part changed
replace_part_for_changed_image<-function(data_frame_color, data_frame_color_changed,
                                         from_i, from_j, size_i, size_j,
                                         Mean_for_ROIs, Mean_for_ROIs_color,
                                         ROIs_size_2,height, width){
  
  color_image_vector<-data_frame_color$value
  color_len<-length(color_image_vector)/3
  R_image_vector<-color_image_vector[1:color_len]
  G_image_vector<-color_image_vector[(color_len+1):(2*color_len)]
  B_image_vector<-color_image_vector[(2*color_len+1):(3*color_len)]
  
  imagematrix_R<-matrix(R_image_vector,nrow=sqrt(color_len),byrow=TRUE)
  imagematrix_G<-matrix(G_image_vector,nrow=sqrt(color_len),byrow=TRUE)
  imagematrix_B<-matrix(B_image_vector,nrow=sqrt(color_len),byrow=TRUE)
  
  imagematrix_R_new<-imagematrix_R
  imagematrix_G_new<-imagematrix_G
  imagematrix_B_new<-imagematrix_B
  
  
  color_image_vector_changed<-data_frame_color_changed$value
  
  R_image_vector_changed<-color_image_vector_changed[1:color_len]
  G_image_vector_changed<-color_image_vector_changed[(color_len+1):(2*color_len)]
  B_image_vector_changed<-color_image_vector_changed[(2*color_len+1):(3*color_len)]
  
  
  imagematrix_R_changed<-matrix(R_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)
  imagematrix_G_changed<-matrix(G_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)
  imagematrix_B_changed<-matrix(B_image_vector_changed,nrow=sqrt(color_len),byrow=TRUE)
  
  imagematrix_R_new[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]<-
    imagematrix_R_changed[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]
  
  imagematrix_G_new[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]<-
    imagematrix_G_changed[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]
  
  imagematrix_B_new[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]<-
    imagematrix_B_changed[from_i:(from_i+(size_i-1)),from_j:(from_j+(size_j-1))]
  
  
  Mean_im_R_new <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_R_new, ROIs_size_2)
  Mean_im_G_new <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_G_new, ROIs_size_2)
  Mean_im_B_new <- Mean_for_ROIs(Matrices_for_ROIS_size_20x20, imagematrix_B_new, ROIs_size_2)
  
  Mean_im_RGB_new<-Mean_for_ROIs_color(height, width, Mean_im_R_new, Mean_im_G_new, Mean_im_B_new)
  
  
  return(Mean_im_RGB_new)
}


#change in whole picture
flag_as_hsi<- RGBtoHSI(colored_flag_image)
flag_as_hsi_changed_hue<-imchange(flag_as_hsi,~ c==1,~ .+(1.5) ) #for hue change .....*exp(-(rho/550)^2), 25 for illustration?
#flag_as_hsi_changed_hue<-imchange(flag_as_hsi_changed_hue,~ c==1,~ .%%(360) )
#flag_as_hsi_changed<-imchange(flag_as_hsi,~ c==2,~ ./(1.01) ) #for saturation change
colored_flag_changed<-HSItoRGB(flag_as_hsi_changed_hue)
plot(colored_flag_changed)
data_frame_color_changed<-as.data.frame(colored_flag_changed)

Mean_im_RGB_changed<-replace_part_for_changed_image(data_frame_color, data_frame_color_changed,
                                                    140, 60, 20, 180,
                                                    Mean_for_ROIs, Mean_for_ROIs_color,
                                                    ROIs_size_2,height, width)


#save(Mean_im_RGB_changed, file = "Mean_im_RGB_changed_hue(_1_5)_stripe.Rdata")

