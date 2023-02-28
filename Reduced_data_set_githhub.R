install.packages("scales")
library(scales)
##################  Reduced data set ####
# Data with first and last marker shared by all individuals ####
All_data <- read.table("all_data.txt", header= T) #four columns with "Parent","LG","Estimated_physicalposition_bp","Genetical_position_cM"
All_data <- cbind(All_data, NA)

LG1<-subset(All_data, All_data$LG == "1" )
LG2<-subset(All_data, All_data$LG == "2" )
LG3<-subset(All_data, All_data$LG == "3" )
LG4<-subset(All_data, All_data$LG == "4" )
LG5<-subset(All_data, All_data$LG == "5" )
LG6<-subset(All_data, All_data$LG == "6" )
LG7<-subset(All_data, All_data$LG == "7" )
LG8<-subset(All_data, All_data$LG == "8" )
LG9<-subset(All_data, All_data$LG == "9" )
LG10<-subset(All_data, All_data$LG == "10" )
LG11<-subset(All_data, All_data$LG == "11" )
LG12<-subset(All_data, All_data$LG == "12" )

df_list <- list (LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10,LG11,LG12)
finaldata <- rep(list(matrix(NA)), length(df_list))
for (k in 1:length(df_list)) { 
  print(k)
  LG11_A<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_TEX114" )
  LG11_B<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_TEX1")
  LG11_C<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_LPB87" )
  LG11_D<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_STM2")
  end <- min(c(LG11_A[which.max(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.max(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.max(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.max(LG11_D$Estimated_physicalposition_bp.),3]))
  debut <- max(c(LG11_A[which.min(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.min(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.min(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.min(LG11_D$Estimated_physicalposition_bp.),3]))
  #end <- LG12[which.max(LG12$Estimated_physicalposition_bp.),3]
  LG <- as.character(LG11_A[1,2])
  if( all(LG11_A$Estimated_physicalposition_bp. != debut )) {
    LG11_A<-rbind(LG11_A, c("Female_TEX114", LG,debut,NA, NA))}
  if( all(LG11_A$Estimated_physicalposition_bp. != end )) {
    LG11_A<-rbind(LG11_A, c("Female_TEX114", LG,end,NA, NA))  }
  LG11_A$Estimated_physicalposition_bp.<- as.numeric(LG11_A$Estimated_physicalposition_bp.)
  LG11_A <- LG11_A[order(LG11_A$Estimated_physicalposition_bp.),]
  LG11_A$Genetical_position_cM. <- as.numeric(LG11_A$Genetical_position_cM.)
  if( all(LG11_B$Estimated_physicalposition_bp. != debut )) {
    LG11_B<-rbind(LG11_B, c("Male_TEX1", LG,debut,NA, NA))}
  if( all(LG11_B$Estimated_physicalposition_bp. != end )) {
    LG11_B<-rbind(LG11_B, c("Male_TEX1", LG,end,NA, NA))}
  LG11_B$Estimated_physicalposition_bp. <- as.numeric(LG11_B$Estimated_physicalposition_bp.)
  LG11_B <- LG11_B[order(LG11_B$Estimated_physicalposition_bp.),]
  LG11_B$Genetical_position_cM. <- as.numeric(LG11_B$Genetical_position_cM.)
  if( all(LG11_C$Estimated_physicalposition_bp. != debut )) {
    LG11_C<-rbind(LG11_C, c("Female_LPB87", LG,debut,NA, NA))}
  if( all(LG11_C$Estimated_physicalposition_bp. != end )) {
    LG11_C<-rbind(LG11_C, c("Female_LPB87", LG,end,NA, NA))}
  LG11_C$Estimated_physicalposition_bp. <- as.numeric(LG11_C$Estimated_physicalposition_bp.)
  LG11_C <- LG11_C[order(LG11_C$Estimated_physicalposition_bp.),]
  LG11_C$Genetical_position_cM. <- as.numeric(LG11_C$Genetical_position_cM.)
  if( all(LG11_D$Estimated_physicalposition_bp. != debut )) {
    LG11_D<-rbind(LG11_D, c("Male_STM2", LG,debut,NA, NA))}
  if( all(LG11_D$Estimated_physicalposition_bp. != end )) {
    LG11_D<-rbind(LG11_D, c("Male_STM2", LG,end,NA, NA))}
  LG11_D$Estimated_physicalposition_bp. <- as.numeric(LG11_D$Estimated_physicalposition_bp.)
  LG11_D <- LG11_D[order(LG11_D$Estimated_physicalposition_bp.),]
  LG11_D$Genetical_position_cM. <- as.numeric(LG11_D$Genetical_position_cM.)
  
# New cM with the beginning and end, approx linearly
  
  for (i in 1:nrow(LG11_A)) {
    if( LG11_A[i,3] != debut & LG11_A[i,3] != end){ 
      LG11_A[i,5] <- LG11_A[i,4]
    }  else if(i != 1 & LG11_A[i,3] == debut) {
      T2_A <- data.frame(bp=c(LG11_A[i-1,3],LG11_A[i+1,3]), cm=(c(LG11_A[i-1,4],LG11_A[i+1,4])))
      T1_A <-data.frame(bp = debut)
      LG11_A[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    } else if( i != nrow(LG11_A) & LG11_A[i,3] == end){
      T2_A <- data.frame(bp=c(LG11_A[i-1,3],LG11_A[i+1,3]), cm=(c(LG11_A[i-1,4],LG11_A[i+1,4])))
      T1_A <-data.frame(bp = end)
      LG11_A[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    }  else {
      LG11_A[i,5] <- LG11_A[i,4]  
    }
  }
  
  for (i in 1:nrow(LG11_B)) {
    if( LG11_B[i,3] != debut & LG11_B[i,3] != end){ 
      LG11_B[i,5] <- LG11_B[i,4]
    }  else if( i != 1 & LG11_B[i,3] == debut) {
      T2_A <- data.frame(bp=c(LG11_B[i-1,3],LG11_B[i+1,3]), cm=(c(LG11_B[i-1,4],LG11_B[i+1,4])))
      T1_A <-data.frame(bp = debut)
      LG11_B[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    } else if( i != nrow(LG11_B) & LG11_B[i,3] == end){
      T2_A <- data.frame(bp=c(LG11_B[i-1,3],LG11_B[i+1,3]), cm=(c(LG11_B[i-1,4],LG11_B[i+1,4])))
      T1_A <-data.frame(bp = end)
      LG11_B[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    }  else {
      LG11_B[i,5] <- LG11_B[i,4]  
    }
  }
  for (i in 1:nrow(LG11_C)) {
    if( LG11_C[i,3] != debut & LG11_C[i,3] != end){ 
      LG11_C[i,5] <- LG11_C[i,4]
    }  else if( i != 1 &LG11_C[i,3] == debut) {
      T2_A <- data.frame(bp=c(LG11_C[i-1,3],LG11_C[i+1,3]), cm=(c(LG11_C[i-1,4],LG11_C[i+1,4])))
      T1_A <-data.frame(bp = debut)
      LG11_C[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    } else if( i != nrow(LG11_C) & LG11_C[i,3] == end){
      T2_A <- data.frame(bp=c(LG11_C[i-1,3],LG11_C[i+1,3]), cm=(c(LG11_C[i-1,4],LG11_C[i+1,4])))
      T1_A <-data.frame(bp = end)
      LG11_C[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    }   else {
      LG11_C[i,5] <- LG11_C[i,4]  
    }
  }
  for (i in 1:nrow(LG11_D)) {
    if( LG11_D[i,3] != debut & LG11_D[i,3] != end){ 
      LG11_D[i,5] <- LG11_D[i,4]
    }  else if(i != 1 & LG11_D[i,3] == debut) {
      T2_A <- data.frame(bp=c(LG11_D[i-1,3],LG11_D[i+1,3]), cm=(c(LG11_D[i-1,4],LG11_D[i+1,4])))
      T1_A <-data.frame(bp = debut)
      LG11_D[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    } else if( i != nrow(LG11_D) & LG11_D[i,3] == end){
      T2_A <- data.frame(bp=c(LG11_D[i-1,3],LG11_D[i+1,3]), cm=(c(LG11_D[i-1,4],LG11_D[i+1,4])))
      T1_A <-data.frame(bp = end)
      LG11_D[i,5] <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y
    }   else {
      LG11_D[i,5] <- LG11_D[i,4]  
    }
  }
  LG11_A_crop <- LG11_A[LG11_A$Estimated_physicalposition_bp.>=debut & LG11_A$Estimated_physicalposition_bp.<= end,]
  LG11_B_crop <- LG11_B[LG11_B$Estimated_physicalposition_bp.>=debut & LG11_B$Estimated_physicalposition_bp.<= end,]
  LG11_C_crop <- LG11_C[LG11_C$Estimated_physicalposition_bp.>=debut & LG11_C$Estimated_physicalposition_bp.<= end,]
  LG11_D_crop <- LG11_D[LG11_D$Estimated_physicalposition_bp.>=debut & LG11_D$Estimated_physicalposition_bp.<= end,]
  finaldata[[k]] <-rbind(LG11_A_crop, LG11_B_crop, LG11_C_crop,LG11_D_crop)
}

reduced_data <- as.data.frame(do.call(rbind, finaldata))

write.table(reduced_data, file = "data_for_comparison_analysis.txt", sep = "\t", row.names = FALSE)

#### WARNING be carefull this data set is ok but we then put each first position of each LG at 0 !!


######### ALL boxplots together : ####
library(reshape2)
library(Rcpp)
setwd("D:/Maxtor/Cecile_Molinier/RAD_seq_cecile/2_librairies_controles/Results/Visualisation/Genetic_map/")
#counts = read.delim("/media/loukesio/Maxtor/Cecile_Molinier/RAD_seq_cecile/2_librairies_controles/Results/counts_controle.txt")
counts = read.delim("D:/Maxtor/Cecile_Molinier/RAD_seq_cecile/2_librairies_controles/Results/counts_controle_reduced.txt", dec=",", stringsAsFactors = F)
#counts_2 = read.delim("/media/loukesio/Maxtor/Cecile_Molinier/RAD_seq_cecile/2_librairies_controles/Results/counts_asexvssex.txt")
counts_2 = read.delim("D:/Maxtor/Cecile_Molinier/RAD_seq_cecile/2_librairies_controles/Results/counts_asexvssex_reduced.txt", dec=",", stringsAsFactors = F)
counts_2 <- subset(counts_2[1:12,c(1,2,4)])
colnames(counts_2) <- c("chrom","CP_female_TEX-114","OP_male_STM-2")
colnames(counts) <- c("chrom","CP_male_TEX-1","CP_female_LPB-87")


  ############################# localization of the CO  #########
  # We will use the data scaled.
  
  ## Determine the bin size fixed for all : 5 bin for each LG
  library(reshape2)
  phys = read.delim("reduced_physical_length.txt", stringsAsFactors = F)
  bin_size <- phys$phys_length/15
  phys <- cbind(phys, bin_size)
  
  #the intial data is not reduced for the terminals markers but this is done in the script
  All_data <- read.table("Data_scaled.txt", header= T)
  
  colnames(All_data) <- c("Parent","LG","Estimated_physicalposition_bp.","Genetical_position_cM_scaled"   )
  LG1<-subset(All_data, All_data$LG == "1" )
  LG2<-subset(All_data, All_data$LG == "2" )
  LG3<-subset(All_data, All_data$LG == "3" )
  LG4<-subset(All_data, All_data$LG == "4" )
  LG5<-subset(All_data, All_data$LG == "5" )
  LG6<-subset(All_data, All_data$LG == "6" )
  LG7<-subset(All_data, All_data$LG == "7" )
  LG8<-subset(All_data, All_data$LG == "8" )
  LG9<-subset(All_data, All_data$LG == "9" )
  LG10<-subset(All_data, All_data$LG == "10" )
  LG11<-subset(All_data, All_data$LG == "11" )
  LG12<-subset(All_data, All_data$LG == "12" )
  
  df_list <- list (LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10,LG11,LG12)
  finaldata <- rep(list(matrix(NA)), length(df_list))
  for (k in 1:length(df_list)) { 
    print(k)
    LG11_A<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_TEX114" )
    LG11_B<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_TEX1")
    LG11_C<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_LPB87" )
    LG11_D<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_STM2")
    end <- min(c(LG11_A[which.max(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.max(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.max(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.max(LG11_D$Estimated_physicalposition_bp.),3]))
    debut <- max(c(LG11_A[which.min(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.min(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.min(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.min(LG11_D$Estimated_physicalposition_bp.),3]))
    wind <- phys[k,3]
    ##### LG sifnificantly lower than other LG
    Boot <- matrix(NA, nrow = 1, ncol = 4)
    Boot<-as.data.frame(Boot)
    colnames(Boot)<-c("Female_TEX114","Male_TEX1","Female_LPB87","Male_STM2")
    bon <-seq(from=debut, to=(end-wind), by=wind)
    
    for(i in bon){ 
      #  print(i)
      i0_A <- LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)),3]
      if( LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)),3] > i){ 
        i0_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))
        i0_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))-1 } else {
          i0_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))+1
          i0_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)) }
      T2_A <- data.frame(bp=c(LG11_A[i0_real_A_down,3],LG11_A[i0_real_A_up,3]), cm=(c(LG11_A[i0_real_A_down,4],LG11_A[i0_real_A_up,4])))
      T1_A <-data.frame(bp = i)
      i0_bon_A <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y 
      fin_est_A <- i + wind
      if( LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A)),3] > fin_est_A){ 
        fin_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))
        fin_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))-1 } else {
          fin_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))+1
          fin_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))}
      Table2_A <- data.frame(bp=c(LG11_A[fin_real_A_down,3],LG11_A[fin_real_A_up,3]), cm=(c(LG11_A[fin_real_A_down,4],LG11_A[fin_real_A_up,4])))
      Table1_A <-data.frame(bp = fin_est_A )
      if ( any(is.na((Table2_A)))) { final_bon_A <- na.omit(Table2_A)$cm } else {
        final_bon_A <- approx(Table2_A$bp, Table2_A$cm, xout = Table1_A$bp, method = "linear")$y }
      length_A <- final_bon_A - i0_bon_A
      Boot[i,1] <- length_A
      
      i0_B <- LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-i)),3]
      if( LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-i)),3] > i){ 
        i0_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))
        i0_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))-1 } else {
          i0_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))+1
          i0_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))   }
      T2_B <- data.frame(bp=c(LG11_B[i0_real_B_down,3],LG11_B[i0_real_B_up,3]), cm=(c(LG11_B[i0_real_B_down,4],LG11_B[i0_real_B_up,4])))
      T1_B <-data.frame(bp = i)
      i0_bon_B <- approx(T2_B$bp, T2_B$cm, xout = T1_B$bp, method = "linear")$y 
      fin_est_B <- i + wind
      if( LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B)),3] > fin_est_B){ 
        fin_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))
        fin_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))-1 } else {
          fin_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))+1
          fin_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))    }
      Table2_B <- data.frame(bp=c(LG11_B[fin_real_B_down,3],LG11_B[fin_real_B_up,3]), cm=(c(LG11_B[fin_real_B_down,4],LG11_B[fin_real_B_up,4])))
      Table1_B <-data.frame(bp = fin_est_B )
      if ( any(is.na((Table2_B)))) { final_bon_B <- na.omit(Table2_B)$cm } else {
        final_bon_B <- approx(Table2_B$bp, Table2_B$cm, xout = Table1_B$bp, method = "linear")$y }
      length_B <- final_bon_B - i0_bon_B
      Boot[i,2] <- length_B
      
      i0_C <- LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)),3]
      if( LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)),3] > i){ 
        i0_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))
        i0_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))-1 } else {
          i0_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))+1
          i0_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)) 
        }
      T2_C <- data.frame(bp=c(LG11_C[i0_real_C_down,3],LG11_C[i0_real_C_up,3]), cm=(c(LG11_C[i0_real_C_down,4],LG11_C[i0_real_C_up,4])))
      T1_C <-data.frame(bp = i)
      i0_bon_C <- approx(T2_C$bp, T2_C$cm, xout = T1_C$bp, method = "linear")$y 
      fin_est_C <- i + wind
      if( LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C)),3] > fin_est_C){ 
        fin_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))
        fin_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))-1 } else {
          fin_real_C_up <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))+1
          fin_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))
        }
      Table2_C <- data.frame(bp=c(LG11_C[fin_real_C_down,3],LG11_C[fin_real_C_up,3]), cm=(c(LG11_C[fin_real_C_down,4],LG11_C[fin_real_C_up,4])))
      Table1_C <-data.frame(bp = fin_est_C )
      if ( any(is.na((Table2_C)))) { final_bon_C <- na.omit(Table2_C)$cm } else {
        final_bon_C <- approx(Table2_C$bp, Table2_C$cm, xout = Table1_C$bp, method = "linear")$y }
      length_C <- final_bon_C - i0_bon_C
      Boot[i,3] <- length_C
      
      i0_D <- LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)),3]
      if( LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)),3] > i){ 
        i0_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))
        i0_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))-1 } else {
          i0_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))+1
          i0_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)) }
      T2_D <- data.frame(bp=c(LG11_D[i0_real_D_down,3],LG11_D[i0_real_D_up,3]), cm=(c(LG11_D[i0_real_D_down,4],LG11_D[i0_real_D_up,4])))
      T1_D <-data.frame(bp = i)
      i0_bon_D <- approx(T2_D$bp, T2_D$cm, xout = T1_D$bp, method = "linear")$y 
      fin_est_D <- i + wind
      if( LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D)),3] > fin_est_D){ 
        fin_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))
        fin_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))-1 } else {
          fin_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))+1
          fin_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))}
      Table2_D <- data.frame(bp=c(LG11_D[fin_real_D_down,3],LG11_D[fin_real_D_up,3]), cm=(c(LG11_D[fin_real_D_down,4],LG11_D[fin_real_D_up,4])))
      Table1_D <-data.frame(bp = fin_est_D )
      if ( any(is.na((Table2_D)))) { final_bon_D <- na.omit(Table2_D)$cm } else {
        final_bon_D <- approx(Table2_D$bp, Table2_D$cm, xout = Table1_D$bp, method = "linear")$y }
      length_D <- final_bon_D - i0_bon_D
      Boot[i,4] <- length_D
    }
    BootLG2 <- na.omit(Boot)
    BootLG2 <- cbind(BootLG2, zone= c(rep("zone1",5),rep("zone2",5),rep("zone3",5)))#,rep("zone4",5)))#,rep("zone5",4)))
  
  #test <- apply(Boot, 1, function(x) all(is.na(x)))
  #Boot2 <- Boot[ !test, ]
  finaldata[[k]] <- cbind(BootLG2, k)
}
  data_LG <- as.data.frame(do.call(rbind, finaldata))

#ANOVA
data_LG$window<-rownames(data_LG)
DATA_all <- melt(data_LG,id.var=c("window", "k", "zone"))
DATA_all$variable <- factor(DATA_all$variable, levels=c("Female_LPB87","Female_TEX114","Male_TEX1","Male_STM2"))
charac_zone <- as.character(unique(DATA_all$zone))
DATA_all$zone <- factor(DATA_all$zone, levels=charac_zone)
DATA_all$k <- as.character(DATA_all$k)
DATA_all$k <- factor(DATA_all$k, levels= unique(DATA_all$k))
DATA_all$window <- factor(DATA_all$window, levels= unique(DATA_all$window))
DATA_all$zone_not_nested <- paste0(DATA_all$zone,DATA_all$k)
colnames(DATA_all)<-c("window","LG","zone","carte","length", "zone_not_nested")

model = lmer(length ~ carte  + zone_not_nested + carte:zone_not_nested + (1|window), data=DATA_all)

# Create a QQ plot of residuals
library(ggpubr)
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro.test(residuals(model)) # NON NORMAL

# Data LG
data_LG$zone_not_nested <- paste0(data_LG$zone,data_LG$k)
# non parametric test :
A <- levels(factor(data_LG$zone_not_nested))
B <- levels(factor(data_LG$k))
library(ZIR)
result <- matrix(NA, ncol=13)
colnames(result) <- c("zone","LPB87vsTEX114","value_LT4","LPB87vsTEX1","value_LT","LPB87vsSTM2","value_LS","TEX114vsTEX1","value_TT","TEX114vsSTM2","value_T4S","TEX1vsSTM2","value_TS")
result <- as.data.frame(result)
for (i in 1: length(A)) {
  zone <- A[[i]]
data <- data_LG[data_LG$zone_not_nested==zone,]
l <- data$Female_LPB87
t4 <- data$Female_TEX114
t1 <- data$Male_TEX1
s <- data$Male_STM2
result[i,1]<- zone
result[i,2] <-  ziw(l, t4, perm = T)$p.value
result[i,3] <-  ziw(l, t4, perm = T)$statistics
result[i,4] <- ziw(l, t1, perm = TRUE)$p.value
result[i,5] <- ziw(l, t1, perm = TRUE)$statistics
result[i,6] <- ziw(l, s, perm = TRUE)$p.value
result[i,7] <- ziw(l, s, perm = TRUE)$statistics
result[i,8] <- ziw(t4, t1, perm = T)$p.value
result[i,9] <- ziw(t4, t1, perm = T)$statistics
result[i,10] <- ziw(t4, s, perm = TRUE)$p.value
result[i,11] <- ziw(t4, s, perm = TRUE)$statistics
result[i,12] <- ziw(t1, s, perm = TRUE)$p.value
result[i,13] <- ziw(t1, s, perm = TRUE)$statistics
}
statist<- result[,-c(2,4,6,8,10,12)]
statist <- melt(statist, id.vars = c("zone"))
adj <- melt(result[,c(1,2,4,6,8,10,12)],id.vars=c("zone"))
new_adj <- sapply(adj$value, function(x) { if (is.na(x)) {x <- "NA"}  else if (x == 0) {x <- 0.0001 } else {x <- x}}) 
new_adj <-as.numeric(new_adj)
adj <- cbind (adj, new_adj)
adj <-cbind(adj,p.adjust(adj$new_adj, method = "BH", n = length(adj$new_adj)))
adj <-cbind(adj,statist)
write.table(adj,"results_zone_scaled.txt")

result_LG <- matrix(NA, ncol=13)
colnames(result_LG) <- c("LG","LPB87vsTEX114","value_LT4","LPB87vsTEX1","value_LT","LPB87vsSTM2","value_LS","TEX114vsTEX1","value_TT","TEX114vsSTM2","value_T4S","TEX1vsSTM2","value_TS")
result_LG <- as.data.frame(result_LG)
#result with LG level
for (i in 1: length(B)) {
  LG <- B[[i]]
  data <- data_LG[data_LG$k==LG,]
  l <- data$Female_LPB87
  t4 <- data$Female_TEX114
  t1 <- data$Male_TEX1
  s <- data$Male_STM2
  result_LG[i,1]<- LG
  result_LG[i,2] <-  ziw(l, t4, perm = T)$p.value
  result_LG[i,3] <-  ziw(l, t4, perm = T)$statistics
  result_LG[i,4] <- ziw(l, t1, perm = TRUE)$p.value
  result_LG[i,5] <- ziw(l, t1, perm = TRUE)$statistics
  result_LG[i,6] <- ziw(l, s, perm = TRUE)$p.value
  result_LG[i,7] <- ziw(l, s, perm = TRUE)$statistics
  result_LG[i,8] <- ziw(t4, t1, perm = T)$p.value
  result_LG[i,9] <- ziw(t4, t1, perm = T)$statistics
  result_LG[i,10] <- ziw(t4, s, perm = TRUE)$p.value
  result_LG[i,11] <- ziw(t4, s, perm = TRUE)$statistics
  result_LG[i,12] <- ziw(t1, s, perm = TRUE)$p.value
  result_LG[i,13] <- ziw(t1, s, perm = TRUE)$statistics
}
statist_LG<- result_LG[,-c(2,4,6,8,10,12)]
statist_LG <- melt(statist_LG, id.vars = c("LG"))
adj_LG <- melt(result_LG[,c(1,2,4,6,8,10,12)],id.var=c("LG"))
new_adj_LG <- sapply(adj_LG$value, function(x) { if (is.na(x)) {x <- "NA"}  else if (x == 0) {x <- 0.0001 } else {x <- x}}) 
new_adj_LG <-as.numeric(new_adj_LG)
adj_LG <- cbind (adj_LG, new_adj_LG)
adj_LG <-cbind(adj_LG,p.adjust(adj_LG$new_adj, method = "BH", n = length(adj_LG$new_adj)))
adj_LG <-cbind(adj_LG,statist_LG)
write.table(adj_LG,"results_LG_scaled.txt")



##############################################################################################################################
#####  with the same data reduced and divided in 3 zones and 5 windows each, non parametric tests for the genetic length ####
##############################################################################################################################
## Determine the bin size fixed for all : 5 bin for each LG
library(reshape2)
phys = read.delim("reduced_physical_length.txt", stringsAsFactors = F)
bin_size <- phys$phys_length/15
phys <- cbind(phys, bin_size)

#the intial data is not reduced for the markers terminals markers but this is done in the script
All_data <- read.table("test.txt", header= T)

LG1<-subset(All_data, All_data$LG == "1" )
LG2<-subset(All_data, All_data$LG == "2" )
LG3<-subset(All_data, All_data$LG == "3" )
LG4<-subset(All_data, All_data$LG == "4" )
LG5<-subset(All_data, All_data$LG == "5" )
LG6<-subset(All_data, All_data$LG == "6" )
LG7<-subset(All_data, All_data$LG == "7" )
LG8<-subset(All_data, All_data$LG == "8" )
LG9<-subset(All_data, All_data$LG == "9" )
LG10<-subset(All_data, All_data$LG == "10" )
LG11<-subset(All_data, All_data$LG == "11" )
LG12<-subset(All_data, All_data$LG == "12" )

df_list <- list (LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10,LG11,LG12)
finaldata <- rep(list(matrix(NA)), length(df_list))
for (k in 1:length(df_list)) { 
  print(k)
  LG11_A<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_TEX114" )
  LG11_B<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_TEX1")
  LG11_C<-subset(df_list[[k]], df_list[[k]]$Parent== "Female_LPB87" )
  LG11_D<-subset(df_list[[k]], df_list[[k]]$Parent== "Male_STM2")
  end <- min(c(LG11_A[which.max(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.max(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.max(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.max(LG11_D$Estimated_physicalposition_bp.),3]))
  debut <- max(c(LG11_A[which.min(LG11_A$Estimated_physicalposition_bp.),3],LG11_B[which.min(LG11_B$Estimated_physicalposition_bp.),3],LG11_C[which.min(LG11_C$Estimated_physicalposition_bp.),3],LG11_D[which.min(LG11_D$Estimated_physicalposition_bp.),3]))
  wind <- phys[k,3]
  ##### LG sifnificantly lower than other LG
  Boot <- matrix(NA, nrow = 1, ncol = 4)
  Boot<-as.data.frame(Boot)
  colnames(Boot)<-c("Female_TEX114","Male_TEX1","Female_LPB87","Male_STM2")
  bon <-seq(from=debut, to=(end-wind), by=wind)
  
  for(i in bon){ 
    #  print(i)
    i0_A <- LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)),3]
    if( LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)),3] > i){ 
      i0_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))
      i0_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))-1 } else {
        i0_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-i))+1
        i0_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-i)) }
    T2_A <- data.frame(bp=c(LG11_A[i0_real_A_down,3],LG11_A[i0_real_A_up,3]), cm=(c(LG11_A[i0_real_A_down,4],LG11_A[i0_real_A_up,4])))
    T1_A <-data.frame(bp = i)
    i0_bon_A <- approx(T2_A$bp, T2_A$cm, xout = T1_A$bp, method = "linear")$y 
    fin_est_A <- i + wind
    if( LG11_A[which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A)),3] > fin_est_A){ 
      fin_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))
      fin_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))-1 } else {
        fin_real_A_up <-which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))+1
        fin_real_A_down <- which.min(abs(LG11_A$Estimated_physicalposition_bp.-fin_est_A))}
    Table2_A <- data.frame(bp=c(LG11_A[fin_real_A_down,3],LG11_A[fin_real_A_up,3]), cm=(c(LG11_A[fin_real_A_down,4],LG11_A[fin_real_A_up,4])))
    Table1_A <-data.frame(bp = fin_est_A )
    if ( any(is.na((Table2_A)))) { final_bon_A <- na.omit(Table2_A)$cm } else {
      final_bon_A <- approx(Table2_A$bp, Table2_A$cm, xout = Table1_A$bp, method = "linear")$y }
    length_A <- final_bon_A - i0_bon_A
    Boot[i,1] <- length_A
    
    i0_B <- LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-i)),3]
    if( LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-i)),3] > i){ 
      i0_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))
      i0_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))-1 } else {
        i0_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))+1
        i0_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-i))   }
    T2_B <- data.frame(bp=c(LG11_B[i0_real_B_down,3],LG11_B[i0_real_B_up,3]), cm=(c(LG11_B[i0_real_B_down,4],LG11_B[i0_real_B_up,4])))
    T1_B <-data.frame(bp = i)
    i0_bon_B <- approx(T2_B$bp, T2_B$cm, xout = T1_B$bp, method = "linear")$y 
    fin_est_B <- i + wind
    if( LG11_B[which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B)),3] > fin_est_B){ 
      fin_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))
      fin_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))-1 } else {
        fin_real_B_up <-which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))+1
        fin_real_B_down <- which.min(abs(LG11_B$Estimated_physicalposition_bp.-fin_est_B))    }
    Table2_B <- data.frame(bp=c(LG11_B[fin_real_B_down,3],LG11_B[fin_real_B_up,3]), cm=(c(LG11_B[fin_real_B_down,4],LG11_B[fin_real_B_up,4])))
    Table1_B <-data.frame(bp = fin_est_B )
    if ( any(is.na((Table2_B)))) { final_bon_B <- na.omit(Table2_B)$cm } else {
      final_bon_B <- approx(Table2_B$bp, Table2_B$cm, xout = Table1_B$bp, method = "linear")$y }
    length_B <- final_bon_B - i0_bon_B
    Boot[i,2] <- length_B
    
    i0_C <- LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)),3]
    if( LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)),3] > i){ 
      i0_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))
      i0_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))-1 } else {
        i0_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-i))+1
        i0_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-i)) 
      }
    T2_C <- data.frame(bp=c(LG11_C[i0_real_C_down,3],LG11_C[i0_real_C_up,3]), cm=(c(LG11_C[i0_real_C_down,4],LG11_C[i0_real_C_up,4])))
    T1_C <-data.frame(bp = i)
    i0_bon_C <- approx(T2_C$bp, T2_C$cm, xout = T1_C$bp, method = "linear")$y 
    fin_est_C <- i + wind
    if( LG11_C[which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C)),3] > fin_est_C){ 
      fin_real_C_up <-which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))
      fin_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))-1 } else {
        fin_real_C_up <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))+1
        fin_real_C_down <- which.min(abs(LG11_C$Estimated_physicalposition_bp.-fin_est_C))
      }
    Table2_C <- data.frame(bp=c(LG11_C[fin_real_C_down,3],LG11_C[fin_real_C_up,3]), cm=(c(LG11_C[fin_real_C_down,4],LG11_C[fin_real_C_up,4])))
    Table1_C <-data.frame(bp = fin_est_C )
    if ( any(is.na((Table2_C)))) { final_bon_C <- na.omit(Table2_C)$cm } else {
      final_bon_C <- approx(Table2_C$bp, Table2_C$cm, xout = Table1_C$bp, method = "linear")$y }
    length_C <- final_bon_C - i0_bon_C
    Boot[i,3] <- length_C
    
    i0_D <- LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)),3]
    if( LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)),3] > i){ 
      i0_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))
      i0_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))-1 } else {
        i0_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-i))+1
        i0_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-i)) }
    T2_D <- data.frame(bp=c(LG11_D[i0_real_D_down,3],LG11_D[i0_real_D_up,3]), cm=(c(LG11_D[i0_real_D_down,4],LG11_D[i0_real_D_up,4])))
    T1_D <-data.frame(bp = i)
    i0_bon_D <- approx(T2_D$bp, T2_D$cm, xout = T1_D$bp, method = "linear")$y 
    fin_est_D <- i + wind
    if( LG11_D[which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D)),3] > fin_est_D){ 
      fin_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))
      fin_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))-1 } else {
        fin_real_D_up <-which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))+1
        fin_real_D_down <- which.min(abs(LG11_D$Estimated_physicalposition_bp.-fin_est_D))}
    Table2_D <- data.frame(bp=c(LG11_D[fin_real_D_down,3],LG11_D[fin_real_D_up,3]), cm=(c(LG11_D[fin_real_D_down,4],LG11_D[fin_real_D_up,4])))
    Table1_D <-data.frame(bp = fin_est_D )
    if ( any(is.na((Table2_D)))) { final_bon_D <- na.omit(Table2_D)$cm } else {
      final_bon_D <- approx(Table2_D$bp, Table2_D$cm, xout = Table1_D$bp, method = "linear")$y }
    length_D <- final_bon_D - i0_bon_D
    Boot[i,4] <- length_D
  }
  BootLG2 <- na.omit(Boot)
  BootLG2 <- cbind(BootLG2, zone= c(rep("zone1",5),rep("zone2",5),rep("zone3",5)))#,rep("zone4",5)))#,rep("zone5",4)))
  
  #test <- apply(Boot, 1, function(x) all(is.na(x)))
  #Boot2 <- Boot[ !test, ]
  finaldata[[k]] <- cbind(BootLG2, k)
}
data_LG <- as.data.frame(do.call(rbind, finaldata))

#ANOVA
data_LG$window<-rownames(data_LG)
DATA_all <- melt(data_LG,id.var=c("window", "k", "zone"))
DATA_all$variable <- factor(DATA_all$variable, levels=c("Female_LPB87","Female_TEX114","Male_TEX1","Male_STM2"))
charac_zone <- as.character(unique(DATA_all$zone))
DATA_all$zone <- factor(DATA_all$zone, levels=charac_zone)
DATA_all$k <- as.character(DATA_all$k)
DATA_all$k <- factor(DATA_all$k, levels= unique(DATA_all$k))
DATA_all$window <- factor(DATA_all$window, levels= unique(DATA_all$window))
DATA_all$zone_not_nested <- paste0(DATA_all$zone,DATA_all$k)
colnames(DATA_all)<-c("window","LG","zone","carte","length", "zone_not_nested")

model = lmer(length ~ carte  + zone_not_nested + carte:zone_not_nested + (1|window), data=DATA_all)
model2 = lmer(length ~ carte  + LG + carte:LG + (1|zone_not_nested), data=DATA_all)
# Create a QQ plot of residuals
library(ggpubr)
ggqqplot(residuals(model))
ggqqplot(residuals(model2))
# Compute Shapiro-Wilk test of normality
shapiro.test(residuals(model))
shapiro.test(residuals(model2)) # NON NORMAL
anova(model2) 

# Data LG
data_LG$zone_not_nested <- paste0(data_LG$zone,data_LG$k)
# non parametric test :
A <- levels(factor(data_LG$zone_not_nested))
B <- levels(factor(data_LG$k))
library(ZIR)
result <- matrix(NA, ncol=13)
colnames(result) <- c("zone","LPB87vsTEX114","value_LT4","LPB87vsTEX1","value_LT","LPB87vsSTM2","value_LS","TEX114vsTEX1","value_TT","TEX114vsSTM2","value_T4S","TEX1vsSTM2","value_TS")
result <- as.data.frame(result)

#result with zones
for (i in 1: length(A)) {
  zone <- A[[i]]
  data <- data_LG[data_LG$zone_not_nested==zone,]
  l <- data$Female_LPB87
  t4 <- data$Female_TEX114
  t1 <- data$Male_TEX1
  s <- data$Male_STM2
  result[i,1]<- zone
  result[i,2] <-  ziw(l, t4, perm = T)$p.value
  result[i,3] <-  ziw(l, t4, perm = T)$statistics
  result[i,4] <- ziw(l, t1, perm = TRUE)$p.value
  result[i,5] <- ziw(l, t1, perm = TRUE)$statistics
  result[i,6] <- ziw(l, s, perm = TRUE)$p.value
  result[i,7] <- ziw(l, s, perm = TRUE)$statistics
  result[i,8] <- ziw(t4, t1, perm = T)$p.value
  result[i,9] <- ziw(t4, t1, perm = T)$statistics
  result[i,10] <- ziw(t4, s, perm = TRUE)$p.value
  result[i,11] <- ziw(t4, s, perm = TRUE)$statistics
  result[i,12] <- ziw(t1, s, perm = TRUE)$p.value
  result[i,13] <- ziw(t1, s, perm = TRUE)$statistics
}
statist<- result[,-c(2,4,6,8,10,12)]
statist <- melt(statist, id.vars = c("zone"))
adj <- melt(result[,c(1,2,4,6,8,10,12)],id.vars=c("zone"))
new_adj <- sapply(adj$value, function(x) { if (is.na(x)) {x <- "NA"}  else if (x == 0) {x <- 0.0001 } else {x <- x}}) 
new_adj <-as.numeric(new_adj)
adj <- cbind (adj, new_adj)
adj <-cbind(adj,p.adjust(adj$new_adj, method = "BH", n = length(adj$new_adj)))
adj <-cbind(adj,statist)

write.table(adj,"results_window_noscaled.txt")


result_LG <- matrix(NA, ncol=13)
colnames(result_LG) <- c("LG","LPB87vsTEX114","value_LT4","LPB87vsTEX1","value_LT","LPB87vsSTM2","value_LS","TEX114vsTEX1","value_TT","TEX114vsSTM2","value_T4S","TEX1vsSTM2","value_TS")
result_LG <- as.data.frame(result_LG)
#result with LG level
for (i in 1: length(B)) {
 LG <- B[[i]]
  data <- data_LG[data_LG$k==LG,]
  l <- data$Female_LPB87
  t4 <- data$Female_TEX114
  t1 <- data$Male_TEX1
  s <- data$Male_STM2
  result_LG[i,1]<- LG
  result_LG[i,2] <-  ziw(l, t4, perm = T)$p.value
  result_LG[i,3] <-  ziw(l, t4, perm = T)$statistics
  result_LG[i,4] <- ziw(l, t1, perm = TRUE)$p.value
  result_LG[i,5] <- ziw(l, t1, perm = TRUE)$statistics
  result_LG[i,6] <- ziw(l, s, perm = TRUE)$p.value
  result_LG[i,7] <- ziw(l, s, perm = TRUE)$statistics
  result_LG[i,8] <- ziw(t4, t1, perm = T)$p.value
  result_LG[i,9] <- ziw(t4, t1, perm = T)$statistics
  result_LG[i,10] <- ziw(t4, s, perm = TRUE)$p.value
  result_LG[i,11] <- ziw(t4, s, perm = TRUE)$statistics
  result_LG[i,12] <- ziw(t1, s, perm = TRUE)$p.value
  result_LG[i,13] <- ziw(t1, s, perm = TRUE)$statistics
}
statist_LG<- result_LG[,-c(2,4,6,8,10,12)]
statist_LG <- melt(statist_LG, id.vars = c("LG"))
adj_LG <- melt(result_LG[,c(1,2,4,6,8,10,12)],id.var=c("LG"))
new_adj_LG <- sapply(adj_LG$value, function(x) { if (is.na(x)) {x <- "NA"}  else if (x == 0) {x <- 0.0001 } else {x <- x}}) 
new_adj_LG <-as.numeric(new_adj_LG)
adj_LG <- cbind (adj_LG, new_adj_LG)
adj_LG <-cbind(adj_LG,p.adjust(adj_LG$new_adj, method = "BH", n = length(adj_LG$new_adj)))
adj_LG <-cbind(adj_LG,statist_LG)
write.table(adj_LG,"results_LG_noscaled.txt")




