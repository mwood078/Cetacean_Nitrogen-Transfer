## Run Sensitivity Analysis ##

## Initialization ----
rm(list=(ls()))

load_libraries<-function() {
  library(PBSmodelling)
  library(snowfall)
  library(parallel)
  library(snow)
  library(foreach)
  library(doSNOW)
  library(PBSadmb)
  library(tidyr)
  library(xlsx)
  library(janitor)
  library(dplyr)
  
}
load_libraries()


## Set global working directory
workdir <- "D:/Cetacean Nutrient Transfer"
setwd(workdir)


## Import datasets with clean_names() function
abundance <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Abundance"))
history <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Life History"))
fish_prox <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Fish Proximate"))
ceph_prox <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Cephalopod Proximate"))
crust_prox <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Crustacean Proximate"))
dive_table <- clean_names(read.xlsx(file=paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheetName = "Dive Table"))

n_boot <-100  #number of runs

#* Dataframe of input parameters
inputs <- data.frame(Simulation = rep(0,n_boot),
                     Sperm_Whale_n = rep(0,n_boot),
                     Rices_Whale_n = rep(0,n_boot),
                     CBW_n = rep(0,n_boot),
                     BBW_n = rep(0,n_boot),
                     GBW_n = rep(0,n_boot),
                     CBDO_n = rep(0,n_boot),
                     PSD_n = rep(0,n_boot),
                     STD_n = rep(0,n_boot),
                     SPD_n = rep(0,n_boot),
                     CD_n = rep(0,n_boot),
                     FD_n = rep(0,n_boot),
                     KW_n = rep(0,n_boot),
                     FKW_n = rep(0,n_boot),
                     PKW_n = rep(0,n_boot),
                     DSW_n = rep(0,n_boot),
                     PSW_n = rep(0,n_boot),
                     MHW_n = rep(0,n_boot),
                     RD_n = rep(0,n_boot),
                     SFPW_n = rep(0,n_boot),
                     fish_prot_prop = rep(0,n_boot),
                     ceph_prot_prop = rep(0,n_boot),
                     crust_prot_prop = rep(0,n_boot),
                     n_pods = rep(0,n_boot))


n_ts <- (24*60)/1 #Number of minutes in a day divided by the length of the actual time step (e.g., X/1 = 1 minute time step);  1 minute makes sense because some surface intervals are less than 5 minutes


#* Protein concentration for fishes
pro_prop_fish_mean_init <- mean(c(fish_prox$mean_protein_concentration_wet_weight))/100 #*Stickney and Torres (1989) %Wet Weight
pro_prop_fish_sd_init <- sd(c(fish_prox$low_protein_concentration_wet_weight))/100 #* Stickney and Torres (1989)

#* Protein concentration for cephalopods
pro_prop_ceph_mean_init <- mean(c(ceph_prox$mean_protein_concentration_wet_weight))/100 #*Sinclair et al. 2015 %Wet Weight
pro_prop_ceph_sd_init <- sd(c(ceph_prox$mean_protein_concentration_wet_weight))/100 #* Sinclair et al. 2015

#* Protein concentratoin for crustaceans
pro_prop_crust_mean_init <- mean(c(na.omit(crust_prox$mean_protein_concentration_wet_weight)))/100 #*Donnelly et al. 1993 %Wet Weight
pro_prop_crust_sd_init <- sd(c(na.omit(crust_prox$mean_protein_concentration_wet_weight)))/100 #* Donnelly et al. 1993


resdir <- paste(workdir,"/Results/Sensitivity/",sep="") #* Results output dataframe

# run the simulations section

#set up parallel 
no_cores <- detectCores()  #determine number of cores on computer
cl<-makeCluster(no_cores)  #setup the size of your parallel cluster based on number of cores
registerDoSNOW(cl)         #register cluster to get it ready for parallel processing

#set up text progress bar....this is personal code just to have a progress bar letting your know % completion
pb <- txtProgressBar(max = n_boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



# following is the main code for running in parallel
#ls is your parallel processing object, essentially the same as a for loop...need to include any packages that are required for use in the loop
#the code after setting up the parallel processing object is the loop you want to run on each core, this code will split the number of runs evenly across cores
# here the number of runs is set by nsim (so run the code following ls nsim times, and divide nsim/ncores)
#if need any info saved from each run, make sure that code explicitly saves that information to hardrive (e.g., as csv, etc.) otherwise values will not be saved into memory

ls=foreach(sim=1:n_boot,.options.snow = opts,.combine='rbind',.packages =c('PBSmodelling')) %dopar% {
  
  ## Base Model (No variation in all tested parameters) ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    #groups$Meso_fish_Abun[group] <- 0.370 ## Average
    #groups$Meso_ceph_Abun[group] <- 0.400 ## Average
    #groups$Meso_crust_Abun[group] <- 0.526 ## Average
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Base/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Base/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Base/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Base/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Base/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-based dataframes ##
  setwd(paste(resdir,"Base/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## Abundance Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  for (spec in 1:n_distinct(abundance$species)){
    while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
      general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
    }
    inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  }
  
  #general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    #groups$Meso_fish_Abun[group] <- 0.370 ## Average
    #groups$Meso_ceph_Abun[group] <- 0.400 ## Average
    #groups$Meso_crust_Abun[group] <- 0.526 ## Average
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Abundance/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Abundance/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Abundance/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Abundance/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Abundance/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-Abundanced dataframes ##
  setwd(paste(resdir,"Abundance/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## N proportion Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  #inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  #inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  #inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    #groups$Meso_fish_Abun[group] <- 0.370 ## Average
    #groups$Meso_ceph_Abun[group] <- 0.400 ## Average
    #groups$Meso_crust_Abun[group] <- 0.526 ## Average
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"N proportion/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"N proportion/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"N proportion/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"N proportion/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"N proportion/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-N proportiond dataframes ##
  setwd(paste(resdir,"N proportion/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## Protein N Concentration Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- rnorm(1,0.17,0.17*0.2) #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    #groups$Meso_fish_Abun[group] <- 0.370 ## Average
    #groups$Meso_ceph_Abun[group] <- 0.400 ## Average
    #groups$Meso_crust_Abun[group] <- 0.526 ## Average
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Protein N Concentration/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Protein N Concentration/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Protein N Concentration/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Protein N Concentration/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Protein N Concentration/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-Protein N Concentrationd dataframes ##
  setwd(paste(resdir,"Protein N Concentration/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## N Excreted Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- rnorm(1,0.8,0.8*0.2) #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    #groups$Meso_fish_Abun[group] <- 0.370 ## Average
    #groups$Meso_ceph_Abun[group] <- 0.400 ## Average
    #groups$Meso_crust_Abun[group] <- 0.526 ## Average
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"N Excreted/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"N Excreted/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"N Excreted/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"N Excreted/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"N Excreted/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-N Excretedd dataframes ##
  setwd(paste(resdir,"N Excreted/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## Meso Abundance Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    
    groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    #groups$Meso_fish_Abun[group] <- 0.369862
    #groups$Meso_ceph_Abun[group] <- 0.3999993
    #groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Meso Abundance/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Meso Abundance/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Meso Abundance/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Meso Abundance/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Meso Abundance/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-Meso Abundanced dataframes ##
  setwd(paste(resdir,"Meso Abundance/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## Intervals Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      #deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      #deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
        shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      }
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
      #  shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      #}
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          #cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Intervals/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Intervals/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Intervals/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Intervals/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Intervals/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-Intervalsd dataframes ##
  setwd(paste(resdir,"Intervals/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
  ## Consumption Rate Variation ----
  
  
  library(xlsx)
  library(janitor)
  library(dplyr)
  
  
  prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
  prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this
  
  
  
  cat("\n**** Iteration ",sim," ****\n",date(),"\n")
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Estimates a random value between minimum and maximum abundance
  
  #for (spec in 1:n_distinct(abundance$species)){
  #  while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
  #    general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
  #  }
  #  inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  #}
  
  general$Abundance_n <- abundance$abundance
  
  
  ## Sample pod sizes and distributions
  cat("Generating Groups\n")
  
  #* Dataframe containing information for each group
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change group size to match reserve to get to population size for final
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  inputs$n_pods[sim] <- n_group
  groups <- groups[1:(n_group-1),] #Remove Excess Rows
  
  ## Resample Mesopelagic N proportion around mean 
  #inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init)
  #inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init)
  #inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init)
  
  inputs$fish_prot_prop[sim] <- pro_prop_fish_mean_init
  inputs$ceph_prot_prop[sim] <- pro_prop_ceph_mean_init
  inputs$crust_prot_prop[sim] <- pro_prop_crust_mean_init
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  cat("Placing groups at Lat and Longs\n")
  groups$Latitude <- c(24)
  groups$Longitude <- c(-90)
  groups$Bottom_Depth <- c(5000)
  
  
  #### Model Run ---
  cat("\nRunning Model\n")
  pb <- txtProgressBar(min=1,max=dim(groups)[1],style=3)
  
  
  
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group)
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) ##
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    ## Counters for time spent in a dive/rest stage
    
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proprtions of mesopelagics at night within the epipelagic zone
    
    #groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    #groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    #groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    groups$Meso_fish_Abun[group] <- 0.369862
    groups$Meso_ceph_Abun[group] <- 0.3999993
    groups$Meso_crust_Abun[group] <- 0.526189
    
    for (ts in 1:n_ts){
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval
      #deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      #deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      deep_dive_int <- round(dive_table$deep_dive_duration_min[species_num],0)
      deep_surface_int <- round(dive_table$deep_surface_interval_min[species_num],0)
      
      
      #if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
      #  shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
      #  shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      #}
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(dive_table$shallow_dive_duration_min[species_num],0)
        shallow_surface_int <- round(dive_table$shallow_surface_interval_min[species_num],0)
      }
      
      
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        ## Initiate deep dive sequence ##
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
        }
        
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group
        
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)
        #cons <- history$mean_consumption_rate_kg_day[spec]
        
        # Calculate Ration (kg/day)
        rat <- (cons*groups$Individuals_n[group])/dive_table$deep_n_dives_per_day[species_num]
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100
        
        ## MESOPELAGIC Nitrogen Consumed = Feces + Urine + Storage
        # Incoporate fish, cephalopod, and crustacean contributions
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
          ## Daytime feeding
          
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
          ## Nighttime feeding below 200m
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        } else {
          ## Daytime feeding in epipelagic zone
          #* Only a proportion of the diet is mesopelagic
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){
        ## Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){
        ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence ##
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 1% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          cons <- rnorm(1,history$mean_consumption_rate_kg_day[spec],history$mean_consumption_rate_kg_day[spec]*0.2)         
          #cons <- history$mean_consumption_rate_kg_day[spec]  
          
          # Literature Consumption Rate (kg/day)
          rat <- (cons*groups$Individuals_n[group])/dive_table$shallow_n_dives_per_day[species_num]
          
          
          ## Percent bodyweight consumed
          bwgt<- (rat/groups$Biomass_kg[group])*100
          
          ## Nitrogen Consumed = Feces + Urine + Storage
          #* Same as deep dive
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
            ## Daytime feeding
            
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else if ((ts/n_ts > (420/1440) || ts/n_ts < (1140/1440)) && dives$Depth[ts] > 200){
            ## Nighttime feeding below 200m
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
            
          } else {
            ## Daytime feeding in epipelagic zone
            #* Only a proportion of the diet is mesopelagic
            con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n
          }   
          
          #* Assumed 80% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"Consumption Rate/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"Consumption Rate/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"Consumption Rate/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"Consumption Rate/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"Consumption Rate/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-Consumption Rated dataframes ##
  setwd(paste(resdir,"Consumption Rate/Sim ",sim,sep=""))
  
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
}


# end parallel job
close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 


