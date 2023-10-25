# Produces a comprehensive list of cazyme families by removing all tertiary data from the columns and getting an intersect from
# both predictors


get_familys = function(df){
  
  # df = df %>% 
  #   filter(DIAMOND != "-" | HMMER != "-" | dbCAN_sub != "-") # removing rows where no family was identified by any tools
  
  df = df %>% 
    filter(X.ofTools >= 2)
  
  df$family = NA # adding a family column
  
  df$contig = sapply(df$Gene.ID, function(x){
    paste0(strsplit(as.character(x), "_", fixed = TRUE)[[1]][-4], collapse = "_")
  })
  
  ## Take left most where multiple exist ##
  df$DIAMOND = sapply(df$DIAMOND, function(x){
    strsplit(as.character(x), "+", fixed = TRUE)[[1]][1] # Splitting the column at "+" and taking the first cazyme family (has the most proteins matched to this family)
  })
  
  df$HMMER = sapply(df$HMMER, function(x){
    strsplit(as.character(x), "+", fixed = TRUE)[[1]][1]
  })
  
  df$dbCAN_sub = sapply(df$dbCAN_sub, function(x){
    strsplit(as.character(x), "+", fixed = TRUE)[[1]][1]
  })
  
  ## Remove tertiary info ##
  
  df$DIAMOND = sapply(df$DIAMOND, function(x){
    strsplit(as.character(x), "_", fixed = TRUE)[[1]][1] # Splitting the column at "_" to remove the protein number and extra stuff dbSUB has (not sure what it is)
  })
  
  df$HMMER = sapply(df$HMMER, function(x){
    strsplit(as.character(x), "_", fixed = TRUE)[[1]][1]
  })
  
  df$dbCAN_sub = sapply(df$dbCAN_sub, function(x){
    strsplit(as.character(x), "_", fixed = TRUE)[[1]][1]
  })

  df$DIAMOND = sapply(df$DIAMOND, function(x){
    strsplit(as.character(x), "(", fixed = TRUE)[[1]][1] # Splitting the column at "(+)" not sure what is in the brackets though
  })
  
  df$HMMER = sapply(df$HMMER, function(x){
    strsplit(as.character(x), "(", fixed = TRUE)[[1]][1]
  })
  
  df$dbCAN_sub = sapply(df$dbCAN_sub, function(x){
    strsplit(as.character(x), "(", fixed = TRUE)[[1]][1]
  })
  
  ## Removing "-" and replacing with NA ##
  df = df %>%
    mutate(HMMER = na_if(HMMER, "-"),
           DIAMOND = na_if(DIAMOND, "-"),
           dbCAN_sub = na_if(dbCAN_sub, "-"))
  
  for(r in 1:nrow(df)){
    votes = c(df[r,]$HMMER, df[r,]$DIAMOND, df[r,]$dbCAN_sub)
    votes[! votes %in% c("GH0", "GT0", "PL0", "CE0", "AA0", "CBM0")] # ensuring we remove where it has been assigned to unclassified from the votes
    df[r,]$family = names(sort(table(votes),decreasing=TRUE)[1])
  }
  
  # df_coalescense = df %>% 
  #   # taking the coalenscense of the two columns with DIAMOND being the "master" column
  #   mutate(DIAMOND = na_if(DIAMOND, "-")) %>%
  #   mutate(family = coalesce(DIAMOND, HMMER))
  
  return(df)
  }
  

get_substrates = function(df){
  
  colnames(df) = c("dbCAN subfam",	"Subfam Composition",	"Subfam EC",	"Substrate",	"Profile Length",	"Gene ID",
                   "Gene Length",	"E Value",	"Profile Start", "Profile End",	"Gene Start",	"Gene End",	"Coverage", "Coverage2")
  df = df %>%
    filter(Substrate != "-")
  
  df$Substrate = sapply(df$Substrate, function(x){
    strsplit(as.character(x), ",", fixed = TRUE)[[1]][1]
  })
  
  df$`Gene ID` = sapply(df$`Gene ID`, function(x){
    paste0(strsplit(as.character(x), "_", fixed = TRUE)[[1]][-4], collapse = "_")
  })
  
  return(df)
}
