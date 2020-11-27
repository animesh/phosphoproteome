#' Extract the phosphorylated proteins and the corresponding phosphorylated sites.
#' Extract the phosphorylated proteins and the corresponding phosphorylated sites from the export .xlsx data created by Proteome Discovery 2.4
#' Note: only the first protein and the first position will be be remained
#' @param x A data.frame
#' @return A data.frame in whcih the first column is the unique phosphorylated protein accession and the second column is the sum of the phosphorylated sites
#' @import stringr
#' @export

PhosphorylatedProteinInformation<-function(x){

  ## x is one sheet of .xlsx file of phosphoproteome exported from Proteome Discovery 2.4.   ##
  ## Purpose: 1. Extract `Modifications in Proteins` information and divide it into two columns, phosphoProtein and Site ##
  ## 2. a. Only the first protein is considered if the protein>2 matches the same peptide; ##
  ## b. The same peptide segment has multiple matches in the same protein, and only the first matching segment is considered; ##
  ## Output: The first column is the non-redundant protein name; the second column is the summary phosphorylation site information, including phosphorylated AA and score. ##  #############################################################################################

  # 1. Divide the phosphoProtein and Site information into two columns;
  modification_1<-str_split(x$`Modifications in Proteins`," ",n=2,simplify = T) %>% data.frame(.,stringsAsFactors = F)

  # 2. When multiple proteins are matched with multiple positions of the same protein, only the first protein and one site are considered. In addition, the sites that are not scored are also deleted, such as S/Y/T;
  modification_1$X2<-str_split(modification_1$X2,"]; ",n=2,simplify = T) %>% data.frame(.,stringsAsFactors = F) %>% .[,1]

  # 3. Organize the data into data.frame, the first column is non-redundant protein, and the second column is phosphorylation site information
  stringPro<-modification_1[,1] %>% unique %>% na.omit
  ProDataFrame<-data.frame()
  for (j in 1:length(stringPro)){
    ProDataFrame[j,1]<-stringPro[j]
    Matrix<-modification_1[grep(paste(stringPro[j]),modification_1[,1]),]
    ProDataFrame[j,2]<-nrowToSingleString(Matrix)
  }
  colnames(ProDataFrame)<-c("protein","PosphoSites")
  return(ProDataFrame)
}

#' Merge some character into one character
#'
#' Merge the information of different phosphorylated sites in one protein into one row in a data.frame
#' @param x A data.frame that contains the same one protein accession and different phosphorylated sites
#' @return a vector contains one character in which all the phosphorylated site information were recorded
#' @export
#'

nrowToSingleString<-function(x){

  ## Input x is: non-facter type Data.frame, the first column is protein accession, the second column is phosphorylation site information; ##
  ## Output: The lines with the same protein are merged into one line.               ##

  y=c()
  if (nrow(x)==1) {y=paste(x[1,2])}
  if (nrow(x)!=1) {for (e in 1:nrow(x)){y=paste(y,x[e,2])}}
  return(y)
}


#' Extract the high confident phosphorylated sites
#'
#' High confident phosphorylated sites: score > 75
#' @param x the result of PhosphorylatedProteinInformation,namely a data.frame containing one column of unique protein and all of the phosphorylated information.
#' @return a list which consist of data.frame, the name of each data.frame is the unique protein accession, the data.frame contains the separated site and phosphoScore and HighConfident
#' @export
#'

sumHighConfidentSites<-function(x){

  ## x       is the result of PhosphorylatedProteinInformation,namely a data.frame containing one column of         ##
  ##         unique protein and all of the phosphorylated information.                                              ##
  ## export  a list which consist of data.frame, the name of each data.frame is the unique protein accession,       ##
  ##         the data.frame contains the separated site and phosphoScore and HighConfident                          ##

  y<-x[,2]
  proList<-str_extract_all(y,"[STY]{1}[0-9]+\\([0-9]+\\.?[0-9]?\\)")
  siteList<-list()
  siteListHighScore<-list()
  for (i in 1:length(proList)){
    siteList[[i]]<-str_split(proList[[i]],"\\(|\\)",n=3,simplify = T)%>% data.frame(.,stringsAsFactors = F)
    for (j in 1:nrow(siteList[[i]])){
      if (isTRUE((as.numeric(siteList[[i]][j,2])<75))){
        siteList[[i]][j,3]<-NA
      }
      if (isTRUE(as.numeric(siteList[[i]][j,2])>=75)){
        siteList[[i]][j,3]<-"HighScore"
      }
    }

    # 1. Delete the sites with PhosphoSites Score < 75
    siteListHighScore[[i]]<-na.omit(siteList[[i]])

    # 2. Change the typeoff of Score: form character to numeric
    siteListHighScore[[i]]$X2<-as.numeric(siteListHighScore[[i]]$X2)

    # 3. Sort the Score in descending order to delete the lower score of repeated site
    siteListHighScore[[i]]<-siteListHighScore[[i]][order(siteListHighScore[[i]][,2],decreasing=T),]
    siteListHighScore[[i]]<-siteListHighScore[[i]][!duplicated(siteListHighScore[[i]]$X1),]
    colnames(siteListHighScore[[i]])<-c("Site","PhosphoScore","HighConfident")
  }
  names(siteListHighScore)<-x[,1]
  return(siteListHighScore)
}


#' Extract the high confident phosphorylated sites
#'
#' The summary of Total phosphoSites and the most phosphorylated protein
#' @param x the result of sumHighConfidentSites,namely a list of all phosphorylated proteins。
#' @return  a data. frame contains Total phosphoSites, The most phosphorylated protein, The number of phosphoSites of most phosphorylated protein.
#' @export
#'

SummaryTotalPhosphosites<-function(x){

  ## x       the result of sumSites,namely a list of all phosphorylated proteins。                 ##
  ## export   a data. frame contains Total phosphoSites,                                           ##
  ##                                 The most phosphorylated protein,                              ##
  ##                                 The number of phosphoSites of most phosphorylated protein.    ##

  count<-vector()
  for (i in 1: length(x)){
    count[i]<-nrow(x[[i]])
  }
  for (j in 1:length(count)){
    if (isTRUE(count[j]==max(count))){protein<-names(x)[j]}
  }
  sumcount<-data.frame('Total phosphoSites'= sum(count),'The most phosphorylated protein'= protein,'The number of phosphoSites of most phosphorylated protein'=max(count))
  return(sumcount)
}


#' Batch tolower specific postion letters
#'
#' To tolower specific postion letters in a vector that contains one character
#' @param x a vector that contains one character
#' @param y a numeric vector
#' @return  a vector contains one character that will be labeled with lowercases in the specific position
#' @export

tolowerSpecificSite<-function(x,y){
  ### x, y   vector，x is character, y is numeric, the position of AA, which will be tolower   ##
  ##  for example: x= c("ABEHUI"), y=c(1,4,5)                                                  ##
  if (isTRUE(is.na(x))){
    newAA<-NA
  }
  if (isFALSE(is.na(x))){
    matrixAA<-str_split(x,"",simplify=T) %>% t %>% data.frame(.,stringsAsFactors = F)
    for (i in y){
      matrixAA[i,]<-tolower(matrixAA[i,])
    }
    newAA<-paste(as.vector(matrixAA$.),collapse="")
  }
  return(newAA)
}

#' To label the phosphorylated amino acid in protein fasta file
#'
#' To tolower the phosphorylated amino acid in protein fasta file
#' @param x  a list, the result of sumHighConfidentSites
#' @param fasta a fasta dataset download form Uniprot
#' @return  a fasta file of identified phosphorylated proteins, which was labeled with phosphorylated sites
#' @export
#'

tolowerPhosphoSitesinFasta<-function(x,fasta){
  ## x       a list, the result of sumHighConfidentSites                                                                    ##
  ## export  a fasta file of identified phosphorylated proteins, which was labeled with phosphorylated sites   ##
  position<-list()
  fastaPro<-list()
  newfastaPro<-list()
  for (i in 1:length(x)){
    if (isTRUE(nrow(x[[i]])!=0)){
      position[[i]]<-str_split(x[[i]][,1],"[STY]",simplify = T)[,2] %>% as.numeric
      if (isTRUE(paste(names(x)[i])%in%str_split(names(fasta),'\\|',n=3,simplify = T)[,2])){
        fastaPro[[i]]<-fasta[[grep(paste(names(x)[i]),names(fasta))]][1]
      }
      if (isFALSE(paste(names(x)[i])%in%str_split(names(fasta),'\\|',n=3,simplify = T)[,2])){
        fastaPro[[i]]<-NA
      }
    }
    if (isTRUE(nrow(x[[i]])==0)){
      fastaPro[[i]]<-NA
      position[[i]]<-NA
    }
    newfastaPro[[i]]<-tolowerSpecificSite(fastaPro[[i]],position[[i]])
  }
  names(newfastaPro)<-names(x)
  return(newfastaPro)
}
