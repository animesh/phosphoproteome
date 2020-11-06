#' Extract the phosphorylated proteins and the corresponding phosphorylated sites.
#'
#' Extract the phosphorylated proteins and the corresponding phosphorylated sites from the export .xlsx data created by Proteome Discovery 2.4
#' Note: only the first protein and the first position will be be remained
#' @param x A data.frame
#' @return A data.frame in whcih the first column is the unique phosphorylated protein accession and the second column is the sum of the phosphorylated sites
#' @import stringr
#' @export
PhosphorylatedProteinInformation<-function(x){

  #############################################################################################
  ## x is one sheet of .xlsx file of phosphoproteome exported from Proteome Discovery 2.4.   ##
  ## 目的：1. 提取`Modifications in Proteins`信息，并将其分成phosphoProtein和Site两列        ##
  ##       2.  a. 同一个肽段匹配的protein>2的仅仅考虑第一个蛋白；                            ##
  ##           b. 同一个肽段，在同一个protein中有多处匹配，仅考虑第一个匹配的区段；          ##
  ## 输出：第一列为非冗余的protein name; 第二列为汇总的磷酸化位点信息，包括磷酸化的AA及打分。##
  #############################################################################################

  # 1. 分成phosphoProtein和Site信息分成两列；
  modification_1<-str_split(x$`Modifications in Proteins`," ",n=2,simplify = T) %>% data.frame(.,stringsAsFactors = F)

  # 2. 多个蛋白和同一个蛋白多个位置匹配的情况，仅考虑第一个protein和滴哟个位点，此外，没有打分的位点也被删除部分，如 S/Y/T；
  modification_1$X2<-str_split(modification_1$X2,"]; ",n=2,simplify = T) %>% data.frame(.,stringsAsFactors = F) %>% .[,1]

  # 3. 将数据整理为data.frame,第一列为非冗余的protein，第二列为磷酸化位点信息
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

  ############################################################################################
  ## 输入x为: 非facter类型的Data.frame，第一列为protein accession,第二列为磷酸化位点信息；  ##
  ## 输出: protein相同的行合并为一行。                                                      ##
  ############################################################################################

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

  ####################################################################################################################
  ## x       is the result of PhosphorylatedProteinInformation,namely a data.frame containing one column of         ##
  ##         unique protein and all of the phosphorylated information.                                              ##
  ## export  a list which consist of data.frame, the name of each data.frame is the unique protein accession,       ##
  ##         the data.frame contains the separated site and phosphoScore and HighConfident                          ##
  ####################################################################################################################

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

  ###################################################################################################
  ## x       the result of sumSites,namely a list of all phosphorylated proteins。                 ##
  ## export   a data. frame contains Total phosphoSites,                                           ##
  ##                                 The most phosphorylated protein,                              ##
  ##                                 The number of phosphoSites of most phosphorylated protein.    ##
  ####################################################################################################

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
  ###############################################################################################
  ### x, y   vector，x is character, y is numeric, the position of AA, which will be tolower   ##
  ##  for example: x= c("ABEHUI"), y=c(1,4,5)                                                  ##
  ###############################################################################################
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
  ###############################################################################################################
  ## x       a list, the result of sumHighConfidentSites                                                                    ##
  ## export  a fasta file of identified phosphorylated proteins, which was labeled with phosphorylated sites   ##
  ###############################################################################################################
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

