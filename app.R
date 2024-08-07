library(shiny)
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggpubr)
library(plotly)
library(processx)
library(DT)
library(shinyBS)
library(data.table)
library(dplyr)
library(reshape2)
library(assertthat)

#if (interactive()) {
ui <- fluidPage(
  titlePanel("Protein Annotator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose data matrix and sample file", accept = c(".txt",".tsv",".csv"),multiple = F),
      helpText("upload the data file in the specified format"),
      uiOutput("selectfile"),
      br(),
      uiOutput("Species"),
      br(),
      uiOutput("sampleGroup"),
      br(),
      uiOutput("vx"),
      br(),
      uiOutput("div"),
      br(),
      uiOutput("imput")
    ),
    mainPanel(
      tableOutput("contents"),
      uiOutput("tb"),
      #tabsetPanel(id = "tabset",tabPanel("Plot3", plotlyOutput("p3", height = "800px"))),
      #tags$p("This app takes the protein intensities from the two group of samples and does differential analysis"),
      #tags$p("The sample data file and sample data key file should be in a spcific format(given below) all samples in the two files should match"),
      #br(),
      #br(),
      #tags$strong("sample data file"),
      #br(),
      #tags$img(src="Sample_data.png",width = "500px", height = "100px"),
      #tags$div("The data matrix file should be named Data.txt",style="color:blue"),
      #tags$strong("sample data-key file"),
      #br(),
      #tags$img(src="Sample_Key.png",width = "300px", height = "100px"),
      #tags$div("The samples-key file should be named Samples.txt",style="color:blue"),
      #tags$div("the sample groups should be of same number",style="color:red")
      
      
    )
  )
)
options(shiny.maxRequestSize = 30*1024^2)

server <- function(input, output) {
  output$fileob <- renderTable({
    file <- input$file
    ext1 <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext1 == c("txt","tsv","csv"), "Please upload a tab seperated file"))
    nfiles = nrow(input$file) 
    inFile = list()
    #for (i in 1 : nfiles)
    #{
    #if(input$file$name[i] == "Data.txt"){
    #inFile[[i]] <- read.table(input$file$datapath[i],sep = "\t",header = TRUE)
    #df <- inFile[[i]]
    #}else{
    #inFile[[i]] = read.table(input$file$datapath[i],sep = "\t",header = TRUE)
    #sam <- inFile[[i]]
    #}                  
    #}
    #inFile = ""
    if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
      main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
      #main <- inFile[[1]]
      for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
    }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
      main <- read.csv(input$file$datapath[1],header = TRUE)
      #main <- inFile[[1]]
      for_sam <- read.csv(input$file$datapath[1],header = FALSE)
    }
    sub_main <- select(main,matches("_|Gene"))
    df <- sub_main[2:nrow(sub_main),]
    firstRow <- sub_main[1,1:ncol(sub_main)]
    df[is.na(df)] <- 0
    mati <- df[,2:ncol(df)]
    mati <- as.data.frame(sapply(mati, as.numeric))
    for(i in seq_along(mati)) {
      newVal <- as.data.frame(mati[,i]) %>% filter(as.data.frame(mati[,i]) > 0)
      mati[[i]][mati[[i]] == 0] <- runif(nrow(mati), min = 0.0001, max=min(newVal))
    }
    NewData <-cbind(df$Gene,mati)
    colnames(NewData)[1] <- "Gene"
    all_new <- rbind(firstRow,NewData)
    write.table(all_new,"Data_imput.txt",sep = "\t",quote = F,row.names = FALSE)
    print(file)
  })
  imputdownload <- reactive({
    imputdat <- read.table("Data_imput.txt",header = TRUE,sep = "\t")
    print(imputdat)
  })
  testData <- reactive({
    df <- read.csv("testData.csv",header = TRUE)
    print(df)
  })
  #testSam <- reactive({
  #sam <- read.table("testSamples.txt",sep = "\t",header = TRUE)
  #print(sam)
  #})
  button1 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      #df <- ""
      #sam <- ""
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      print(sub_match_sam)
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      #samGroup <- names(sam)
      print(names(sam))
    }
  })
  summary <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df <- as.data.frame(sapply(df, as.numeric))
      print(head(df))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      type <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
        type <- read.table("Type_scoring_Hum.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
        type <- read.table("Type_scoring_Mus.txt",sep = "\t",header = T)
      }
      #ann <- read.table("Class_annot.txt",sep = "\t",header = T)
      #type <- read.table("Type_scoring.txt",sep = "\t",header = T)
      Nsam <- input$Group
      print(Nsam)
      sumlist <- list()
      dat <- 0
      rbc_sub_col1Mean <- 0
      rbc_sub_col2Mean <- 0
      plasma_sub_col1Mean <- 0
      plasma_sub_col2Mean <- 0
      coreMat_sub_col1Mean <- 0
      coreMat_sub_col2Mean <- 0
      glycECM_sub_col1Mean <- 0
      glycECM_sub_col2Mean <- 0
      protGly_sub_col1Mean <- 0
      protGly_sub_col2Mean <- 0
      Coll_sub_col1Mean <- 0
      Coll_sub_col2Mean <- 0
      regECM_sub_col1Mean <- 0
      regECM_sub_col2Mean <- 0
      affECM_sub_col1Mean <- 0
      affECM_sub_col2Mean <- 0
      secfac_sub_col1Mean <- 0
      secfac_sub_col2Mean <- 0
      for(i in 1:length(Nsam)){
        sam1 <- as.data.frame(sam[[Nsam[[i]]]])
        colnames(sam1) <- "samp"
        print(Nsam[[i]])
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp) | sam1$samp == ""), ])
        colnames(sam1) <- "samp"
        grp <- as.data.frame(df %>% select(one_of(dput(as.character(sam1$samp)))))
        grp[is.na(grp)] = 0
        n1 <- nrow(sam1)
        print(n1)
        val <- 0
        for(j in 1:n1) {
          val[j] <- sum(grp[,j])
        }
        High <- max(val)
        for(k in 1:ncol(grp)) {
          tot <- sum(grp[,k])
          grp[,k] <- (grp[,k]*High)/tot
        }
        mean <- as.data.frame(rowMeans(grp))
        df1 <- cbind(dfGene,mean,grp)
        print(head(dfGene))
        df2 <- filter(df1,df1$`rowMeans(grp)` > 0)
        #annot <- merge(df2, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE
        annot <- merge(df2, type, by.x = "df$Gene", by.y = "Gene.names",all.x = TRUE)
        #ecm <- filter(annot,annot$Division == "ECM")
        #matri <- filter(annot,annot$Division == "Secreted-matrisome associated")
        #plas <- filter(annot,annot$Division == "Secreted plasma")
        #df3 <- as.data.frame(paste(c(Nsam[[i]],length(dat[[i]]),nrow(df2),nrow(ecm),nrow(matri),nrow(plas)),sep = "\t"))
        #val <- data.frame(c("Sample","# of replicates","# of genes","Total ECM proteins","Total matrisome-secreted","Total secreted-plasma"))
        #colnames(df3) <- "stat"
        #colnames(val) <- "summary"
        #df4 <- cbind(val,df3)
        #sumlist[[i]] <- df4
        hem_chain <- ""
        if(input$species == "Human"){
          hem_chain <- filter(annot,annot$`df$Gene` == "HBB" | annot$`df$Gene` == "HBA1" | annot$`df$Gene` == "HBG1" | annot$`df$Gene` == "HBG2" | annot$`df$Gene` == "HBA2" | annot$`df$Gene` == "HBD" | annot$`df$Gene` == "HBE1" | annot$`df$Gene` == "HBZ" | annot$`df$Gene` == "HBQ1" | annot$`df$Gene` == "HBM")
        }else if(input$species == "Mouse"){
          hem_chain <- filter(annot,annot$`df$Gene` == "Hbb-b1" | annot$`df$Gene` == "Hbb-b2" | annot$`df$Gene` == "Hbb" | annot$`df$Gene` == "Hbb-y" | annot$`df$Gene` == "Hbb-bh1" | annot$`df$Gene` == "Hba-a1" | annot$`df$Gene` == "Hba" | annot$`df$Gene` == "Hba-x" | annot$`df$Gene` == "Hbb-bt" | annot$`df$Gene` == "Hba-a2" | annot$`df$Gene` == "Hbb-bs")
        }
        #rbc <- filter(annot,annot$RBC == "RBC")
        rbc <- filter(annot,annot$RBC == "RBC" & annot$`df$Gene` != "Hbb-b1" & annot$`df$Gene` != "Hbb-b2" & annot$`df$Gene` != "Hbb" & annot$`df$Gene` != "Hbb-y" & annot$`df$Gene` != "Hbb-bh1" & annot$`df$Gene` != "Hba-a1" & annot$`df$Gene` != "Hba" & annot$`df$Gene` != "Hba-x" & annot$`df$Gene` != "Hbb-bt" & annot$`df$Gene` != "Hba-a2" & annot$`df$Gene` != "Hbb-bs")
        plasma <- filter(annot,annot$Plasma == "Plasma")
        coreMat <- filter(annot,annot$Matrisome.1 == "Core matrisome")
        glycECM <- filter(annot,annot$Matrisome.2 == "ECM Glycoproteins")
        protGly <- filter(annot,annot$Matrisome.2 == "Proteoglycans")
        Coll <- filter(annot,annot$Matrisome.2 == "Collagens")
        regECM <- filter(annot,annot$Matrisome.2 == "ECM Regulators")
        affECM <- filter(annot,annot$Matrisome.2 == "ECM-affiliated Proteins")
        secfac <- filter(annot,annot$Matrisome.2 == "Secreted Factors")
        bm <- filter(annot,annot$Plasma.1 == "BM")
        protease <- filter(annot,annot$Matrisome.3 == "Protease")
        cross <- filter(annot,annot$Matrisome.3 == "Crosslinking Enzyme")
        comp <- filter(annot,annot$Plasma.1 == "comp")
        Col1 <- ""
        Col2 <- ""
        Col3 <- ""
        Col5 <- ""
        if(input$species == "Human"){
          Col1 <-  filter(annot,annot$`df$Gene` == "COL1A1")
          Col2 <- filter(annot,annot$`df$Gene` == "COL1A2")
          Col3 <- filter(annot,annot$`df$Gene` == "COL3A1")
          Col5 <- filter(annot,annot$`df$Gene` == "COL6A5")  
        }else if(input$species == "Mouse"){
          Col1 <-  filter(annot,annot$`df$Gene` == "Col1a1")
          Col2 <- filter(annot,annot$`df$Gene` == "Col1a2")
          Col3 <- filter(annot,annot$`df$Gene` == "Col3a1")
          Col5 <- filter(annot,annot$`df$Gene` == "Col6a5")
        }
        rat_col1_2 <- 0
        rat_col3_1 <- 0
        rat_col5_1 <- 0
        if(nrow(Col1) > 0 && Col1$`rowMeans(grp)`> 0){
          if(nrow(Col2) > 0 && Col2$`rowMeans(grp)`> 0){
            rat_col1_2 <- Col1$`rowMeans(grp)`/Col2$`rowMeans(grp)`
          }
        }
        if(nrow(Col3) > 0 && Col3$`rowMeans(grp)`> 0){
          if(nrow(Col1) > 0 && Col1$`rowMeans(grp)`> 0){
            rat_col3_1 <- Col3$`rowMeans(grp)`/Col1$`rowMeans(grp)`
          }
        }
        if(nrow(Col5) > 0 && Col5$`rowMeans(grp)`> 0){
          if(nrow(Col1) > 0 && Col1$`rowMeans(grp)`> 0){
            rat_col5_1 <- Col5$`rowMeans(grp)`/Col1$`rowMeans(grp)`
          }
        }
        #Srbc <- round(nrow(rbc)/2,0)
        Srbc <- 20
        Smat <- round(nrow(coreMat)/2,0)
        SglyE <- round(nrow(glycECM)/2,0)
        SprotG <- round(nrow(protGly)/2,0)
        Scol <- round(nrow(Coll)/2,0)
        SregE <- round(nrow(regECM)/2,0)
        SaffE <- round(nrow(affECM)/2,0)
        Ssec <- round(nrow(secfac)/2,0)
        bmS <- round(nrow(bm)/2,0)
        protS <- round(nrow(protease)/2,0)
        Scross <- round(nrow(cross)/2,0)
        compS <- round(nrow(comp)/2,0)
        Splas <- 150
        rbc_sub <- head(rbc[order(rbc$`rowMeans(grp)`, decreasing=TRUE), ], Srbc)
        plasma_sub <- head(plasma[order(plasma$`rowMeans(grp)`, decreasing=TRUE), ], Splas)
        coreMat_sub <- head(coreMat[order(coreMat$`rowMeans(grp)`, decreasing=TRUE), ], Smat)
        glycECM_sub <-head(glycECM[order(glycECM$`rowMeans(grp)`, decreasing=TRUE), ], SglyE)
        protGly_sub <- head(protGly[order(protGly$`rowMeans(grp)`, decreasing=TRUE), ], SprotG)
        Coll_sub <- head(Coll[order(protGly$`rowMeans(grp)`, decreasing=TRUE), ], Scol)
        regECM_sub <- head(regECM[order(regECM$`rowMeans(grp)`, decreasing=TRUE), ], SregE)
        affECM_sub <- head(affECM[order(affECM$`rowMeans(grp)`, decreasing=TRUE), ], SaffE)
        secfac_sub <- head(secfac[order(secfac$`rowMeans(grp)`, decreasing=TRUE), ], Ssec)
        bm_sub <- head(bm[order(bm$`rowMeans(grp)`, decreasing=TRUE), ], bmS)
        Prot_sub <- head(protease[order(protease$`rowMeans(grp)`, decreasing=TRUE), ], protS)
        cross_sub <- head(cross[order(cross$`rowMeans(grp)`, decreasing=TRUE), ], Scross)
        comp_sub <- head(comp[order(comp$`rowMeans(grp)`, decreasing=TRUE), ], compS)
        #rbc_sub_col1Mean <- 0;
        #rbc_sub_col2Mean <- 0;
        if(i == 1){
          for(l in 1:nrow(rbc_sub)) {
            rbc_sub_col1 <- rbc_sub[,3:(ncol(rbc_sub)-7)]
            rbc_sub_col1Mean <- colMeans(rbc_sub_col1)
            rbc_sub_col1Mean[is.na(rbc_sub_col1Mean)] = 0
            rbc_sub_col1Mean <- as.numeric(rbc_sub_col1Mean)
          }
          for(l in 1:nrow(plasma_sub)) {
            plasma_sub_col1 <- plasma_sub[,3:(ncol(plasma_sub)-7)]
            plasma_sub_col1Mean <- colMeans(plasma_sub_col1)
            plasma_sub_col1Mean[is.na(plasma_sub_col1Mean)] = 0
            plasma_sub_col1Mean <- as.numeric(plasma_sub_col1Mean)
          }
          for(l in 1:nrow(coreMat_sub)) {
            coreMat_sub_col1 <- coreMat_sub[,3:(ncol(coreMat_sub)-7)]
            coreMat_sub_col1Mean <- colMeans(coreMat_sub_col1)
            coreMat_sub_col1Mean[is.na(coreMat_sub_col1Mean)] = 0
            coreMat_sub_col1Mean <- as.numeric(coreMat_sub_col1Mean)
          }
          for(l in 1:nrow(glycECM_sub)) {
            glycECM_sub_col1 <- glycECM_sub[,3:(ncol(glycECM_sub)-7)]
            glycECM_sub_col1Mean <- colMeans(glycECM_sub_col1)
            glycECM_sub_col1Mean[is.na(glycECM_sub_col1Mean)] = 0
            glycECM_sub_col1Mean <- as.numeric(glycECM_sub_col1Mean)
          }
          for(l in 1:nrow(protGly_sub)) {
            protGly_sub_col1 <- protGly_sub[,3:(ncol(protGly_sub)-7)]
            protGly_sub_col1Mean <- colMeans(protGly_sub_col1)
            protGly_sub_col1Mean[is.na(protGly_sub_col1Mean)] = 0
            protGly_sub_col1Mean <- as.numeric(protGly_sub_col1Mean)
          }
          for(l in 1:nrow(Coll_sub)) {
            Coll_sub_col1 <- Coll_sub[,3:(ncol(Coll_sub)-7)]
            Coll_sub_col1Mean <- colMeans(Coll_sub_col1)
            Coll_sub_col1Mean[is.na(Coll_sub_col1Mean)] = 0
            Coll_sub_col1Mean <- as.numeric(Coll_sub_col1Mean)
          }
          for(l in 1:nrow(regECM_sub)) {
            regECM_sub_col1 <- regECM_sub[,3:(ncol(regECM_sub)-7)]
            regECM_sub_col1Mean <- colMeans(regECM_sub_col1)
            regECM_sub_col1Mean[is.na(regECM_sub_col1Mean)] = 0
            regECM_sub_col1Mean <- as.numeric(regECM_sub_col1Mean)
          }
          for(l in 1:nrow(affECM_sub)) {
            affECM_sub_col1 <- affECM_sub[,3:(ncol(affECM_sub)-7)]
            affECM_sub_col1Mean <- colMeans(affECM_sub_col1)
            affECM_sub_col1Mean[is.na(affECM_sub_col1Mean)] = 0
            affECM_sub_col1Mean <- as.numeric(affECM_sub_col1Mean)
          }
          for(l in 1:nrow(secfac_sub)) {
            secfac_sub_col1 <- secfac_sub[,3:(ncol(secfac_sub)-7)]
            secfac_sub_col1Mean <- colMeans(secfac_sub_col1)
            secfac_sub_col1Mean[is.na(secfac_sub_col1Mean)] = 0
            secfac_sub_col1Mean <- as.numeric(secfac_sub_col1Mean)
          }
          for(l in 1:nrow(bm_sub)) {
            bm_sub_col1 <- bm_sub[,3:(ncol(bm_sub)-7)]
            bm_sub_col1Mean <- colMeans(bm_sub_col1)
            bm_sub_col1Mean[is.na(bm_sub_col1Mean)] = 0
            bm_sub_col1Mean <- as.numeric(bm_sub_col1Mean)
          }
          for(l in 1:nrow(Prot_sub)) {
            Prot_sub_col1 <- Prot_sub[,3:(ncol(Prot_sub)-7)]
            Prot_sub_col1Mean <- colMeans(Prot_sub_col1)
            Prot_sub_col1Mean[is.na(Prot_sub_col1Mean)] = 0
            Prot_sub_col1Mean <- as.numeric(Prot_sub_col1Mean)
          }
          for(l in 1:nrow(cross_sub)) {
            cross_sub_col1 <- cross_sub[,3:(ncol(cross_sub)-7)]
            cross_sub_col1Mean <- colMeans(cross_sub_col1)
            cross_sub_col1Mean[is.na(cross_sub_col1Mean)] = 0
            cross_sub_col1Mean <- as.numeric(cross_sub_col1Mean)
          }
          for(l in 1:nrow(comp_sub)) {
            comp_sub_col1 <- comp_sub[,3:(ncol(comp_sub)-7)]
            comp_sub_col1Mean <- colMeans(comp_sub_col1)
            comp_sub_col1Mean[is.na(comp_sub_col1Mean)] = 0
            comp_sub_col1Mean <- as.numeric(comp_sub_col1Mean)
          }
        }else if(i == 2){
          for(l in 1:nrow(rbc_sub)) {
            rbc_sub_col2 <- rbc_sub[,3:(ncol(rbc_sub)-7)]
            rbc_sub_col2Mean <- colMeans(rbc_sub_col2)
            rbc_sub_col2Mean[is.na(rbc_sub_col2Mean)] = 0
            rbc_sub_col2Mean <- as.numeric(rbc_sub_col2Mean)
          }
          for(l in 1:nrow(plasma_sub)) {
            plasma_sub_col2 <- plasma_sub[,3:(ncol(plasma_sub)-7)]
            plasma_sub_col2Mean <- colMeans(plasma_sub_col2)
            plasma_sub_col2Mean[is.na(plasma_sub_col2Mean)] = 0
            plasma_sub_col2Mean <- as.numeric(plasma_sub_col2Mean)
          }
          for(l in 1:nrow(coreMat_sub)) {
            coreMat_sub_col2 <- coreMat_sub[,3:(ncol(coreMat_sub)-7)]
            coreMat_sub_col2Mean <- colMeans(coreMat_sub_col2)
            coreMat_sub_col2Mean[is.na(coreMat_sub_col2Mean)] = 0
            coreMat_sub_col2Mean <- as.numeric(coreMat_sub_col2Mean)
          }
          for(l in 1:nrow(glycECM_sub)) {
            glycECM_sub_col2 <- glycECM_sub[,3:(ncol(glycECM_sub)-7)]
            glycECM_sub_col2Mean <- colMeans(glycECM_sub_col2)
            glycECM_sub_col2Mean[is.na(glycECM_sub_col2Mean)] = 0
            glycECM_sub_col2Mean <- as.numeric(glycECM_sub_col2Mean)
          }
          for(l in 1:nrow(protGly_sub)) {
            protGly_sub_col2 <- protGly_sub[,3:(ncol(protGly_sub)-7)]
            protGly_sub_col2Mean <- colMeans(protGly_sub_col2)
            protGly_sub_col2Mean[is.na(protGly_sub_col2Mean)] = 0
            protGly_sub_col2Mean <- as.numeric(protGly_sub_col2Mean)
          }
          for(l in 1:nrow(Coll_sub)) {
            Coll_sub_col2 <- Coll_sub[,3:(ncol(Coll_sub)-7)]
            Coll_sub_col2Mean <- colMeans(Coll_sub_col2)
            Coll_sub_col2Mean[is.na(Coll_sub_col2Mean)] = 0
            Coll_sub_col2Mean <- as.numeric(Coll_sub_col2Mean)
          }
          for(l in 1:nrow(regECM_sub)) {
            regECM_sub_col2 <- regECM_sub[,3:(ncol(regECM_sub)-7)]
            regECM_sub_col2Mean <- colMeans(regECM_sub_col2)
            regECM_sub_col2Mean[is.na(regECM_sub_col2Mean)] = 0
            regECM_sub_col2Mean <- as.numeric(regECM_sub_col2Mean)
          }
          for(l in 1:nrow(affECM_sub)) {
            affECM_sub_col2 <- affECM_sub[,3:(ncol(affECM_sub)-7)]
            affECM_sub_col2Mean <- colMeans(affECM_sub_col2)
            affECM_sub_col2Mean[is.na(affECM_sub_col2Mean)] = 0
            affECM_sub_col2Mean <- as.numeric(affECM_sub_col2Mean)
          }
          for(l in 1:nrow(secfac_sub)) {
            secfac_sub_col2 <- secfac_sub[,3:(ncol(secfac_sub)-7)]
            secfac_sub_col2Mean <- colMeans(secfac_sub_col2)
            secfac_sub_col2Mean[is.na(secfac_sub_col2Mean)] = 0
            secfac_sub_col2Mean <- as.numeric(secfac_sub_col2Mean)
          }
          for(l in 1:nrow(bm_sub)) {
            bm_sub_col2 <- bm_sub[,3:(ncol(bm_sub)-7)]
            bm_sub_col2Mean <- colMeans(bm_sub_col2)
            bm_sub_col2Mean[is.na(bm_sub_col2Mean)] = 0
            bm_sub_col2Mean <- as.numeric(bm_sub_col2Mean)
          }
          for(l in 1:nrow(Prot_sub)) {
            Prot_sub_col2 <- Prot_sub[,3:(ncol(Prot_sub)-7)]
            Prot_sub_col2Mean <- colMeans(Prot_sub_col2)
            Prot_sub_col2Mean[is.na(Prot_sub_col2Mean)] = 0
            Prot_sub_col2Mean <- as.numeric(Prot_sub_col2Mean)
          }
          for(l in 1:nrow(cross_sub)) {
            cross_sub_col2 <- cross_sub[,3:(ncol(cross_sub)-7)]
            cross_sub_col2Mean <- colMeans(cross_sub_col2)
            cross_sub_col2Mean[is.na(cross_sub_col2Mean)] = 0
            cross_sub_col2Mean <- as.numeric(cross_sub_col2Mean)
          }
          for(l in 1:nrow(comp_sub)) {
            comp_sub_col2 <- comp_sub[,3:(ncol(comp_sub)-7)]
            comp_sub_col2Mean <- colMeans(comp_sub_col2)
            comp_sub_col2Mean[is.na(comp_sub_col2Mean)] = 0
            comp_sub_col2Mean <- as.numeric(comp_sub_col2Mean)
          }
        }
        totSum <- sum(df1$`rowMeans(grp)`,na.rm = TRUE)
        hem_chain_score <- sum(hem_chain$`rowMeans(grp)`,na.rm = "TRUE")*10
        #rbc_score <- sum(rbc_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5/totSum
        rbc_score <- (hem_chain_score + sum(rbc_sub$`rowMeans(grp)`,na.rm = "TRUE")*5)/totSum
        plasma_score <- (sum(plasma_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        coreMat_score <- (sum(coreMat_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        glycECM_score <- (sum(glycECM_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        protGly_score <- (sum(protGly_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        Coll_score <- (sum(Coll_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        regECM_score <- (sum(regECM_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        affECM_score <- (sum(affECM_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        secfac_score <- (sum(secfac_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        bm_score <- (sum(bm_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        Prot_score <- (sum(Prot_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        cross_score <- (sum(cross_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        comp_score <- (sum(comp_sub$`rowMeans(grp)`,na.rm = "TRUE")*2.5)/totSum
        print(Coll_score)
        df3 <- as.data.frame(paste(c(Nsam[[i]],round(rbc_score,5),round(plasma_score,5),round(coreMat_score,5),round(glycECM_score,5),round(protGly_score,5),round(Coll_score,5),round(regECM_score,5),round(affECM_score,5),round(secfac_score,5),round(bm_score,5),round(Prot_score,5),round(cross_score,5),round(comp_score,5),round(rat_col1_2,5),round(rat_col3_1,5),round(rat_col5_1,5)),sep = "\t"))
        #val <- data.frame(c("Sample","# of replicates","# of genes","RBC score","Plasma score","Core Matrisome score"))
        colnames(df3) <- Nsam[[i]]
        #colnames(val) <- "summary"
        #df4 <- cbind(val,df3)
        sumlist[[i]] <- t(df3)
      }
      print(protGly_sub_col1Mean)
      rbc_sub_stats <- t.test(rbc_sub_col1Mean,rbc_sub_col2Mean,paired = FALSE)
      rbc_sub_pval <- rbc_sub_stats$p.value
      plasma_sub_stats <- t.test(plasma_sub_col1Mean,plasma_sub_col2Mean,paired = FALSE)
      plasma_sub_pval <- plasma_sub_stats$p.value
      coreMat_sub_stats <- t.test(coreMat_sub_col1Mean,coreMat_sub_col2Mean,paired = FALSE)
      coreMat_sub_pval <- coreMat_sub_stats$p.value
      glycECM_sub_stats <- t.test(glycECM_sub_col1Mean,glycECM_sub_col2Mean,paired = FALSE)
      glycECM_sub_pval <- glycECM_sub_stats$p.value
      protGly_sub_stats <- t.test(protGly_sub_col1Mean,protGly_sub_col2Mean,paired = FALSE)
      protGly_sub_pval <- protGly_sub_stats$p.value
      Coll_sub_stats <- t.test(Coll_sub_col1Mean,Coll_sub_col2Mean,paired = FALSE)
      Coll_sub_pval <- Coll_sub_stats$p.value
      regECM_sub_stats <- t.test(regECM_sub_col1Mean,regECM_sub_col2Mean,paired = FALSE)
      regECM_sub_pval <- regECM_sub_stats$p.value
      affECM_sub_stats <- t.test(affECM_sub_col1Mean,affECM_sub_col2Mean,paired = FALSE)
      affECM_sub_pval <- affECM_sub_stats$p.value
      secfac_sub_stats <- t.test(secfac_sub_col1Mean,secfac_sub_col2Mean,paired = FALSE)
      secfac_sub_pval <- secfac_sub_stats$p.value
      bm_sub_stats <- t.test(bm_sub_col1Mean,bm_sub_col2Mean,paired = FALSE)
      bm_sub_pval <- bm_sub_stats$p.value
      Prot_sub_stats <- t.test(Prot_sub_col1Mean,Prot_sub_col2Mean,paired = FALSE)
      Prot_sub_pval <- Prot_sub_stats$p.value
      cross_sub_stats <- t.test(cross_sub_col1Mean,cross_sub_col2Mean,paired = FALSE)
      cross_sub_pval <- cross_sub_stats$p.value
      comp_sub_stats <- t.test(comp_sub_col1Mean,comp_sub_col2Mean,paired = FALSE)
      comp_sub_pval <- comp_sub_stats$p.value
      pval <- data.frame(c("pvalue",round(rbc_sub_pval,10),round(plasma_sub_pval,10),round(coreMat_sub_pval,10),round(glycECM_sub_pval,10),round(protGly_sub_pval,10),round(Coll_sub_pval,10),round(regECM_sub_pval,10),round(affECM_sub_pval,10),round(secfac_sub_pval,10),round(bm_sub_pval,10),round(Prot_sub_pval,10),round(cross_sub_pval,10),round(comp_sub_pval,10),"NA","NA","NA"))
      big_data = do.call(rbind, sumlist)
      val <- data.frame(c("Sample","RBC score","Plasma score","Core ECM score","ECM glycoprotein score","Proteoglycan score","Collagen score","ECM regulator score","ECM-affiliated protein score","Secretor factor score","Basement Membrane score","Protease score","Crosslinking enzyme score","Complement score","COL1A1/COL1A2 ratio","COL3A1/COL1A1 ratio","COL6A5/COL1A1 ratio"))
      #first <- paste(c("Total samples: ",length(Nsam)),sep = "\t")
      #print(rbind(first,big_data))
      all_val <- cbind(val,t(big_data),pval)
      names(all_val) <- as.matrix(all_val[1, ])
      all_val <- all_val[-1, ]
      all_val[] <- lapply(all_val, function(x) type.convert(as.character(x)))
      row.names(all_val) <- NULL
      print(all_val)
    }
  })
  var <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      #df <- inFile[[1]]
      #sam <- inFile[[2]]
      df <- arrange(df, Gene, desc(Gene))
      print(df$Gene)
    }
  })
  division <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      #df <- ""
      #sam <- ""
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      df <- arrange(df, Gene, desc(Gene))
      df1 <- merge(df, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
      df1[is.na(df1)] = "Unknown"
      divUniq <- as.data.frame(unique(df1$Division))
      colnames(divUniq) <- "Div"
      divUniq <- arrange(divUniq,Div,desc(Div))
      print(divUniq)
    }
  })
  
  output$Species <- renderUI({
    if(is.null(input$file)) {
      return() 
    }else{
      radioButtons("species","Select the species",choices = c("Human" = "Human","Mouse" = "Mouse"),selected = "Human")
    }
  })
  output$imput <- renderUI({
    if(is.null(input$file)) {
      return() 
    }else{
      radioButtons("imputate","Is imputation needed?",choices = c("Yes" = "Yes","No" = "No"),selected = "No")
    }
  })
  output$sampleGroup <- renderUI({
    if(is.null(input$file)) {
      return() 
    }else{
      checkboxGroupInput("Group","Select the two Sample group", choices = button1())
    }
  })
  Analysis <- reactive({
    
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      print(sam1)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(dfGene,mean1,mean2)
      attach(input$avg_val)
      if(input$imputate == "No"){
        df1 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        df1[is.na(df1)] = "Unknown"
        print(nrow(df1))
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        print(nrow(df))
        log2FC <- cbind(df1,log2Rat)
        df2 <- cbind(dfGene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df2)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else if(n1 != n2){
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value       
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        names(log2FC)[names(log2FC) == 'df$Gene'] <- 'Protein'
        logFC1 <- cbind(log2FC$Protein,log2FC$Division,log2FC$Category,log2FC$Secretome,log2FC$log2Rat)
        colnames(logFC1) <- c("Protein","Division","Category","Secretome","log2Ratio")
        tab <- cbind(logFC1,pval)
        print(tab)
      }else if(input$imputate == "Yes"){
        sub_main <- read.table("Data_imput.txt",sep = "\t",header = TRUE)
        df <- sub_main[2:nrow(sub_main),]
        dfGene <- as.data.frame(df$Gene)
        df <- arrange(df, Gene, desc(Gene))
        #Idf <- ndf[,1:ncol(ndf)]
        samN <- input$Group
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        print(sam1)
        df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        mean2 <- rowMeans(grp2)
        mean1 <- as.data.frame(mean1)
        mean2 <- as.data.frame(mean2)
        avg_val <- cbind(dfGene,mean1,mean2)
        df1 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        df1[is.na(df1)] = "Unknown"
        print(nrow(df1))
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        print(nrow(df))
        log2FC <- cbind(df1,log2Rat)
        df2 <- cbind(dfGene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df2)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else if(n1 != n2){
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value       
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        names(log2FC)[names(log2FC) == 'df$Gene'] <- 'Protein'
        logFC1 <- cbind(log2FC$Protein,log2FC$Division,log2FC$Category,log2FC$Secretome,log2FC$log2Rat)
        colnames(logFC1) <- c("Protein","Division","Category","Secretome","log2Ratio")
        tab <- cbind(logFC1,pval)
        print(tab)
      }
    }
  })
  output$newdata1 <- DT::renderDataTable({
    DT::datatable(Analysis(),options = list(lengthMenu = c(10, 50, 100), pageLength = 10))
  })
  Sgene <- reactive({
    req(input$newdata1_rows_selected)
    selRow <- Analysis()[input$newdata1_rows_selected,]
    print(selRow)
    #write.table(selRow,"selected_gene_labels.txt",sep = "\t",quote = FALSE,row.names = FALSE)
  })
  output$div <- renderUI({
    if(is.null(input$file)) {
      return() 
    }else{
      selectInput("level","Select the division for Volcano Plot",choices = division())
    }
  })
  Plotting1 <- reactive({
    
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(dfGene,mean1,mean2)
      attach(input$avg_val)
      if(input$imputate == "No"){
        status <- exists("sel")
        #Prot <- sel$protein
        #print(Prot)
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(dfGene,log2Rat)
        df1 <- cbind(dfGene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        #volc <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        #volc+geom_text_repel(data=subset(newTab1,log10pvalue>=1.3 & abs(log2FC)>=0.6),aes(label=Gene),max.overlaps = 100)
        volc1 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        #volc3 <- volc1+geom_text_repel(data=subset(newTab1,log10pvalue>=1.3 & abs(log2FC)>=0.6),aes(label=Gene),max.overlaps = 1000)
        #volc3 <- volc1+geom_text_repel(data=subset(newTab1,Gene=sel),aes(label=Gene),max.overlaps = 1000)
        newTab1$group = ifelse(newTab1$Division == input$level, "Selected Division","Other Division" )
        volc2 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=group))+scale_color_manual(values=c("black","green"),labels = c("Other divisions", "Selected Division"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        volc4 <- volc2+geom_text_repel(data=subset(newTab1,Division == input$level),aes(label=Gene),max.overlaps = 1000)
        ggarrange(volc1, volc4, 
                  ncol = 2, nrow = 1)
      }else if(input$imputate == "Yes"){
        sub_main <- read.table("Data_imput.txt",sep = "\t",header = TRUE)
        df <- sub_main[2:nrow(sub_main),]
        dfGene <- as.data.frame(df$Gene)
        df <- arrange(df, Gene, desc(Gene))
        samN <- input$Group
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        mean2 <- rowMeans(grp2)
        mean1 <- as.data.frame(mean1)
        mean2 <- as.data.frame(mean2)
        avg_val <- cbind(dfGene,mean1,mean2)
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(dfGene,log2Rat)
        df1 <- cbind(dfGene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        #volc <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        #volc+geom_text_repel(data=subset(newTab1,log10pvalue>=1.3 & abs(log2FC)>=0.6),aes(label=Gene),max.overlaps = 100)
        volc1 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        #volc3 <- volc1+geom_text_repel(data=subset(newTab1,log10pvalue>=1.3 & abs(log2FC)>=0.6),aes(label=Gene),max.overlaps = 1000)
        newTab1$group = ifelse(newTab1$Division == input$level, "Selected Division","Other Division" )
        volc2 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=group))+scale_color_manual(values=c("black","green"),labels = c("Other divisions", "Selected Division"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        volc4 <- volc2+geom_text_repel(data=subset(newTab1,Division == input$level),aes(label=Gene),max.overlaps = 1000)
        ggarrange(volc1, volc4, 
                  ncol = 2, nrow = 1)
      }
      
    }
  })
  Plotting1_1 <- reactive({
    
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      df <- arrange(df, Gene, desc(Gene))
      df[is.na(df)] <- 0
      dfGene <- as.data.frame(df$Gene)
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(df$Gene,mean1,mean2)
      attach(input$avg_val)
      if(input$imputate == "No"){
        status <- exists("sel")
        #Prot <- sel$protein
        #print(Prot)
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(dfGene,log2Rat)
        df1 <- cbind(dfGene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        sel <- as.data.frame(Sgene())
        Prot <- sel$Protein
        print(Prot)
        newTab1$genelabels <- factor(newTab1$Gene,levels = Prot)
        volc1 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        volc3 <- volc1+geom_text_repel(data=newTab1,aes(label=genelabels),col = "black", na.rm = TRUE, box.padding = unit(0.45, "lines"), hjust = 1,max.overlaps = 1000)
        #newTab1$group = ifelse(newTab1$Division == input$level, "Selected Division","Other Division" )
        #volc2 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=group))+scale_color_manual(values=c("black","green"),labels = c("Other divisions", "Selected Division"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        #volc4 <- volc2+geom_text_repel(data=subset(newTab1,Division == input$level),aes(label=Gene),max.overlaps = 1000)
        ggarrange(volc1, volc3,ncol = 2, nrow = 1)
        
      }else if(input$imputate == "Yes"){
        sub_main <- read.table("Data_imput.txt",sep = "\t",header = TRUE)
        df <- sub_main[2:nrow(sub_main),]
        dfGene <- as.data.frame(df$Gene)
        df <- arrange(df, Gene, desc(Gene))
        samN <- input$Group
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        mean2 <- rowMeans(grp2)
        mean1 <- as.data.frame(mean1)
        mean2 <- as.data.frame(mean2)
        avg_val <- cbind(df$Gene,mean1,mean2)
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(df$Gene,log2Rat)
        df1 <- cbind(df$Gene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        sel <- as.data.frame(Sgene())
        Prot <- sel$Protein
        print(Prot)
        newTab1$genelabels <- factor(newTab1$Gene,levels = Prot)
        #volc <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        #volc+geom_text_repel(data=subset(newTab1,log10pvalue>=1.3 & abs(log2FC)>=0.6),aes(label=Gene),max.overlaps = 100)
        volc1 <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=threshold))+scale_color_manual(values=c("black","red"),labels = c("Non significant", "Significant"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))+ggtitle(paste(samN[[2]],"_vs_",samN[[1]],sep = ""))
        volc3 <- volc1+geom_text_repel(data=newTab1,aes(label=genelabels),col = "black", na.rm = TRUE, box.padding = unit(0.45, "lines"), hjust = 1)
        ggarrange(volc1, volc3, ncol = 2, nrow = 1)
        
      }
      
    }
  })
  Plotting2 <- reactive({
    
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(df$Gene,mean1,mean2)
      if(input$imputate == "No"){
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(df$Gene,log2Rat)
        df1 <- cbind(df$Gene,grp1,grp2)
        pval <- 0
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        #newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        newTab1$group = ifelse(newTab1$Division == input$level, "Selected Division","Other Division" )
        volc <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=group))+scale_color_manual(values=c("black","green"),labels = c("Other divisions", "Selected Division"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        volc+geom_text_repel(data=subset(newTab1,Division == input$level),aes(label=Gene),max.overlaps = 1000)
      }else if(input$imputate == "Yes"){
        newMean1 <- mean1 %>% filter(mean1 > 0)
        newMean2 <- mean2 %>% filter(mean2 > 0)
        #sumStat1 <- summary(newMean1$mean1)
        #sumStat2 <- summary(newMean2$mean2)
        setDT(avg_val)[mean1 == 0, mean1 := runif(.N, min=0.0001, max=min(newMean1))]
        setDT(avg_val)[mean2 == 0, mean2 := runif(.N, min=0.0001, max=min(newMean2))]
        avg_val <- as.data.frame(avg_val)
        log2Rat <- log2(avg_val$mean2/avg_val$mean1)
        log2Rat[log2Rat == "NaN"] <- 0
        log2Rat[log2Rat == "-Inf"] <- 0
        log2Rat[log2Rat == "Inf"] <- 0
        log2Rat <- as.data.frame(log2Rat)
        log2FC <- cbind(df$Gene,log2Rat)
        df1 <- cbind(df$Gene,grp1,grp2)
        pval <- 0
        for(i in seq_along(grp1)) {
          grp1[[i]][grp1[[i]] == 0] <- runif(nrow(grp1), min = 0.0001, max=min(newMean1))
        }
        for(i in seq_along(grp2)) {
          grp2[[i]][grp2[[i]] == 0] <- runif(nrow(grp2), min = 0.0001, max=min(newMean2))
        }
        for(i in 1:nrow(df1)) {
          x <- c(grp1[i,1:ncol(grp1)])
          y <- c(grp2[i,1:ncol(grp2)])
          x <- as.numeric(x)
          y <- as.numeric(y)
          if(n1 == n2){
            stats <- (t.test(x,y,paired = TRUE))
            pval[i] <- stats$p.value
          }else{
            stats <- (t.test(x,y,paired = FALSE))
            pval[i] <- stats$p.value
          }
        }
        pval <- as.data.frame(pval)
        pval[pval == "NaN"] <- 0
        tab <- cbind(log2FC,pval)
        newTab <- data.frame(Gene=tab$`df$Gene`,log2FC=tab$log2Rat,log10pvalue=-log10(tab$pval))
        newTab$log10pvalue[newTab$log10pvalue == "Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "-Inf"] <- 0
        newTab$log10pvalue[newTab$log10pvalue == "NaN"] <- 0
        #newTab$threshold = as.factor(abs(newTab$log2FC) >= 0.6 & newTab$log10pvalue >= 1.3)
        newTab1 <- merge(newTab, ann, by.x = "Gene", by.y = "Gene_Symbol",all.x = TRUE)
        newTab1[is.na(newTab1)] = "Unknown"
        attach(input$newTab1)
        newTab1$group = ifelse(newTab1$Division == input$level, "Selected Division","Other Division" )
        volc <- ggplot(newTab1,aes(x=log2FC,y=log10pvalue))+geom_point(aes(col=group))+scale_color_manual(values=c("black","green"),labels = c("Other divisions", "Selected Division"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
        volc+geom_text_repel(data=subset(newTab1,Division == input$level),aes(label=Gene),max.overlaps = 1000)
      }
    }
  })
  output$sum <- DT::renderDataTable({
    DT::datatable(summary(),options = list(lengthMenu = c(20, 50, 100), pageLength = 20))
  })
  output$Plot1 <- renderPlot({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    Plotting1()
  })
  output$Plot1_1 <- renderPlot({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    Plotting1_1()
  })
  output$Plot5 <- renderPlot({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    Plotting2()
  })
  output$vx <- renderUI({
    if(is.null(input$file)) {
      return() 
    }else{
      selectizeInput("Protein","Select the Protein for boxplot",choices = var(),options = list(maxOptions = 25000))
    }
  })
  box <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      if(input$imputate == "No"){
        grp1_G <- cbind(df$Gene,grp1)
        grp2_G <- cbind(df$Gene,grp2)
        tgrp1 <- t(grp1_G)
        tgrp2 <- t(grp2_G)
        attach(input$df)
        for(i in 1:ncol(tgrp2)) {
          if(tgrp1[1,i] == input$Protein){
            norm_inten <- tgrp1[2:nrow(tgrp1),i]
            norm1 <- as.data.frame(norm_inten)
            norm1$grp <- samN[[1]]
          }
        }
        for(i in 1:ncol(tgrp2)) {
          if(tgrp2[1,i] == input$Protein){
            norm_inten <- tgrp2[2:nrow(tgrp2),i]
            norm2 <- as.data.frame(norm_inten)
            norm2$grp <- samN[[2]]
          }
        }
        comb <- rbind(norm1,norm2)
        comb$norm_inten <- as.numeric(comb$norm_inten)
        print(comb)
        #print(input$Protein) 
        boxp <- ggboxplot(comb, x="grp", y="norm_inten",color="grp",add = "jitter",palette = "jco") + stat_compare_means(label = "p.signif", hide.ns = TRUE)+ggtitle(input$Protein)+ theme_bw() + theme(legend.position = "right",plot.title = element_text(hjust=0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("group")+ylab("Normalized intensity")
        boxp
      }else if(input$imputate == "Yes"){
        sub_main <- read.table("Data_imput.txt",sep = "\t",header = TRUE)
        df <- sub_main[2:nrow(sub_main),]
        dfGene <- as.data.frame(df$Gene)
        df <- arrange(df, Gene, desc(Gene))
        df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
        firstRow <- sub_main[1,1:ncol(sub_main)]
        df[is.na(df)] <- 0
        samN <- input$Group
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        grp1_G <- cbind(df$Gene,grp1)
        grp2_G <- cbind(df$Gene,grp2)
        tgrp1 <- t(grp1_G)
        tgrp2 <- t(grp2_G)
        attach(input$df)
        for(i in 1:ncol(tgrp2)) {
          if(tgrp1[1,i] == input$Protein){
            norm_inten <- tgrp1[2:nrow(tgrp1),i]
            norm1 <- as.data.frame(norm_inten)
            norm1$grp <- samN[[1]]
          }
        }
        for(i in 1:ncol(tgrp2)) {
          if(tgrp2[1,i] == input$Protein){
            norm_inten <- tgrp2[2:nrow(tgrp2),i]
            norm2 <- as.data.frame(norm_inten)
            norm2$grp <- samN[[2]]
          }
        }
        comb <- rbind(norm1,norm2)
        comb$norm_inten <- as.numeric(comb$norm_inten)
        print(comb)
        #print(input$Protein) 
        boxp <- ggboxplot(comb, x="grp", y="norm_inten",color="grp",add = "jitter",palette = "jco") + stat_compare_means(label = "p.signif",hide.ns = TRUE)+ggtitle(input$Protein)+ theme_bw() + theme(legend.position = "right",plot.title = element_text(hjust=0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("group")+ylab("Normalized intensity")
        boxp
        
      }
    }
    
  })
  output$Plot2 <- renderPlot({
    box()
  })
  pie1 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      print(length(samN))
      
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(df$Gene,mean1,mean2)
      tot1 <- sum(mean1)
      tot2 <- sum(mean2)
      df1 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
      df1[is.na(df1)] = "Unknown"
      grp1int <- aggregate(df1$mean1, by=list(Category=df1$Division), FUN=sum)
      grp2int <- aggregate(df1$mean2, by=list(Category=df1$Division), FUN=sum)
      grp1int$percent <- (grp1int$x/tot1)*100
      grp2int$percent <- (grp2int$x/tot2)*100
      df2 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
      grp1int1 <- aggregate(df2$mean1, by=list(Category=df2$Category), FUN=sum)
      grp2int1 <- aggregate(df2$mean2, by=list(Category=df2$Category), FUN=sum)
      matritot1 <- sum(df2$mean1)
      matritot2 <- sum(df2$mean2)
      grp1int1$percent <- (grp1int1$x/matritot1)*100
      grp2int1$percent <- (grp2int1$x/matritot2)*100
      grp1Tot = paste(samN[[1]],"total",sep = "_")
      grp2Tot = paste(samN[[2]],"total",sep = "_")
      grp1mat = paste(samN[[1]],"Matri",sep = "_")
      grp2mat = paste(samN[[2]],"Matri",sep = "_")
      
      if(length(samN) == 2){
        fig <- plot_ly(textinfo='percent',height = 1200) %>%config(toImageButtonOptions = list(width = NULL, height = NULL))%>%
          add_pie(data = grp1int, name = samN[[1]], labels = ~Category, values = ~x,domain = list(row = 0, column = 0))%>%
          #add_pie(data = grp1int1, name = grp1mat, labels = ~Category, values = ~x,domain = list(row = 1, column = 0))%>%
          add_pie(data = grp2int, name = samN[[2]], labels = ~Category, values = ~x,domain = list(row = 0, column = 1))%>%
          #add_pie(data = grp2int1, name = grp2mat, labels = ~Category, values = ~x,domain = list(row = 1, column = 1))%>%
          layout(title = "", showlegend = T,
                 grid=list(rows=2, columns=2),
                 xaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 yaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 annotations = list(x = c(.001, .95),
                                    y = c(.98, .98),
                                    text = c(samN[[1]],samN[[2]]),
                                    xref = "papper",
                                    yref = "papper",
                                    showarrow =F
                 ),margin=list(r=100, l=70, t=40, b=100))
        print(fig)
      }else{
        print("select two sample groups!!!")
      }
      
    }
  })
  
  pie2 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      print(length(samN))
      
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(df$Gene,mean1,mean2)
      tot1 <- sum(mean1)
      tot2 <- sum(mean2)
      df1 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
      df1[is.na(df1)] = "Unknown"
      grp1int <- aggregate(df1$mean1, by=list(Category=df1$Division), FUN=sum)
      grp2int <- aggregate(df1$mean2, by=list(Category=df1$Division), FUN=sum)
      grp1int$percent <- (grp1int$x/tot1)*100
      grp2int$percent <- (grp2int$x/tot2)*100
      df2 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
      grp1int1 <- aggregate(df2$mean1, by=list(Category=df2$Category), FUN=sum)
      grp2int1 <- aggregate(df2$mean2, by=list(Category=df2$Category), FUN=sum)
      matritot1 <- sum(df2$mean1)
      matritot2 <- sum(df2$mean2)
      grp1int1$percent <- (grp1int1$x/matritot1)*100
      grp2int1$percent <- (grp2int1$x/matritot2)*100
      print(grp1int1)
      grp1Tot = paste(samN[[1]],"total",sep = "_")
      grp2Tot = paste(samN[[2]],"total",sep = "_")
      grp1mat = paste(samN[[1]],"Matri",sep = "_")
      grp2mat = paste(samN[[2]],"Matri",sep = "_")
      
      if(length(samN) == 2){
        fig <- plot_ly(textinfo='percent',height = 1200) %>%config(toImageButtonOptions = list(width = NULL, height = NULL))%>%
          add_pie(data = grp1int1, name = samN[[1]], labels = ~Category, values = ~x,domain = list(row = 0, column = 0))%>%
          add_pie(data = grp2int1, name = samN[[2]], labels = ~Category, values = ~x,domain = list(row = 0, column = 1))%>%
          layout(title = "", showlegend = T,
                 grid=list(rows=2, columns=2),
                 xaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 yaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 annotations = list(x = c(.001, .95),
                                    y = c(.98, .98),
                                    text = c(samN[[1]],samN[[2]]),
                                    xref = "papper",
                                    yref = "papper",
                                    showarrow =F
                 ),margin=list(r=100, l=70, t=40, b=100))
        print(fig)
      }else{
        print("select two sample groups!!!")
      }
      
    }
  })
  pie3 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      df[is.na(df)] <- 0
      samN <- input$Group
      print(length(samN))
      
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      mean2 <- rowMeans(grp2)
      mean1 <- as.data.frame(mean1)
      mean2 <- as.data.frame(mean2)
      avg_val <- cbind(df$Gene,mean1,mean2)
      tot1 <- sum(mean1)
      tot2 <- sum(mean2)
      df1 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
      df1[is.na(df1)] = "Unknown"
      grp1int <- aggregate(df1$mean1, by=list(Category=df1$Secretome), FUN=sum)
      grp2int <- aggregate(df1$mean2, by=list(Category=df1$Secretome), FUN=sum)
      grp1int$percent <- (grp1int$x/tot1)*100
      grp2int$percent <- (grp2int$x/tot2)*100
      df2 <- merge(avg_val, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
      grp1int1 <- aggregate(df2$mean1, by=list(Category=df2$Category), FUN=sum)
      grp2int1 <- aggregate(df2$mean2, by=list(Category=df2$Category), FUN=sum)
      matritot1 <- sum(df2$mean1)
      matritot2 <- sum(df2$mean2)
      grp1int1$percent <- (grp1int1$x/matritot1)*100
      grp2int1$percent <- (grp2int1$x/matritot2)*100
      print(grp1int1)
      grp1Sec = paste(samN[[1]],"Secretome",sep = "_")
      grp2Sec = paste(samN[[2]],"Secretome",sep = "_")
      grp1mat = paste(samN[[1]],"Matri",sep = "_")
      grp2mat = paste(samN[[2]],"Matri",sep = "_")
      
      if(length(samN) == 2){
        fig <- plot_ly(textinfo='percent',height = 900) %>%config(toImageButtonOptions = list(width = NULL, height = NULL))%>%
          add_pie(data = grp1int, name = samN[[1]], labels = ~Category, values = ~x,domain = list(row = 0, column = 0))%>%
          add_pie(data = grp2int, name = samN[[2]], labels = ~Category, values = ~x,domain = list(row = 0, column = 1))%>%
          layout(title = "", showlegend = T,
                 grid=list(rows=2, columns=2),
                 xaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 yaxis = list(showgrid = TRUE, zeroline = TRUE, showticklabels = TRUE),
                 annotations = list(x = c(.001, .95),
                                    y = c(.98, .98),
                                    text = c(samN[[1]],samN[[2]]),
                                    xref = "papper",
                                    yref = "papper",
                                    showarrow =F
                 ),margin=list(r=100, l=70, t=40, b=100))
        print(fig)
      }else{
        print("select two sample groups!!!")
      }
      
    }
  })
  
  output$Plot3 <- renderPlotly({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    pie1()
  })
  
  output$Plot4 <- renderPlotly({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    pie2()
  })
  output$Plot5 <- renderPlotly({
    dat <- data.frame(x = numeric(0), y = numeric(0))
    withProgress(message = 'Making plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
    })
    plot(dat$x, dat$y)
    pie3()
  })
  pie_table1 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      sam1 <- as.data.frame(sam[[samN[[1]]]])
      colnames(sam1) <- "samp1"
      sam2 <- as.data.frame(sam[[samN[[2]]]])
      colnames(sam2) <- "samp2"
      sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
      sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
      colnames(sam1) <- "samp1"
      colnames(sam2) <- "samp2"
      n1 <- nrow(sam1)
      n2 <- nrow(sam2)
      grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
      grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
      val1 <- 0
      val2 <- 0
      for(i in 1:n1) {
        val1[i] <- sum(grp1[,i])
      }
      High1 <- max(val1)
      for(i in 1:n2) {
        val2[i] <- sum(grp2[,i])
      }
      High2 <- max(val2)
      
      for(i in 1:ncol(grp1)) {
        tot <- sum(grp1[,i])
        grp1[,i] <- (grp1[,i]*High1)/tot
      }
      mean1 <- rowMeans(grp1)
      
      for(i in 1:ncol(grp2)) {
        tot <- sum(grp2[,i])
        grp2[,i] <- (grp2[,i]*High2)/tot
      }
      grp1_G <- cbind(df$Gene,grp1)
      grp2_G <- cbind(df$Gene,grp2)
      df1 <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
      print(head(df1))
      #df1b <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
      df2 <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
      print(head(df2))
      #df2b <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
      df1[is.na(df1)] = "Unknown"
      df2[is.na(df2)] = "Unknown"
      grp1divsp <- split(df1,df1$Division)
      grp2divsp <- split(df2,df2$Division)
      #grp1catsp <- split(df1b,df1b$Category)
      #grp2catsp <- split(df2b,df2b$Category)
      name1div <- names(grp1divsp)
      name2div <- names(grp2divsp)
      #name1cat <- names(grp1catsp)
      #name2cat <- names(grp2catsp)
      pvald <- 0
      pvalc <- 0
      for(i in 1:length(name1div)){
        temp1 <- grp1divsp[[i]]
        temp2 <- grp2divsp[[i]]
        dat1 <- temp1[,2:(ncol(temp1)-4)]
        dat2 <- temp2[,2:(ncol(temp2)-4)]
        col1meanDiv <- colMeans(dat1)
        col2meanDiv <- colMeans(dat2)
        x <- c(as.numeric(col1meanDiv))
        y <- c(as.numeric(col2meanDiv))
        if(n1 == n2){
          statsd <- (t.test(x,y,paired = TRUE))
          pvald[i] <- statsd$p.value
        }else{
          statsd <- (t.test(x,y,paired = FALSE))
          pvald[i] <- statsd$p.value
        }
      }
      divp <- as.data.frame(cbind(pvald,name1div))
      #catp <- as.data.frame(cbind(pvalc,name1cat))
      colnames(divp) <- c("pvalue","Division")
      #colnames(catp) <- c("pvalue","Category")
      divp$pvalue <- round(as.numeric(divp$pvalue),digits = 4)
      divp$pvalue <- as.character(divp$pvalue)
      divp$pvalue[is.na(divp$pvalue)] <- 'NA'
      print(divp)
      #print(catp)
    }
  })
  pie_table2 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      if(length(samN) == 2){
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        grp1_G <- cbind(df$Gene,grp1)
        grp2_G <- cbind(df$Gene,grp2)
        #df1 <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        df1b <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
        #df2 <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        df2b <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
        #df1[is.na(df1)] = "Unk"
        #df2[is.na(df2)] = "Unk"
        #grp1divsp <- split(df1,df1$Division)
        #grp2divsp <- split(df2,df2$Division)
        grp1catsp <- split(df1b,df1b$Category)
        grp2catsp <- split(df2b,df2b$Category)
        #name1div <- names(grp1divsp)
        #name2div <- names(grp2divsp)
        name1cat <- names(grp1catsp)
        name2cat <- names(grp2catsp)
        pvald <- 0
        pvalc <- 0
        for(i in 1:length(name1cat)){
          temp1a <- grp1catsp[[i]]
          temp2a <- grp2catsp[[i]]
          dat1a <- temp1a[,2:(ncol(temp1a)-4)]
          dat2a <- temp2a[,2:(ncol(temp2a)-4)]
          col1meanCat <- colMeans(dat1a)
          col2meanCat <- colMeans(dat2a)
          x <- c(as.numeric(col1meanCat))
          y <- c(as.numeric(col2meanCat))
          if(n1 == n2){
            statsc <- (t.test(x,y,paired = TRUE))
            pvalc[i] <- statsc$p.value
            pvalc[is.na(i)] <- "NA"
          }else{
            statsc <- (t.test(x,y,paired = FALSE))
            pvalc[i] <- statsc$p.value
            pvalc[is.na(i)] <- "NA"
          }
        }
        #divp <- as.data.frame(cbind(pvald,name1div))
        catp <- as.data.frame(cbind(pvalc,name1cat))
        #colnames(divp) <- c("pvalue","Division")
        colnames(catp) <- c("pvalue","Category")
        catp$pvalue <- round(as.numeric(catp$pvalue),digits = 4)
        #print(divp)
        catp$pvalue <- as.character(catp$pvalue)
        catp$pvalue[is.na(catp$pvalue)] <- 'NA'
        print(catp)
      }else{
        print("select two sample groups!!!")
      }
    }
  })
  pie_table3 <- reactive({
    if(is.null(input$file))
      return()
    else 
    {
      nfiles = nrow(input$file) 
      inFile = list()
      if(see_if(has_extension(input$file$datapath, 'txt')) == TRUE){
        main <- read.table(input$file$datapath[1],sep = "\t",header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.table(input$file$datapath[1],sep = "\t",header = FALSE)
      }else if(see_if(has_extension(input$file$datapath, 'csv')) == TRUE){
        main <- read.csv(input$file$datapath[1],header = TRUE)
        #main <- inFile[[1]]
        for_sam <- read.csv(input$file$datapath[1],header = FALSE)
      }
      sub_main <- select(main,matches("_|Gene"))
      df <- sub_main[2:nrow(sub_main),]
      dfGene <- as.data.frame(df$Gene)
      df[is.na(df)] <- 0
      df <- arrange(df, Gene, desc(Gene))
      df[,2:ncol(df)] <- as.data.frame(sapply(df[,2:ncol(df)], as.numeric))
      firstRow <- sub_main[1,1:ncol(sub_main)]
      df[is.na(df)] <- 0
      All_sam <- for_sam[1:2,]
      sub_sam <- select(All_sam,matches("_|Gene"))
      tsam_sub <- t(All_sam)
      tsam_sub <- as.data.frame(tsam_sub)
      colnames(tsam_sub) <- c("sample","Group")
      match_sam <- dcast(tsam_sub,sample~Group,margins = TRUE)
      sub_match_sam <- match_sam[grep("_", match_sam$sample), ]
      sub_match_sam <- sub_match_sam[,1:ncol(sub_match_sam)-1]
      sub_match_sam <- sub_match_sam[, colSums(sub_match_sam != 0) > 0]
      write.table(sub_match_sam,"proc_samples.txt",sep = "\t",quote = FALSE,row.names = FALSE)
      arg1 <- "proc_samples.txt"
      arg2 <- "Samples.txt"
      cmd <- paste("perl", "process_samplefile.pl", arg1,arg2)
      system(cmd)
      sam_init <- read.table("Samples.txt",sep = "\t",header = TRUE)
      sam_init[sam_init == 0] = ""
      sam <- sam_init[,2:ncol(sam_init)]
      ann <- ""
      if(input$species == "Human"){
        ann <- read.table("Class_annot_v3.txt",sep = "\t",header = T)
      }else if(input$species == "Mouse"){
        ann <- read.table("Class_annot_m_v3.txt",sep = "\t",header = T)
      }
      samN <- input$Group
      if(length(samN) == 2){
        sam1 <- as.data.frame(sam[[samN[[1]]]])
        colnames(sam1) <- "samp1"
        sam2 <- as.data.frame(sam[[samN[[2]]]])
        colnames(sam2) <- "samp2"
        sam1 <- as.data.frame(sam1[!(is.na(sam1$samp1) | sam1$samp1 == ""), ])
        sam2 <- as.data.frame(sam2[!(is.na(sam2$samp2) | sam2$samp2 == ""), ])
        colnames(sam1) <- "samp1"
        colnames(sam2) <- "samp2"
        n1 <- nrow(sam1)
        n2 <- nrow(sam2)
        grp1 <- df %>% select(one_of(dput(as.character(sam1$samp1))))
        grp2 <- df %>% select(one_of(dput(as.character(sam2$samp2))))
        val1 <- 0
        val2 <- 0
        for(i in 1:n1) {
          val1[i] <- sum(grp1[,i])
        }
        High1 <- max(val1)
        for(i in 1:n2) {
          val2[i] <- sum(grp2[,i])
        }
        High2 <- max(val2)
        
        for(i in 1:ncol(grp1)) {
          tot <- sum(grp1[,i])
          grp1[,i] <- (grp1[,i]*High1)/tot
        }
        mean1 <- rowMeans(grp1)
        
        for(i in 1:ncol(grp2)) {
          tot <- sum(grp2[,i])
          grp2[,i] <- (grp2[,i]*High2)/tot
        }
        grp1_G <- cbind(df$Gene,grp1)
        grp2_G <- cbind(df$Gene,grp2)
        df1 <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        print(head(df1))
        #df1b <- merge(grp1_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
        df2 <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol",all.x = TRUE)
        print(head(df2))
        #df2b <- merge(grp2_G, ann, by.x = "df$Gene", by.y = "Gene_Symbol")
        df1[is.na(df1)] = "Unknown"
        df2[is.na(df2)] = "Unknown"
        grp1divsp <- split(df1,df1$Secretome)
        grp2divsp <- split(df2,df2$Secretome)
        #grp1catsp <- split(df1b,df1b$Category)
        #grp2catsp <- split(df2b,df2b$Category)
        name1div <- names(grp1divsp)
        name2div <- names(grp2divsp)
        #name1cat <- names(grp1catsp)
        #name2cat <- names(grp2catsp)
        pvald <- 0
        pvalc <- 0
        for(i in 1:length(name1div)){
          temp1 <- grp1divsp[[i]]
          temp2 <- grp2divsp[[i]]
          dat1 <- temp1[,2:(ncol(temp1)-4)]
          dat2 <- temp2[,2:(ncol(temp2)-4)]
          col1meanDiv <- colMeans(dat1)
          col2meanDiv <- colMeans(dat2)
          x <- c(as.numeric(col1meanDiv))
          y <- c(as.numeric(col2meanDiv))
          if(n1 == n2){
            statsd <- (t.test(x,y,paired = TRUE))
            pvald[i] <- statsd$p.value
          }else{
            statsd <- (t.test(x,y,paired = FALSE))
            pvald[i] <- statsd$p.value
          }
        }
        divp <- as.data.frame(cbind(pvald,name1div))
        #catp <- as.data.frame(cbind(pvalc,name1cat))
        colnames(divp) <- c("pvalue","Division")
        #colnames(catp) <- c("pvalue","Category")
        divp$pvalue <- round(as.numeric(divp$pvalue),digits = 4)
        divp$pvalue <- as.character(divp$pvalue)
        divp$pvalue[is.na(divp$pvalue)] <- 'NA'
        print(divp)
      }else{
        print("select two sample groups!!!")
      }
    }
  })
  output$PieData1 <- DT::renderDataTable({
    DT::datatable(pie_table1(),options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  },digits = 4)
  output$PieData2 <- DT::renderDataTable({
    DT::datatable(pie_table2(),options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  },digits = 4)
  output$PieData3 <- DT::renderDataTable({
    DT::datatable(pie_table3(),options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  },digits = 4)
  
  output$download1 <- downloadHandler(
    filename = function() {
      Nsam <- input$Group
      paste("Diff_Analysis", Sys.Date(),"_",Nsam[[2]],"_vs_",Nsam[[1]],".txt",sep = "")
    },
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      
      # Write to a file specified by the 'file' argument
      write.table(Analysis(), file, sep = "\t",quote = FALSE,row.names = FALSE)
    }
  )
  output$download2 <- downloadHandler(
    
    filename = function() {
      Nsam <- input$Group
      paste("VolcanoPlot", Sys.Date(),"_",Nsam[[2]],"_vs_",Nsam[[1]],".pdf",sep = "")
    },
    
    content = function(file) {
      # Write to a file specified by the 'file' argument
      pdf(file=file,width = 15,height = 10)
      print(Plotting1())
      dev.off()
    }
  )
  output$download3 <- downloadHandler(
    
    filename = function() {
      paste("BoxPlot", Sys.Date(),".pdf",sep = "")
    },
    
    content = function(file) {
      # Write to a file specified by the 'file' argument
      pdf(file=file)
      print(box())
      dev.off()
    }
  )
  output$download4 <- downloadHandler(
    
    filename = function() {
      paste("testData", Sys.Date(),".txt",sep = "")
    },
    
    content = function(file) {
      # Write to a file specified by the 'file' argument
      write.table(testData(), file, sep = "\t",quote = FALSE,row.names = FALSE)
    }
  )
  
  #output$download5 <- downloadHandler(
  
  #  filename = function() {
  #paste("testSamples", Sys.Date(),".txt",sep = "")
  #},
  
  #content = function(file) {
  # Write to a file specified by the 'file' argument
  #write.table(testSam(), file, sep = "\t",quote = FALSE,row.names = FALSE)
  #}
  #)
  output$download6 <- downloadHandler(
    
    filename = function() {
      paste("Imputed_data", Sys.Date(),".txt",sep = "")
    },
    
    content = function(file) {
      # Write to a file specified by the 'file' argument
      write.table(imputdownload(), file, sep = "\t",quote = FALSE,row.names = FALSE)
    }
  )
  output$download7 <- downloadHandler(
    filename = function() {
      Nsam <- input$Group
      paste("VolcanoPlot_wlables", Sys.Date(),"_",Nsam[[2]],"_vs_",Nsam[[1]],".pdf",sep = "")
    },
    
    content = function(file) {
      # Write to a file specified by the 'file' argument
      pdf(file=file,width = 15,height = 10)
      print(Plotting1_1())
      dev.off()
    }
  )
  output$tb <- renderUI({
    if(is.null(input$file)) {
      h6(tags$p("This app takes the protein intensities from two groups of samples selected and performs annotation and differential analysis"),
         tags$p("The sample data file should be in a spcific format(given below) all sample replicates and the corresponding samples should be in consecutive rows in the file"),
         tags$strong("sample data file"),
         br(),
         tags$img(src="Sample_data.png",width = "550px", height = "200px"),
         #tags$div("All the sample replicates column should have an underscore (_) ",style="color:blue"),
         tags$div("Columns colored RED are the mandatory Gene-Sample columns used for downstream processing, the column name for genes should be 'Gene' and for samples there should be an underscore(_) which should not be there in the non-sample column",style = "color:red"),
         tags$div("Download the testdata as a reference for data formatting and test run",style = "color:blue"),
         br(),
         downloadButton("download4", "Download test data"))
      
      #downloadButton("download5", "Download test Samples"),
    }
    else
      tabsetPanel(
        tabPanel("Files Uploaded", tableOutput("fileob"),downloadButton("download6", "Download Imputed Data" )),
        tabPanel("Overview", DT::dataTableOutput("sum")),
        tabPanel(("Pie Chart"),bsCollapse(open="",multiple = TRUE,bsCollapsePanel("Secretome",plotlyOutput("Plot5", height = "700px",width = "1000px"),(DT::dataTableOutput("PieData3",width = "500px")),style = "success"),bsCollapsePanel("Division",plotlyOutput("Plot3", height = "900px",width = "1000px"),(DT::dataTableOutput("PieData1",width = "500px")),style = "success"),bsCollapsePanel("Category",(plotlyOutput("Plot4", height = "900px",width = "1000px")),(DT::dataTableOutput("PieData2",width = "600px")),style = "success"))),
        #tabPanel("Output data",downloadButton("download1", "Download Data" ), DT::dataTableOutput("newdata1")),
        tabPanel(("Differential Analysis"),bsCollapse(open="",multiple = TRUE,bsCollapsePanel("Volcano Plot",downloadButton("download2", "Download Plot" ),plotOutput("Plot1",height = "700px",width = "1200px"),downloadButton("download1", "Download Data" ), DT::dataTableOutput("newdata1"),style = "success"),bsCollapsePanel("Volcano Plot with labels",downloadButton("download7", "Download Plot" ),plotOutput("Plot1_1",height = "700px",width = "1200px"),style = "success"))),
        #tabPanel("Volcano Plot",bsCollapse(open="Volcano 1",multiple = TRUE,bsCollapsePanel("Volcano 1",downloadButton("download2", "Download Plot" ),plotOutput("Plot1",height = "500px",width = "900px"),style = "success"),bsCollapsePanel("Volcano 2",downloadButton("download4", "Download Plot" ),plotOutput("Plot5",height = "500px",width = "900px"),style = "success"))),
        tabPanel("Box Plot",downloadButton("download3", "Download Plot"),plotOutput("Plot2"))
        #,column(11,offset = 11,(tableOutput("PieData1"))),column(11,offset = 11,(tableOutput("PieData2")))
      )
  })
}
shinyApp(ui, server)
#}