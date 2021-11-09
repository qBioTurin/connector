
UpdatingFCMall<-function(S.cl,CONNECTORList_now,input, output, session,rResult,Flags)
{
  IndexesInfo<- IndexesPlot.Extrapolation(stability.list =S.cl)
  TableGeneration<-function(df){
    Tight = 
    Tight_rep = lapply(1:length(df[, 1]), function(i) {
      freq <- df[i, "Freq"]
      do.call("rbind", lapply(1:freq, function(j) df[i,]))
    })
    df_rep <- do.call("rbind", Tight_rep)
    dt.fr.max <- aggregate(V~ClusterH, data = df_rep,
                           FUN = "max")
    colnames(dt.fr.max)<-c("Model","Max")
    dt.fr.min<- aggregate(V~ClusterH, data = df_rep,
                          FUN = "min")
    colnames(dt.fr.min)<-c("Model","Min")
    dt.fr.median <- aggregate(V~ClusterH, data = df_rep,
                              FUN = "median")
    colnames(dt.fr.median)<-c("Model","Median")
    dt.fr.mean <- aggregate(V~ClusterH, data = df_rep,
                            FUN = "mean")
    colnames(dt.fr.mean)<-c("Model","Mean")
    OutTable<-Reduce(merge, list(dt.fr.mean, dt.fr.median, dt.fr.min, dt.fr.max))
    
    rownames(OutTable) <- OutTable$Model
    return(OutTable)
  }
  
  output$table_tight_all<- renderTable({ 
    OutTable<-TableGeneration(IndexesInfo$IndexesValues$Tight)
    OutTable[,-1]
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb_all <- renderTable({ 
    OutTable<-TableGeneration(IndexesInfo$IndexesValues$fDB)
    OutTable[,-1]
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb1_all<- renderTable({ 
    OutTable<-TableGeneration(IndexesInfo$IndexesValues$fDB1)
    OutTable[,-1]
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  output$table_fdb2_all<- renderTable({ 
    OutTable<-TableGeneration(IndexesInfo$IndexesValues$fDB2)
    OutTable[,-1]
  },bordered = TRUE,rownames = TRUE, colnames = TRUE)
  
  rResult$IndexBoxPlot <- IndexesInfo
  rResult$CONNECTORListFCM_all <- S.cl
  output$visualIndexes <- renderPlot(IndexesInfo$Plot)
  
  ## Selecting a number of clusters:
  # observeEvent(input$GConsMat,{
  #   ConsMatrices <- ConsMatrix.Extrapolation(stability.list = S.cl,
  #                                             data = CONNECTORList_now)
  #   output$visualConsMatrix <- renderPlot({ConsMatrices[[paste0("G",input$GConsMat)]]$ConsensusPlot})
  #   rResult$ConsMatrices <- ConsMatrices
  # })
  
  observeEvent(input$GExtrap,{
    if(!is.null(rResult$CONNECTORList.FCM))
    {
      session$sendCustomMessage(type = 'testmessage', 
                                message = 'The number of clusters choise has been updated!') 
    }
    
    CONNECTORList.FCM <- MostProbableClustering.Extrapolation(stability.list = S.cl,
                                                                 G = input$Gvalue )
    rResult$CONNECTORList.FCM  <- CONNECTORList.FCM 
    rResult$G <- input$Gvalue
    rResult$h <- length(CONNECTORList.FCM$FCM$fit$parameters$Lambda[1,])
    Flags$ContinueCLplot = F
  })
  
  Flags$Continue <- FALSE
  
  output$IndexesClusterError <- renderText({
      validate(need(!is.null(rResult$CONNECTORList.FCM), "Please select a number of clusters!!!") )
 
    Flags$Continue = TRUE
    "Done!"
  })
  observeEvent(Flags$Continue, {
    if(Flags$Continue)
    {
      ConnFCM <- rResult$CONNECTORList.FCM
      output$table_fdb_erre <- renderTable({
        ConnFCM$Cl.Info$Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb_esse <- renderTable({
        ConnFCM$Cl.Info$Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb_emme <- renderTable({
        ConnFCM$Cl.Info$Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_erre <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_esse <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb1_emme <- renderTable({
        ConnFCM$Cl.Info$Deriv.Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_erre <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$erre
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_esse <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$esse
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
      output$table_fdb2_emme <- renderTable({
        ConnFCM$Cl.Info$Deriv2.Coefficents$emme
      },bordered = TRUE,rownames = TRUE, colnames = TRUE)
    }
  })
  
}

PlotModify <- function(plot,title,x,y,feature = NULL)
{
    if(is.null(feature)){
      plot = plot + labs( x = x, y= y, title = title)
    }else{
      plot = plot + labs( x = x, y= y, title = title)
    }
    
    return(plot)
}

withConsoleRedirect <- function(expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr )
  if (length(txt) > 0) {
  #   # insertUI(paste0("#", containerId), where = "beforeEnd",
  #   #          ui = paste0(txt, "\n", collapse = "") )
     return(txt)
  }
  return("")
}

Visual_1Step <- function(CONNECTORList,input, output, session,rResult,Flags,From=NULL)
{
  TimeGridPlot <- TimeGridDensity(CONNECTORList)
  TimeGridPlot <- TimeGridPlot$TimeGrid_plot
  #output$visualGrowthEND <- renderPlot({OutGrowPlot})
  #rResult$GrowPlot <- OutGrowPlot
  output$visualTimes <- renderPlot({TimeGridPlot})
  
  # updating the feature to select:
  Feat = colnames(CONNECTORList$LabCurv)
  updateSelectInput(session, "feature",
                    choices = Feat )
  updateSliderInput(session,"truncTime",
                    min = min(CONNECTORList$TimeGrid),
                    max = max(CONNECTORList$TimeGrid),
                    value = c(min(CONNECTORList$TimeGrid), max(CONNECTORList$TimeGrid) ) )
  
  ## Growth plot:
  #1) change the feature
  observeEvent(c(input$feature,input$truncTime,input$GrowthCurvesLegend),{
    ### preparing dataframe for ggplot
    dataplot <- data.frame(merge(CONNECTORList$Dataset,CONNECTORList$LabCurv,by="ID"))
    col <- as.character(unique(dataplot[,input$feature]) )
    colFetaure <- rainbow(dim(unique(CONNECTORList$Lab[input$feature]))[1])
    
    ### Set growth curve plot with ggplot
    GrowPlot <- ggplot(data=dataplot, aes(x=Time, y=Vol,group=ID,col=as.factor(dataplot[,input$feature]) )) +
      geom_line() +
      geom_point() +
      labs(col=input$feature)+
      #labs(title=title,x=axes.x, y = axes.y,col=input$feature)+
      theme(plot.title = element_text(hjust = 0.5),
            title =element_text(size=10, face='bold'))+
      scale_colour_manual(values = colFetaure,limits = col, breaks = sort(col), name = input$feature)
    
    ###############################
    #GrowPlot<- GrowthCurve(CONNECTORList, feature = "ID")
    
    ## Time Grid plot
    MinTrunc<-input$truncTime[1]
    MaxTrunc<-input$truncTime[2]
    
    GrowPlot <- GrowPlot+
      geom_vline(xintercept = input$truncTime) 
    
    if(input$GrowthCurvesLegend)
      GrowPlot <- GrowPlot + theme(legend.position = "right")
    else 
      GrowPlot <- GrowPlot + theme(legend.position = "none")
      
    TimeGridPlot <- TimeGridPlot +
      geom_vline(xintercept = input$truncTime)+
      geom_hline(yintercept = input$truncTime)
      
    output$visualTimes <- renderPlot({TimeGridPlot})
    output$visualGrowth <- renderPlot({GrowPlot})
    rResult$GrowPlot <- GrowPlot
  })
  
  #2) truncate the curves
  observeEvent(input$TruncateData,{
    MinTrunc<-input$truncTime[1]
    MaxTrunc<-input$truncTime[2]
    if(MaxTrunc != 0 &
       MinTrunc != 0 &
       (MaxTrunc != max(CONNECTORList$TimeGrid) || MinTrunc != min(CONNECTORList$TimeGrid)) )
    {
      txt <- withConsoleRedirect(
        CONNECTORList_trunc<- DataTruncation(CONNECTORList,
                                             feature= input$feature,
                                             truncTime = input$truncTime)
        )
      output$SummaryCutting <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
      
      rResult$CONNECTORList_trunc <- CONNECTORList_trunc

    }else{
      output$SummaryCutting <- renderUI({HTML(paste("Please select better cutting times!"))} ) 
    } 
  })
  
  ## Now it is possible to save something:
  if(is.null(From))
  {
    Flags$save1 = 0
    Flags$GoS1toS2 = 0
  }else{
    Flags$save1 = 0
    Flags$GoS1toS2_RDS = 0
  }
  
}

server <- function(input, output, session) {
  options(shiny.maxRequestSize=500*1024^2)
  
  CONNECTORList<-NULL
  CONNECTORList_trunc<-NULL
  
  # DECLARE REACTIVEVALUES FUNCTION HERE
  rResult <- reactiveValues(CONNECTORList_trunc = CONNECTORList_trunc,
                            CONNECTORList = CONNECTORList,
                            p = NULL,
                            h = NULL,
                            G = NULL,
                            GrowPlot=NULL,
                            TimeGrid=NULL,
                            CrossLL = NULL,
                            IndexBoxPlot = NULL,
                            FCMplots = NULL,
                            DiscriminantPlot = NULL,
                            DiscriminantPlotF = NULL,
                            ConsMatrices = NULL
                            )
  
  Flags  <- reactiveValues(save1 = 1, GoS1toS2 = 1,GoS1toS2_RDS=1, GoS3toS4 = 1)
## Data import:

  Flags$ContinueLoadingData <- FALSE
  Flags$ContinueEstimP<- FALSE
  Flags$ContinueEstimG<- FALSE
  Flags$ContinueUploadingRDS <- FALSE
  Flags$Discrplot <- FALSE
  Flags$ContinueCLplot <- FALSE
  
  observeEvent(input$LoadingData, {
      output$LoadingError1 <- renderText({
    validate(
      need(!is.null(input$fileGrowth) & !is.null(input$fileFeat) , 
           if(is.null(input$fileGrowth)) "Please select a xlsx file!!"
           else "Please select a csv/txt file!!" )
    )
    
    Flags$ContinueLoadingData <- TRUE
    ""
  })
    
    if(Flags$ContinueLoadingData)
    {
      txt <- withConsoleRedirect(
        CONNECTORList <- DataImport(input$fileGrowth$datapath,input$fileFeat$datapath)
      )
      rResult$CONNECTORList <- CONNECTORList
      output$SummaryLoading <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
      
      Visual_1Step(CONNECTORList,input, output, session,rResult,Flags)
      }
  })

## Run CLLik:
  observeEvent(input$Minp, {
    if(input$Minp >= input$Maxp) updateNumericInput(session,"Maxp",value = input$Minp+1, min = input$Minp+1 )
    updateSliderInput(session,"pValue",
                      min = input$Minp,
                      max = input$Maxp,
                      value = input$Minp)
  }) 
  observeEvent(input$RunForP, {
    output$ErrorInP <- renderText({
      validate(
        need(!is.null(rResult$CONNECTORList_trunc) || !is.null(rResult$CONNECTORList) , 
             "Please upload the RData with a ConnectorList object or run the preprocessing step!" )
      )
      Flags$ContinueEstimP<- TRUE
      ""
    })
        
    if(Flags$ContinueEstimP)
    { 
      if(!is.null(rResult$CONNECTORList_trunc)) CONNECTORList_now <- rResult$CONNECTORList_trunc
      else if(!is.null(rResult$CONNECTORList)) CONNECTORList_now <- rResult$CONNECTORList
      cat("#### The CrossLogLikelihood Method is running... \n")
      CrossLogLike<-BasisDimension.Choice(CONNECTORList_now,input$Minp:input$Maxp)
      cat("#### The CrossLogLikelihood Method has successfully finished \n")
      CrossLogLikePlot<-CrossLogLike$CrossLogLikePlot
      output$visualKnots<-renderPlot({CrossLogLike$KnotsPlot})
      
      rResult$CrossLL <- CrossLogLike
      
      observeEvent(input$pValue,{
        if(input$pSelection == input$Minp)
          CrossLogLikePlot_out <- CrossLogLikePlot
        else CrossLogLikePlot_out <- CrossLogLikePlot + geom_vline(xintercept = input$pValue)
        
        output$visualCLLik <- renderPlot({CrossLogLikePlot_out})
      })
    }
  })
  observeEvent(input$pSelection,{
    rResult$p <- input$pValue
    Flags$GoS3toS4 <- 0
    
    mess <- paste("The value of p selected is",input$pValue)
    output$ErrorInP <- renderText({mess})
  })
## Run FCM:
  observeEvent(c(input$MinG,input$MaxG), {
    if(input$MinG >= input$MaxG)
      updateNumericInput(session,"MaxG",value = input$MinG, min = input$MinG+1 )
  }) 
  observeEvent(input$RunForG, {
    output$FCMallError <- renderText({
      validate(need(!(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList)) ,
                    "Please upload the RData with a ConnectorList object or run the preprocessing step!") )
      validate(need(!(is.null(rResult$p)) ,
                    "Please select the value of p!") )
      
      Flags$ContinueEstimG <- T
      ""
      })
    
    if(Flags$ContinueEstimG)
    {
      if(!is.null(rResult$CONNECTORList_trunc)) CONNECTORList_now <- rResult$CONNECTORList_trunc
      else if(!is.null(rResult$CONNECTORList)) CONNECTORList_now <- rResult$CONNECTORList
      
      updateSliderInput(session, "GConsMat",
                        min =  input$MinG,max = input$MaxG)
      updateNumericInput(session, "Gvalue",
                         value =  input$MinG, min = input$MinG,max = input$MaxG)
      
      cat("#### The FCM Method is running... \n")
      S.cl <-ClusterAnalysis(CONNECTORList_now,
                             G=input$MinG:input$MaxG,
                             p=rResult$p,
                             runs=input$Nruns,
                             seed=2404,
                             Cores=input$Cores)
      cat("#### The FCM Method has successfully finished \n")
      
      UpdatingFCMall(S.cl,CONNECTORList_now,input, output, session,rResult,Flags)
    }
  })
## FCM visual
# Selecting a number of clusters:
  observeEvent(input$GConsMat,{
    validate(
      need(expr = !is.null(rResult$CONNECTORListFCM_all), "Run the FCM" )
      )
    
    if(!is.null(rResult$CONNECTORList_trunc)) {
      CONNECTORList_now <- rResult$CONNECTORList_trunc
    }else{
      CONNECTORList_now <- rResult$CONNECTORList
    }
    
    rResult$CONNECTORListFCM_all -> S.cl
    if(is.null(rResult$ConsMatrices)){
      ConsMatrices <- ConsMatrix.Extrapolation(stability.list = S.cl,
                                             data = CONNECTORList_now)
      rResult$ConsMatrices <- ConsMatrices
    }else{
      rResult$ConsMatrices -> ConsMatrices
    }
    
    output$visualConsMatrix <- renderPlot({ConsMatrices[[paste0("G",input$GConsMat)]]$ConsensusPlot})
    
  })
  
## Run CLuster Visualization
  observeEvent(rResult$CONNECTORList.FCM,{
    
    if(!is.null(rResult$CONNECTORList_trunc)) {
      CONNECTORList_now <- rResult$CONNECTORList_trunc
    }else{
      CONNECTORList_now <- rResult$CONNECTORList
    }
    
    output$IndexesClusterErrorCLplot <- renderText({
      validate(need(!is.null(rResult$CONNECTORList.FCM) ,
                    if((is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList)) ){
                      "Please upload the RData with a ConnectorList object or run the preprocessing step!"
                      }else{
                        "Please run the FCM and select a number of clusters!!!"
                        }
                    ) 
               )
      
            if(!Flags$ContinueCLplot){
                          FCMplots<- ClusterWithMeanCurve(clusterdata = rResult$CONNECTORList.FCM,
                                            data= CONNECTORList_now,
                                            feature = input$featureCLplot)
                          rResult$FCMplots <- FCMplots
                          Flags$ContinueCLplot = TRUE
            
            
            updateSelectInput(session, "SplineID",
                              choices = names(FCMplots$spline.plots),
                              selected = names(FCMplots$spline.plots)[1]
            )
            clNames <- rResult$CONNECTORList.FCM$FCM$cluster$cluster.names
            Feat <- names(CONNECTORList_now$LabCurv)
    
            updateSelectInput(session, "featureCLplot",
                              choices = Feat,
                              selected = "ID" )
            updateSelectInput(session, "FeatInsp",
                              choices = Feat,
                              selected = "ID" ) 
            updateSelectInput(session, "typeCLplot",
                              choices = c("All clusters","Cluster means",paste(clNames,"Cluster") ),
                              selected = "All clusters" )
    }
    ""
  })
})
  observeEvent(c(input$SplineID,Flags$ContinueCLplot),{
    rResult$FCMplots -> FCMplots
    output$visualSplinePlot <- renderPlot({FCMplots$spline.plots[[input$SplineID]]})
    rResult$visualSplinePlot <- FCMplots$spline.plots[[input$SplineID]]
  })
  observeEvent(c(input$featureCLplot,Flags$ContinueCLplot),{
    output$G_not_selected <- renderText({
      validate(need(!is.null(rResult$G) ,
                      "Please select a number of clusters!"
                    )
               )
      ""
      })
    if(Flags$ContinueCLplot)
    {
      observeEvent(c(input$typeCLplot,input$ClusterPlotLegend),{
        if(!is.null(rResult$FCMplots)){
          rResult$FCMplots -> FCMplots
          data = FCMplots$plotsCluster$ALL$plot_env$data
          curves <- data.frame(Times=data$Dataset$Time,
                               Vol=data$Dataset$Vol,
                               ID=data$Dataset$ID,
                               Cluster=FCMplots$plotsCluster$ALL$plot_env$classificate,
                               Info=rep(t(data$LabCurv[input$featureCLplot]),data$LenCurv))
          curves$Cluster <- FCMplots$plotsCluster$ALL$plot_env$symbols[curves$Cluster]
          col<- as.character(unique(curves$Info))
          colFetaure <- rainbow(dim(unique(data$Lab[input$featureCLplot]))[1])
          if(input$typeCLplot == "All clusters"){
            visualClplot <- ggplot()+
              scale_linetype_manual(values =1:FCMplots$plotsCluster$ALL$plot_env$G ,
                                    limits=sort(FCMplots$plotsCluster$ALL$plot_env$symbols),
                                    breaks=sort(FCMplots$plotsCluster$ALL$plot_env$symbols),name="Cluster") +
              facet_wrap(~Cluster)+
              geom_line(data = curves,
                        aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
              scale_colour_manual(values = colFetaure,
                                  limits=col,
                                  breaks=sort(col),
                                  name=input$featureCLplot)+
              theme(plot.title = element_text(hjust = 0.5),
                    axis.line = element_line(colour = "black"),
                    panel.background = element_blank())+
              geom_line(data=FCMplots$plotsCluster$ALL$plot_env$plot_data,
                        aes(x=Times,y=means,linetype= as.factor(Cluster)),size = 1.2 )+
              ylim(FCMplots$plotsCluster$ALL$plot_env$ymin,FCMplots$plotsCluster$ALL$plot_env$ymax)+
              xlim(FCMplots$plotsCluster$ALL$plot_env$xmin,FCMplots$plotsCluster$ALL$plot_env$xmax)+
              labs(title=paste("Other parameters: p = ", FCMplots$plotsCluster$ALL$plot_env$p, ", h = ", FCMplots$plotsCluster$ALL$plot_env$h, ", G = ", FCMplots$plotsCluster$ALL$plot_env$G  ,sep = ""),
                   x = FCMplots$plotsCluster$ALL$plot_env$axis.x,
                   y = FCMplots$plotsCluster$ALL$plot_env$axis.y)
            }else if(input$typeCLplot == "Cluster means"){
              visualClplot <- rResult$FCMplots$plotMeanCurve
              }else{
                gsub(" Cluster", "",input$typeCLplot )->symb
                visualClplot <- ggplot()+
                scale_linetype_manual(values =1:FCMplots$plotsCluster$ALL$plot_env$G ,
                                      limits=sort(FCMplots$plotsCluster$ALL$plot_env$symbols),
                                      breaks=sort(FCMplots$plotsCluster$ALL$plot_env$symbols),
                                      name="Cluster") +
                labs(title=paste("Cluster",symb),
                     x = FCMplots$plotsCluster$ALL$plot_env$axis.x,
                     y = FCMplots$plotsCluster$ALL$plot_env$axis.y)+
                geom_line(data = curves[curves$Cluster==symb,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
                scale_colour_manual(values = colFetaure,limits=col,breaks=sort(col),name=input$featureCLplot)+
                theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+
                geom_line(data=FCMplots$plotsCluster$ALL$plot_env$plot_data[FCMplots$plotsCluster$ALL$plot_env$plot_data$Cluster==symb,],
                          aes(x=Times,y=means,linetype= as.factor(Cluster)),size = 1.2 )+
                ylim(FCMplots$plotsCluster$ALL$plot_env$ymin,FCMplots$plotsCluster$ALL$plot_env$ymax)+
                xlim(FCMplots$plotsCluster$ALL$plot_env$xmin,FCMplots$plotsCluster$ALL$plot_env$xmax)
                }
          if(!input$ClusterPlotLegend){
            visualClplot <- visualClplot + theme(legend.position = "none")
            }
          output$visualClplot <- renderPlot({visualClplot})
          rResult$visualClplot<- visualClplot
          }
      })
      }
    })
  observeEvent(c(input$ColorDiscrPlot,input$FeatDiscrCol,Flags$ContinueCLplot),{
    if(Flags$ContinueCLplot)
    {
      if(!is.null(rResult$CONNECTORList_trunc)) {
        CONNECTORList_now <- rResult$CONNECTORList_trunc
      }else{
        CONNECTORList_now <- rResult$CONNECTORList
      }
      
      if(input$ColorDiscrPlot == "Features" &  Flags$Discrplot){
        output$FeatSelection <- renderUI({
          selectInput("FeatDiscrCol", "Features:",
                      choices =  names(CONNECTORList_now$LabCurv))
        })
        Flags$Discrplot = F
      }else if(input$ColorDiscrPlot == "Clusters"){
        Flags$Discrplot = T
      }
      
      output$DiscrPlot <- renderPlot({
        validate(need( rResult$h <= 2,   "h value must be lower or equal to 2."   ) )
        if(is.null(input$FeatDiscrCol)) 
          FeatDiscrCol = "ID"
        else
          FeatDiscrCol = input$FeatDiscrCol
          
        DiscrPlot <- DiscriminantPlot(clusterdata = rResult$CONNECTORList.FCM,
                                        data = CONNECTORList_now,
                                        feature = FeatDiscrCol)
        
        if(input$ColorDiscrPlot == "Clusters")
          visualDiscrPlot <- DiscrPlot$ColCluster
        else
          visualDiscrPlot <- DiscrPlot$ColFeature
          
        output$DiscrPlot <- renderPlot({visualDiscrPlot})
        rResult$visualDiscrPlot<- visualDiscrPlot
        rResult$DiscriminantPlotF <- DiscrPlot$ColFeature
        rResult$DiscriminantPlot <- DiscrPlot$ColCluster
      })
      
    }
  })
  observeEvent(c(input$FeatInsp,Flags$ContinueCLplot),{
    if(Flags$ContinueCLplot){
      
      if(!is.null(rResult$CONNECTORList_trunc)) {
        CONNECTORList_now <- rResult$CONNECTORList_trunc
      }else{
        CONNECTORList_now <- rResult$CONNECTORList
      }
      
      NumberSamples<-CountingSamples(clusterdata= rResult$CONNECTORList.FCM,
                                     data = CONNECTORList_now,
                                     feature = input$FeatInsp )
      output$tableCounting<- renderTable({ 
        NumberSamples$Counting
        },bordered = TRUE, colnames = TRUE)
      output$tableClusterNames<- renderTable({ 
        NumberSamples$ClusterNames
      },bordered = TRUE, colnames = TRUE) 
      
    }

  })
  
## Costume Plots
# Data Visual
  observeEvent(input$ChangeGrowth,{
    pl<-rResult$GrowPlot
    output$visualGrowthEND <- renderPlot( PlotModify(pl,input$titleGrowth,input$xlabGrowth,input$ylabGrowth) )
  }) 
  output$PlotDownloadGrowth <- downloadHandler(
    filename = function() { paste("GrowthCurve.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$GrowPlot , device = "pdf",width = 10,height = 10,units = "cm")
    }
  )
# Cross Log Likelihood
  observeEvent(input$ChangeCLL,{
    pl<-rResult$CrossLL
    output$visualCLLikEND <- renderPlot( PlotModify(pl,input$titleCLL,input$xlabCLL,input$ylabCLL) )
  }) 
  output$PlotDownloadCLL <- downloadHandler(
    filename = function() { paste("CrossLogLike.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$CrossLL , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cons Matrix
  output$PlotDownloadBoxPlotIndexes <- downloadHandler(
    filename = function() { paste("BoxPlotIndexes.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$BoxPlotIndexes , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cons Matrix
  output$PlotDownloadConsMatrix <- downloadHandler(
    filename = function() { paste("ConsMatrix.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$ConsMatrix , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# Cluster plot
  output$PlotDownloadClplot <- downloadHandler(
    filename = function() { paste("Clplot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualClplot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# DiscrPlot
  output$PlotDownloadDiscrPlot <- downloadHandler(
    filename = function() { paste("DiscrPlot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualDiscrPlot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  )
# spline plot
  output$PlotDownloadSplinePlot <- downloadHandler(
    filename = function() { paste("SplinePlot.pdf") },
    content = function(file) {
      ggsave(file, plot = rResult$visualSplinePlot , device = "pdf",width = 10,height = 8,units = "cm")
    }
  ) 
## Change panel
  observe({
    if(Flags$GoS1toS2 == 0){
      shinyjs::enable("FromS1toS2")
    }else{
      shinyjs::disable("FromS1toS2")
    }
  })
  observe({
    if(Flags$GoS1toS2_RDS == 0){
      shinyjs::enable("FromS1toS2_RDs")
    }else{
      shinyjs::disable("FromS1toS2_RDs")
    }
  })
  observe({
    if(Flags$GoS3toS4 == 0){
      shinyjs::enable("FromS3toS4")
    }else{
      shinyjs::disable("FromS3toS4")
    }
  })

  observeEvent(input$FromS1toS2, {
    updateTabsetPanel(session, "tabs",
                      selected = "PreProc"
    )
  })
  observeEvent(input$FromS1toS2_RDs, {
    updateTabsetPanel(session, "tabs",
                      selected = "PreProc"
    )
  })
  observeEvent(input$FromS2toS3, {
    updateTabsetPanel(session, "tabs",
                      selected = "pSelection"
    )
  })
  observeEvent(input$FromS3toS4, {
    updateTabsetPanel(session, "tabs",
                      selected = "FCM"
    )
  })
  
## Upload/ Save / Download:
  observeEvent(input$LoadingConnectorList1, {
    
    Flags$ContinueUploading <- FALSE
    
    output$LoadingError2 <- renderText({
      validate(
        need(!is.null(input$RDataImport1$datapath), "Please select a Connector.RData!!!")
      )
      Flags$ContinueUploading = TRUE
      ""
    })
    
    observeEvent(Flags$ContinueUploading, {

      if(Flags$ContinueUploading)
      {
        if(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList) )
        {
          load  = TRUE
          load(input$RDataImport1$datapath)
          
        }else{ ### alert!!! if it is aleready present! 
          showModal(modalDialog(
            title = "Important message",
            "Do you want to update the ConnectorList already present?",
            easyClose = TRUE,
            footer= tagList(actionButton("confirmUpload", "Update"),
                            modalButton("Cancel")
                            )
          )) 
          
        }
        
        observeEvent(input$confirmUpload & Flags$ContinueUploading, {
          removeModal()
          load(input$RDataImport1$datapath)
          
          mess = "Value uploaded:  "
          
          # Check if the list Info is present in the loaded file,
          # otherwise it is checked if the mandatory files are present
          # and in case the list Info is created
          
          
          if(!exists("Info")){
            validate(
              need(exists("CONNECTORList"),
                   "The CONNECTORList object is mandatory. If the object is present, please check its name which must correspond to CONNECTORList. "))
            
            Info = list()
            Info$CONNECTORList <- CONNECTORList
            rResult$CONNECTORList <- Info$CONNECTORList
             
            if(exists("CONNECTORList_trunc")) 
              Info$CONNECTORList_trunc <- CONNECTORList_trunc
            else Info$CONNECTORList_trunc <- NULL
            
            if(exists("CONNECTORListFCM_all")){
              Info$CONNECTORListFCM_all <- CONNECTORListFCM_all
              Info$p <- length(CONNECTORListFCM_all$Clusters.List[[1]]$ClusterAll[[1]]$FCM$fit$FullS[1,])
              }else{
                Info$CONNECTORListFCM_all <- NULL
                Info$p <- NULL
              }
            }
          
          if(is.null(Info$CONNECTORList_trunc)){
            rResult$CONNECTORList_trunc <- NULL
          } else {
            rResult$CONNECTORList_trunc <- Info$CONNECTORList_trunc
          }

          if(is.null(Info$p)){  
            mess <- c(mess,paste(" -No p value; "))
          }else{ 
            mess <- c(mess,paste(" -p = ",Info$p,"; "))
            rResult$p <- Info$p
            updateNumericInput(session,"pValueInHpanel",value = Info$p, min =1 )
            updateNumericInput(session,"pValueInGpanel",value = Info$p, min =1 )
          }

          if(is.null(Info$G)){
            mess <- c(mess,paste(" -No G value. "))
          }else{
            mess <- c(mess,paste(" -G = ",Info$G,". "))
            rResult$G <- Info$G
          }
          
          if(!is.null(Info$CONNECTORListFCM_all)){
            rResult$CONNECTORListFCM_all <- Info$CONNECTORListFCM_all
            UpdatingFCMall(S.cl = rResult$CONNECTORListFCM_all,
                           CONNECTORList_now = rResult$CONNECTORList,
                           input, output, session,rResult,Flags)
            
            Cl<-names(Info$CONNECTORListFCM_all$Clusters.List)
            Cl.number<-as.numeric(gsub("G", "", Cl))
            
            updateSliderInput(session, "GConsMat",
                              min =  min(Cl.number),max = max(Cl.number))
            updateNumericInput(session, "Gvalue",
                               value =  min(Cl.number), min =  min(Cl.number),max = max(Cl.number))
            updateNumericInput(session,"MaxG",value =  max(Cl.number) )
            updateNumericInput(session,"MinG",value =  min(Cl.number) )
          }
          if(!is.null(Info$CONNECTORList.FCM)){
            rResult$CONNECTORList.FCM <- Info$CONNECTORList.FCM
          }
          
          if(!is.null(Info$GrowPlot))
          {
            output$visualGrowth <- renderPlot({Info$GrowPlot})
            output$visualGrowthEND <- renderPlot({Info$GrowPlot})
          }else{
            Visual_1Step(rResult$CONNECTORList,input, output, session,rResult,Flags)
          }
          
          if(!is.null(Info$TimeGrid))
          {
            output$visualTimes <- renderPlot({Info$TimeGrid})
          }

          output$SummaryLoadingRdata <- renderUI({HTML(paste(mess, collapse = "<br/>"))} ) 
        })
        
      }
    })
    
    ## Now it is possible to save something:
    Flags$save1 = 0
  }) 
  observeEvent(input$LoadingConnectorList2, {
    ### è sbagliato, sto passando dei rdata, invece che dataframe!!!
    
    Flags$ContinueUploading <- FALSE
    
    output$LoadingError3 <- renderText({
      validate(
        need(!is.null(input$RDataImportGrowthData$datapath), "Please select an RDS file!!!")
      )
      Flags$ContinueUploading = TRUE
      ""
    })
    
    observeEvent(Flags$ContinueUploading, {
      
      if(Flags$ContinueUploading)
      {
        if(is.null(rResult$CONNECTORList_trunc) & is.null(rResult$CONNECTORList) )
        {
          Flags$ContinueUploadingRDS = TRUE
        }else{ ### alert!!! if it is aleready present! 
          showModal(modalDialog(
            title = "Important message",
            "Do you want to update the ConnectorList already present?",
            easyClose = TRUE,
            footer= tagList(actionButton("confirmUpload2", "Update"),
                            modalButton("Cancel")
            )
          ))
          
        }
        
        observeEvent(input$confirmUpload2, {
          Flags$ContinueUploadingRDS = TRUE
        })
        
        GrowthData <- readRDS(input$RDataImportGrowthData$datapath)
        cat(input$RDataImportGrowthData$datapath)
        if(Flags$ContinueUploadingRDS ){
    if(!is.null(input$RDataImportFeatureData$datapath)){
            FeatureData <- readRDS(input$RDataImportFeatureData$datapath)
            txt <- withConsoleRedirect(
              CONNECTORList <- DataFrameImport(GrowthData,FeatureData)
            )
          }else{
            txt <- withConsoleRedirect(
              CONNECTORList <- DataFrameImport(GrowthData)
            )
          }
          rResult$CONNECTORList <- CONNECTORList
          output$SummaryLoading <- renderUI({HTML(paste(txt[-1], collapse = "<br/>"))} ) 
          Visual_1Step(CONNECTORList,input, output, session,rResult,Flags,"FromRDS")
  }
        
      }
    })
    
    ## Now it is possible to save something:
    Flags$save1 = 0
  }) 
  observeEvent(input$saveRData1, {
    cat("save.image 1\n") 
    
    if( is.null(rResult$CONNECTORList) & is.null(rResult$CONNECTORList_trunc))
    {
      session$sendCustomMessage(type = 'testmessage',
                                message = 'There are no Connector objects to save!') 
    }else{
      CONNECTORList<-rResult$CONNECTORList
      CONNECTORList_trunc<-rResult$CONNECTORList_trunc
      Info <- reactiveValuesToList(rResult)
      
      save(Info, file = "Connector.RData")
      
      Flags$save1 = 0
    }
  })
  
  observe({
    if(Flags$save1 == 0){
      shinyjs::enable("download")
    }else{
      shinyjs::disable("download")
    }
  })
  
  output$download <- downloadHandler(
    filename = function(){
      paste("Connector.RData")
    },
    content = function(file) {
      #if(file_ext(file)!="RData") paste0(file_path_sans_ext(file),".RData")
      cat("file.copy\n")
      file.copy(from = "Connector.RData", to = file)
    },
    contentType = NULL
  )
  
## Report
  output$reportGen <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      report.path <- system.file("Shiny","report.Rmd", package = "connector")
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy(report.path, tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      infoReport <-  reactiveValuesToList(rResult)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = list(infoReport = infoReport),
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
## end
}


