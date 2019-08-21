library(shiny)
library(dplyr)
library(Seurat)
library(purrr)
library(cowplot)
library(parallel)
library(roxygen2)
library(reshape2)
library(DT)
library(tibble)
source('../DE_analysis.R')
source('../cell_cytometry.R')


seurat.out = readRDS('./input/WNIKINI45POS_out.rds')

# Define UI ----
ui <- fluidPage(titlePanel("Cell Feature Inspection"),
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      "n.cluster",
                      label = h3("The Cluster to be inspected in Feature Analysis"),
                      choices = unique(seurat.out$stim.markers$cluster),
                      selected = 0
                    ),
                    textInput(
                      "feature",
                      label = h3("Feature Inspected"),
                      value = "Enter text..."
                    ),

                    hr(),

                    selectInput(
                      "n.cluster.diff",
                      label = h3("The Cluster to be inspected in Differential Analysis"),
                      choices = unique(seurat.out$stim.markers$cluster),
                      selected = 0
                    ),
                    textInput(
                      "feature.diff",
                      label = h3("Differential Feature Inspected"),
                      value = "Enter text..."
                    ),



                    plotOutput("dimplot"),
                    plotOutput("percent.plot")

                  ),

                  ###################################################################################
                  mainPanel(
                    tabsetPanel(
                      id = 'feature genes',
                      tabPanel("stim.genes", DT::dataTableOutput("stimtab")),
                      tabPanel("ctrl.genes", DT::dataTableOutput("ctrltab")),
                      tabPanel("conserved.genes", DT::dataTableOutput("conservedtab"))
                    ),

                    plotOutput("feature.plot"),
                    hr(),

                    ###########################################################################
                    tabsetPanel(id = 'differential genes',
                                tabPanel(
                                  "diff.genes", DT::dataTableOutput("difftab")
                                )),

                    plotOutput("violin.plot")

                  )
                ))

# Define server logic ----
server <- function(input, output) {
  # You can access the value of the widget with input$select, e.g.
  output$dimplot <- renderPlot({
    DimPlot(
      seurat.out$ob,
      reduction = "umap",
      split.by = "stim",
      label = T
    )
  })

  output$stimtab <- DT::renderDataTable({
    #seurat.out$stim.markers[seurat.out$stim.markers$cluster == input$n.cluster,][1:input$num.gene,]
    seurat.out$stim.markers[seurat.out$stim.markers$cluster == input$n.cluster, ]
  })

  output$ctrltab <- DT::renderDataTable({
    #seurat.out$ctrl.markers[seurat.out$ctrl.markers$cluster == input$n.cluster,][1:input$num.gene,]
    seurat.out$ctrl.markers[seurat.out$ctrl.markers$cluster == input$n.cluster, ]
  })

  output$conservedtab <- DT::renderDataTable({
    #seurat.out$ctrl.markers[seurat.out$ctrl.markers$cluster == input$n.cluster,][1:input$num.gene,]
    seurat.out$conserved.markers[seurat.out$conserved.markers$cluster == input$n.cluster, ]
  })

  output$difftab <- DT::renderDataTable({
    #seurat.out$ctrl.markers[seurat.out$ctrl.markers$cluster == input$n.cluster,][1:input$num.gene,]
    seurat.out$diff[seurat.out$diff$cluster == input$n.cluster.diff, ]
  })

  output$feature.plot <- renderPlot({
    FeaturePlot(
      seurat.out$ob,
      features = input$feature,
      split.by = "stim",
      max.cutoff = 3,
      cols = c("grey", "red")
    )
  })

  output$percent.plot = renderPlot({
    p = Plot.Cluster.Percentage(seurat.out$ob, grouped = T)
    p
  })

  output$violin.plot = renderPlot({
    plots <-
      VlnPlot(
        seurat.out$ob,
        features = input$feature.diff,
        split.by = "stim",
        pt.size = 0,
        combine = FALSE
      )
    p = CombinePlots(plots = plots, ncol = 1)
    p
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
