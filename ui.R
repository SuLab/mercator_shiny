library(shiny)
library(plotly)
library(scatterD3)

## tsne.y <- read.table('data/tsne_points.tsv',sep='\t',header=T)

tsne.y <- readRDS('data/plotData.RDS')

meta.choices <- as.list(names(tsne.y)[3:length(tsne.y)])
names(meta.choices) <- names(tsne.y)[3:length(tsne.y)]
meta.choices <- meta.choices[c('project','tissue_general','tissue_detail','tumor_id','tumor_stage','tumor_sample_type','living','age','sex','ethnicity','extraction_kit','seq_platform','sources')]

## gene.choices <- as.list(c('None',names(tsne.y)[5:length(tsne.y)]))
## names(gene.choices) <- c('None',names(tsne.y)[5:length(tsne.y)])


shinyUI(fluidPage(

    titlePanel("t-SNE test"),
    fluidRow(

        ## column(4,
        ##     radioButtons('genevsgroup',
        ##                  label = 'Choose method of coloring',
        ##                  choices = list('No Coloring' = 1,'Metadata' = 2,'Gene' = 3),
        ##                  selected = 1)),
        column(6,
            radioButtons('genevsgroup',
                         label = 'Choose method of coloring',
                         choices = list('No Coloring' = 1,'Metadata' = 2),
                         selected = 1)),

        column(6,
               checkboxGroupInput('colorfactors',
                               label = 'What to color by?',
                               choices = meta.choices,
                               selected = 'project'))
        ),

        ## column(4,
        ##        selectInput("genealpha",
        ##                    label = 'Gene to alpha',
        ##                    choices = gene.choices,
        ##                    selected='None'))
               

        ## column(4,
        ##        fileInput('gene.vec','Color by euclidian distance to sample',
        ##                  accept = c(
        ##                      'text/tsv',
        ##                      'text/tab-separated-values',
        ##                      'text/plain',
        ##                      '.tsv')
        ##                  )
        ##        )
        ## ),
    
    ## fluidRow(

    ##     column(6,
    ##            selectInput('genecolor',
    ##                        label = 'Gene to color by',
    ##                        choices = gene.choices[2:length(gene.choices)])),

    ##     column(6,
    ##            radioButtons('groupShape',
    ##                               label = 'Change shapes based off',
    ##                               choices = list('Tissue' = 'tissue',
    ##                                              'K-Means Group' = 'kmeans.cluster'),
    ##                               selected = 1))
    ##     ),
    
    fluidRow(
        column(12,
               ## plotOutput('tsne'))
               ## scatterD3Output('tsne',width='1400px',height='1400px')
               plotlyOutput('tsne')
               )
    )
   
))
