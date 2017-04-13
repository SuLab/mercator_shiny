library(shiny)
library(plotly)
library(scatterD3)


## jswidth <-
##     '$(document).on("shiny:conntected", function(e) {
##     var jsWidth = screen.width;
##     Shiny.onInputChange("GetScreenWidth",jsWidth)
## });
## '

## jsheight <-
##     '$(document).on("shiny:conntected", function(e) {
##     var jsHeight = screen.height;
##     Shiny.onInputChange("GetScreenHeight",jsHeight)
## });
## '

tsne.y <- readRDS('data/plotData.RDS')

meta.choices <- as.list(names(tsne.y)[3:length(tsne.y)])
names(meta.choices) <- names(tsne.y)[3:length(tsne.y)]
## meta.choices <- meta.choices[c('project','tissue_general','tissue_detail','tumor_id','tumor_stage','tumor_sample_type','living','age','sex','ethnicity','extraction_kit','seq_platform','sources')]


div(class="outer",    
    tags$head(includeCSS("styles.css")),
    ## tags$script(jswidth),
    ## tags$script(jsheight),
    plotlyOutput('tsne',width='100%',height='100%'),
    ## plotlyOutput('tsne',width='1400px',height='1400px'),
    ## plotlyOutput('tsne',width=jsWidth,height=jsHeight),
    absolutePanel(id = "controls", class = "panel panel-default", fixed=TRUE,
                  draggable=FALSE, top=60, left = 20, right="auto", bottom = "auto",
                  width=330, height="auto",

                  h2("Color Controls"),
    
                  ## selectInput('genevsgroup',
                  ##              label = 'Choose method of coloring',
                  ##              choices = list('No Coloring' = 1,'Metadata' = 2),
                  ##              selected = 1),

              
                  ## radioButtons('genevsgroup',
                  ##              label = 'Choose method of coloring',
                  ##              choices = list('No Coloring' = 1,'Metadata' = 2),
                  ##              selected = 1),

                  selectInput('colorfactors',
                              label = 'What to color by?',
                              choices = meta.choices,
                              multiple=TRUE),

                  fileInput('euclid_input',
                            label='Euclidian coloring',
                            accept = c(
                                "text/tsv",
                                "text/tab-separated-values,text/plain",
                                ".tsv")
                            ),

                  fileInput('spearman_input',
                            label='Spearman coloring',
                            accept = c(
                                'text/tsv',
                                'text/tab-separated-values,text/plain',
                                '.tsv')
                            )

                  ## checkboxGroupInput('colorfactors',
                  ##                    label = 'What to color by?',
                  ##                    choices = meta.choices,
                  ##                    selected = 'project')
                  )
    )


## shinyUI(fluidPage(

##     titlePanel("t-SNE test"),
##     fluidRow(

##         column(6,
##             radioButtons('genevsgroup',
##                          label = 'Choose method of coloring',
##                          choices = list('No Coloring' = 1,'Metadata' = 2),
##                          selected = 1)),

##         column(6,
##                checkboxGroupInput('colorfactors',
##                                label = 'What to color by?',
##                                choices = meta.choices,
##                                selected = 'project'))
##         ),

##     fluidRow(
##         column(12,
##                plotlyOutput('tsne',width='1400px',height='1400px')
##                )
##     )
   
## ))



## tsne.y <- read.table('data/tsne_points.tsv',sep='\t',header=T)

## gene.choices <- as.list(c('None',names(tsne.y)[5:length(tsne.y)]))
## names(gene.choices) <- c('None',names(tsne.y)[5:length(tsne.y)])

        ## column(4,
        ##     radioButtons('genevsgroup',
        ##                  label = 'Choose method of coloring',
        ##                  choices = list('No Coloring' = 1,'Metadata' = 2,'Gene' = 3),
        ##                  selected = 1)),

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

               ## plotOutput('tsne'))
               ## scatterD3Output('tsne',width='1400px',height='1400px')
