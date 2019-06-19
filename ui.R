library(shiny)
library(plotly)
library(scatterD3)
library(shinyTree)
library(DT)
library(shinyjs)

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

## tsne.y <- readRDS('data/all_recount_metasra_summarized.RDS')

## meta.choices <- as.list(names(tsne.y)[5:length(tsne.y)])
## names(meta.choices) <- names(tsne.y)[5:length(tsne.y)]

## meta.choices <- readRDS('data/metadata_fields.RDS')
meta.choices <- c('No Coloring' = 'No Coloring',
                  'Project ID' = 'proj_id',
                  'GTEx/TCGA Gross Tissue' = 'tissue_general',
                  'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                  'Sample Type' = 'sample_type',
                  'Mesh: Anatomy' = 'Mesh',
                  'Tissue' = 'Tissue',
                  'DOID' = 'DOID',
                  'EFO: Cultured Cells'='efo',
                  'Gene' = 'Gene',
                  ## 'Kmeans' = 'KMeans',
                  'Louvain' = 'Louvain',
                  'Projection' = 'Projection')



## tsne.choices <- readRDS('data/tsne_pca_list_names.RDS')
## tsne.choices <- as.list(tsne.choices.vec)
## tsne.choices <- names(readRDS('data/recount_tsne_pca_list_noNAs_over50_noSingle_protein_coding.RDS'))

## names(tsne.choices) <- sapply(tsne.choices,function(x) {
##     y <- strsplit(x,'[.]')[[1]]
##     sprintf('Perplexity = %s, PCs = %s',y[1],y[2])
## })
## tsne.choices <- as.list(tsne.choices)

louvain.vec <- readRDS('data/leiden_pc3sd_r25e-3_vec.RDS')
## louvain.vec <- readRDS('data/louvain_vec_pca_over50_noSingle_k100_entrez.RDS')
naked.louvain.choices <- sort(unique(louvain.vec))
louvain.choices <- sapply(naked.louvain.choices,function(x) sprintf('Louvain Cluster %s',x))
names(louvain.choices) <- louvain.choices

## names(tsne.choices) <- tsne.choices
## meta.choices <- meta.choices[c('project','tissue_general','tissue_detail','tumor_id','tumor_stage','tumor_sample_type','living','age','sex','ethnicity','extraction_kit','seq_platform','sources')]


fluidPage(
    includeCSS('styles.css'),
    includeScript('jssrc/panel_click.js'),
    navbarPage('Mercator',id='topNavBar',
               tabPanel('t-SNE',value='tsnePanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left=0,right='auto',bottom='auto',width='100%',height='95%',value='tsnePanel',
                                      plotlyOutput('tsne',width='100%',height='100%')
                                      )
                        ),
               tabPanel('Violin',value='violinPanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=70,left=25,right='auto',bottom='auto',width='97%',height='90%',value='violinPanel',
                                      plotOutput('violin',width='100%',height='100%')
                                      )
                        ),
               ## tabPanel('violin',
               ##          absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left='10%',right='auto',bottom='auto',width='90%',height='87%',
               ##                        plotlyOutput('violin',width='100%',height='100%')
               ##                        )

               ##          ),
               tabPanel('Bar Plot',value='barPanel',
                        absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=75,left=15,right='auto',bottom='auto',width='99%',height='85%',value='barPanel',
                                      plotOutput('metadataBar',width='100%',height='100%')
                                      )
                        ),

               tabPanel('Gene Table',value='genePanel',
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'Cluster # for marker table')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=620,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('geneGroup',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  choices = c('All'='all',naked.louvain.choices),
                                                  multiple=FALSE,
                                                  selected='All')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=62,left=975,right='auto',bottom='auto',width='auto',height='auto',
                                      downloadButton('geneTableDownload',
                                                     label = 'Download gene table'
                                                     )
                                      ),
                        absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=50,right='auto',bottom='auto',width='99%',height='auto',value='genePanel',

                                    ## DT::DTOutput('geneTable',width='300px')
                                    DT::DTOutput('geneTable',width='90%')
                                    )
                        ),

               tabPanel('Sample Table',value='samplePanel',
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=72,left=435,right='auto',bottom='auto',width='auto',height='auto',
                                      tags$label(class='control-label', 'Subset data')
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=42,left=550,right='auto',bottom='auto',width='auto',height='auto',
                                      selectInput('sampleTableControl',
                                                  label='',
                                                  ## label='Cluster # for marker table',
                                                  choices = c('All','Selection'),
                                                  multiple=FALSE,
                                                  selected=0)
                                      ),
                        absolutePanel(fixed=FALSE,draggable=FALSE,top=62,left=875,right='auto',bottom='auto',width='auto',height='auto',
                                      downloadButton('sampleTableDownload',
                                                     label = 'Download sample table'
                                                     )
                                      ),
                        absolutePanel(id='results',fixed=FALSE,draggable=FALSE,top=120,left=50,right='auto',bottom='auto',width='99%',height='auto',
                                      DT::DTOutput('sampleTable',width='90%')
                                      )
                        
                        )
               ## ),
               ),
    
    absolutePanel(id='controlPanel', fixed=FALSE, draggable=TRUE, top=60, left=20, right='auto', bottom='auto', width=400, height='auto',
                  
                  ## wellPanel(

                  tags$div(class='panel panel-primary',
                           tags$div(class='panel-heading',id='controlPanelHeading',
                                    
                                    useShinyjs(),                           
                                    h3("Plot Controls",class='panel-title'),
                                    tags$span(class='pull-right clickable',id='controlPanelClickable',tags$i(class='glyphicon glyphicon-chevron-up'))
                                    ),
                           tags$div(class='panel-body',

                                    actionButton('colorButton','Redraw Plot'),
                                    
                                    tags$br(),
                                    tags$br(),

                                    sliderInput('markerSize','Adjust Marker Size',min=1,max=12,value=4),

                                    ## selectInput('whichTSNE',
                                    ##             label = 'What t-SNE to plot?',
                                    ##             choices = tsne.choices,
                                    ##             multiple=FALSE,
                                    ##             selected=tsne.choices[10]),

                                    ## tags$br(),
                                    ## tags$br(),

                                    uiOutput('plotlyClick'),
                                    
                                    tags$br(),

                                    selectizeInput('whichGene',
                                                   label='Select a gene',
                                                   choices=NULL),


                                    fileInput('gene.vec','Color by euclidian distance to sample',
                                              accept = c(
                                                  'text/tsv',
                                                  'text/tab-separated-values',
                                                  'text/plain',
                                                  '.tsv')
                                              ),

                                    DT::DTOutput('sampleInputTable',width='300px'),


                                    ## selectInput('whichGene',
                                    ##             label='Select a gene',
                                    ##             choices = gene.choices,
                                    ##             multiple=FALSE,
                                    ##             selectize=FALSE,
                                    ##             size=10),

                                    tags$br(),

                                    selectInput('colorfactors',
                                                label = 'What to color by?',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    ## selectInput('kmeansChoice',
                                    ##             label='Choose k for kmeans',
                                    ##             ## choices=1:100,
                                    ##             choices=200:400,
                                    ##             multiple=FALSE
                                    ##             ),

                                    selectInput('violinXFactors',
                                                label= 'Violin X-axis',
                                                choices = meta.choices,
                                                multiple=FALSE
                                                ),

                                    selectInput('violinYFactors',
                                                label='Violin Y-axis',
                                                choices = c('Projection'='projection','Gene'='gene'),
                                                multiple=FALSE,
                                                selected=c('Gene'='gene')
                                                ),

                                    selectInput('barPlotFactor',
                                                label='Bar plot group',
                                                ## choices = c('All','Selections','Kmeans Cluster 1',sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))),
                                                choices = c('All','Selections',louvain.choices),
                                                multiple=FALSE
                                                ),

                                    selectInput('barPlotXaxis',
                                                label='X axis for Bar Plot',
                                                choices=c('No Coloring' = 'No Coloring',
                                                          'Project ID' = 'proj_id',
                                                          'GTEx/TCGA Gross Tissue' = 'tissue_general',
                                                          'GTEx/TCGA Detailed Tissue' = 'tissue_detail',
                                                          'Sample Type' = 'sample_type',
                                                          'Mesh: Anatomy' = 'Mesh',
                                                          'Tissue' = 'Tissue',
                                                          'DOID' = 'DOID',
                                                          'EFO: Cultured Cells'='efo',
                                                          'Louvain' = 'Louvain'),
                                                multiple=FALSE
                                                ),

                                    shinyTree('meshTree',theme='proton'),

                                    shinyTree('tissueTree',theme='proton'),

                                    shinyTree('doidTree',theme='proton'),

                                    shinyTree('efoTree',theme='proton'),
                                    

                                    ## textInput('markerSearch','Search Marker Table'),
                                    
                                    textInput('selectionName','Selection Name',value=''),

                                    actionButton('saveSelection','Save Selection'),

                                    DT::DTOutput('selectionList',width='17%')

                                    ## selectInput('geneGroup',
                                    ##             label='Cluster # for marker table',
                                    ##             choices = naked.louvain.choices,
                                    ##             multiple=FALSE,
                                    ##             selected=0),

                                    ## DT::DTOutput('geneTable',width='300px')
                                    
                                    
                                    )
                           )
                  )
)

## fluidPage(
##     tags$div(class='outer',
##              titlePanel('Mercator'),
##              actionButton('hideButton','',icon=icon('angle-double-up')),

##              tabsetPanel(
##                  tabPanel('t-SNE',
##                           absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='95%',
##                                         plotlyOutput('tsne',width='100%',height='100%')
##                                         )
##                           ),
##                  tabPanel('violin',
##                           absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='85%',
##                                         plotOutput('violin',width='100%',height='100%')
##                                         )

##                           ),
##                  tabPanel('barPanel',
##                           absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='85%',
##                                         plotOutput('metadataBar',width='100%',height='100%')
##                                         )
##                           )
##              ),


##              absolutePanel(id='controlPanel', draggable=FALSE, top=80, left = 20, right= 'auto', bottom = 'auto', width=330, height='auto',
                           
##                            useShinyjs(),
##                            ## tags$div(id='controlPanel'),

##                            h3("Color Controls"),

##                            actionButton('colorButton','Redraw Plot'),
                           
##                            ## selectInput('whichTSNE',
##                            ##             label = 'What t-SNE to plot?',
##                            ##             choices = tsne.choices,
##                            ##             multiple=FALSE,
##                            ##             selected=tsne.choices[10]),

##                            selectizeInput('whichGene',
##                                           label='Select a gene',
##                                           choices=NULL),

##                            fileInput('gene.vec','Color by euclidian distance to sample',
##                                      accept = c(
##                                          'text/tsv',
##                                          'text/tab-separated-values',
##                                          'text/plain',
##                                          '.tsv')
##                                      ),

##                            ## selectInput('whichGene',
##                            ##             label='Select a gene',
##                            ##             choices = gene.choices,
##                            ##             multiple=FALSE,
##                            ##             selectize=FALSE,
##                            ##             size=10),

##                            selectInput('colorfactors',
##                                        label = 'What to color by?',
##                                        choices = meta.choices,
##                                        multiple=FALSE
##                                        ),

##                            selectInput('kmeansChoice',
##                                        label='Choose k for kmeans',
##                                        choices=1:100,
##                                        multiple=FALSE
##                                        ),

##                            selectInput('violinXFactors',
##                                        label= 'Violin X-axis',
##                                        choices = meta.choices,
##                                        multiple=FALSE
##                                        ),

##                            selectInput('violinYFactors',
##                                        label='Violin Y-axis',
##                                        choices = c('Projection'='projection','Gene'='gene'),
##                                        multiple=FALSE,
##                                        selected=c('Gene'='gene')
##                                        ),

##                            selectInput('barPlotFactor',
##                                        label='Bar plot group',
##                                        ## choices = c('All','Selections','Kmeans Cluster 1',sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))),
##                                        choices = c('All','Selections','Kmeans Cluster 1',louvain.choices),
##                                        multiple=FALSE
##                                        ),

##                            selectInput('barPlotXaxis',
##                                        label='X axis for Bar Plot',
##                                        choices=meta.choices,
##                                        multiple=FALSE
##                                        ),

##                            shinyTree('tree',theme='proton'),
                           
##                            selectInput('markerGroup',
##                                        label='KMeans Cluster for markers',
##                                        choices = 1:30,
##                                        multiple=FALSE,
##                                        selected=1),

##                            DT::DTOutput('markerTable',width='17%'),

##                            textInput('selectionName','Selection Name',value=''),

##                            actionButton('saveSelection','Save Selection'),

##                            DT::DTOutput('selectionList',width='17%')
##                            )
##              )        
## )
    
## fluidPage(
##     titlePanel('Mercator'),
##     actionButton('hideButton','',icon=icon('angle-double-up')),
##     fluidRow(
##         tags$div(id='controlPanel',
##                  column(2,

##                         useShinyjs(),
##                         ## tags$div(id='controlPanel'),

##                         h3("Color Controls"),

##                         actionButton('colorButton','Redraw Plot'),
                        
##                         ## selectInput('whichTSNE',
##                         ##             label = 'What t-SNE to plot?',
##                         ##             choices = tsne.choices,
##                         ##             multiple=FALSE,
##                         ##             selected=tsne.choices[10]),

##                         selectizeInput('whichGene',
##                                        label='Select a gene',
##                                        choices=NULL),

##                         fileInput('gene.vec','Color by euclidian distance to sample',
##                                   accept = c(
##                                       'text/tsv',
##                                       'text/tab-separated-values',
##                                       'text/plain',
##                                       '.tsv')
##                                   ),

##                         ## selectInput('whichGene',
##                         ##             label='Select a gene',
##                         ##             choices = gene.choices,
##                         ##             multiple=FALSE,
##                         ##             selectize=FALSE,
##                         ##             size=10),

##                         selectInput('colorfactors',
##                                     label = 'What to color by?',
##                                     choices = meta.choices,
##                                     multiple=FALSE
##                                     ),

##                         selectInput('kmeansChoice',
##                                     label='Choose k for kmeans',
##                                     choices=1:100,
##                                     multiple=FALSE
##                                     ),

##                         selectInput('violinXFactors',
##                                     label= 'Violin X-axis',
##                                     choices = meta.choices,
##                                     multiple=FALSE
##                                     ),

##                         selectInput('violinYFactors',
##                                     label='Violin Y-axis',
##                                     choices = c('Projection'='projection','Gene'='gene'),
##                                     multiple=FALSE,
##                                     selected=c('Gene'='gene')
##                                     ),

##                         selectInput('barPlotFactor',
##                                     label='Bar plot group',
##                                     ## choices = c('All','Selections','Kmeans Cluster 1',sapply(louvain.choices,function(x) sprintf('Louvain Cluster %s',x))),
##                                     choices = c('All','Selections','Kmeans Cluster 1',louvain.choices),
##                                     multiple=FALSE
##                                     ),

##                         selectInput('barPlotXaxis',
##                                     label='X axis for Bar Plot',
##                                     choices=meta.choices,
##                                     multiple=FALSE,
##                                     ),

##                         shinyTree('tree',theme='proton'),
                        
##                         selectInput('markerGroup',
##                                     label='KMeans Cluster for markers',
##                                     choices = 1:30,
##                                     multiple=FALSE,
##                                     selected=1),

##                         DT::DTOutput('markerTable',width='17%'),

##                         textInput('selectionName','Selection Name',value=''),

##                         actionButton('saveSelection','Save Selection'),

##                         DT::DTOutput('selectionList',width='17%')
##                         )
##                  ),
##         column(10,
##                tabsetPanel(
##                    tabPanel('t-SNE',
##                             absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='95%',
##                                           plotlyOutput('tsne',width='100%',height='100%')
##                                           )
##                             ),
##                    tabPanel('violin',
##                             absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='85%',
##                                           plotOutput('violin',width='100%',height='100%')
##                                           )

##                             ),
##                    tabPanel('barPanel',
##                             absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=150,left='15%',right='auto',bottom='auto',width='85%',height='85%',
##                                           plotOutput('metadataBar',width='100%',height='100%')
##                                           )
##                             )
##                )
               
##                )
##     )
## )



                          

## navbarPage('Mercator',
##            tabPanel('t-SNE',
##                     tags$head(includeCSS("styles.css")),
##                     fluidPage(
##                         fluidRow(

##                             ## absolutePanel(id = "controls", class = "panel panel-default", fixed=FALSE, draggable=FALSE, top=60, left = 20, right="auto", bottom = "auto", width='10%', height="auto",
##                             column(2,
                                   
##                                    h2("Color Controls"),                   
##                                    actionButton('colorButton','Redraw Plot'),
                                   
##                                    selectInput('whichTSNE',
##                                                label = 'What t-SNE to plot?',
##                                                choices = tsne.choices,
##                                                multiple=FALSE,
##                                                selected=tsne.choices[10]),

##                                    selectizeInput('whichGene',
##                                                   label='Select a gene',
##                                                   choices=NULL),

##                                    ## selectInput('whichGene',
##                                    ##             label='Select a gene',
##                                    ##             choices = gene.choices,
##                                    ##             multiple=FALSE,
##                                    ##             selectize=FALSE,
##                                    ##             size=10),
                                   

##                                    selectInput('colorfactors',
##                                                label = 'What to color by?',
##                                                choices = meta.choices,
##                                                multiple=FALSE
##                                                ),

##                                    selectInput('kmeansChoice',
##                                                label='Choose k for kmeans',
##                                                choices=1:100,
##                                                multiple=FALSE
##                                                ),

##                                    ## selectInput('violinXFactors',
##                                    ##             label= 'Violin X-axis',
##                                    ##             choices = meta.choices,
##                                    ##             multiple=FALSE
##                                    ##             ),

##                                    shinyTree('tree',theme='proton')

                                   
##                                    ),
##                             column(10,
##                                    absolutePanel(id = 'results',fixed=TRUE,draggable=FALSE,top=55,left='17%',right='auto',bottom='auto', width='83%',height='100%',
##                                                  plotlyOutput('tsne',width='100%',height='100%')
##                                                  )
##                                    )
                            
##                         )
##                     )
##                     ),
##            tabPanel('Violin',
##                     tags$head(includeCSS('styles.css')),
##                     fluidPage(
##                         fluidRow(
##                             column(2,

##                                    h2('Color Controls'),
##                                    actionButton('violinButton','Redraw Plot'),

##                                    selectInput('violinGene',
##                                                label='Select a gene',
##                                                choices=NULL),

##                                    selectInput('violinXFactors',
##                                                label= 'Violin X-axis',
##                                                choices = meta.choices,
##                                                multiple=FALSE
##                                                )
##                                    ),
##                             column(10,
                                   
##                                    absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=55,left='17%',right='auto',bottom='auto',width='83%',height='100%',
##                                                  plotOutput('violin',width='100%',height='90%')
##                                                  )
##                                    )
##                         )
##                         )
##                     ## absolutePanel(id = 'controls', class = 'panel panel-default', fixed=TRUE,
##                     ##               draggable=FALSE, top = 60, left = 20, right = 'auto', bottom = 'auto', width=330, height='auto',

##                     )
                                                                                                                                   
##            )





## fluidPage(
##     tags$head(includeCSS("styles.css")),

##     fluidRow(
##         column(3,
##                wellPanel(
##                    h2("Color Controls"),                   
##                    actionButton('colorButton','Redraw Plot'),
                   
##                    selectInput('whichTSNE',
##                                label = 'What t-SNE to plot?',
##                                choices = tsne.choices,
##                                multiple=FALSE),

##                    selectizeInput('whichGene',
##                                   label='Select a gene',
##                                   choices=NULL),

##                    ## selectInput('whichGene',
##                    ##             label='Select a gene',
##                    ##             choices = gene.choices,
##                    ##             multiple=FALSE,
##                    ##             selectize=FALSE,
##                    ##             size=10),
                   

##                    selectInput('colorfactors',
##                                label = 'What to color by?',
##                                choices = meta.choices,
##                                multiple=FALSE
##                                ),

##                    selectInput('violinXFactors',
##                                label= 'Violin X-axis',
##                                choices = meta.choices,
##                                multiple=FALSE
##                                ),

##                    shinyTree('tree',theme='proton')
##                )
##                ),
##         column(9,
##                tabsetPanel(
##                    tabPanel('t-SNE',plotlyOutput('tsne',height='100%',width='100%')),
##                    tabPanel('violin',plotOutput('violin')),
##                    tabPanel('bar plot',plotOutput('metadataBar'))
##                )
##                )
        
##     )
## )
    

    ## sidebarLayout(

    ##     sidebarPanel(
    ##         h2("Color Controls"),

    ##         actionButton('colorButton','Redraw Plot'),

    ##               selectInput('whichTSNE',
    ##                           label = 'What t-SNE to plot?',
    ##                           choices = tsne.choices,
    ##                           multiple=FALSE),
                           

    ##               selectInput('colorfactors',
    ##                           label = 'What to color by?',
    ##                           choices = meta.choices,
    ##                           multiple=FALSE
    ##                           ),

    ##         shinyTree('tree',theme='proton')
    ##     ),

    ##     mainPanel(
    ##         tabsetPanel(
    ##             tabPanel('t-SNE',plotlyOutput('tsne',width='100%',height='100%'))
    ##         )
    ##     )
    ## )


## div(class="outer",    
##     tags$head(includeCSS("styles.css")),
##     ## tags$script(jswidth),
##     ## tags$script(jsheight),
##     ## plotlyOutput('tsne',width='1400px',height='1400px'),
##     ## plotlyOutput('tsne',width=jsWidth,height=jsHeight),
##     absolutePanel(id = 'results',fixed=TRUE,draggable=FALSE,top=0,left=0,right='auto',bottom='auto', width='100%',height='100%',
                  
##                   tabsetPanel(type='tabs',
##                               tabPanel('t-SNE',plotlyOutput('tsne',width='100%',height='100%'))
##                               ## tabPanel('t-SNE 2',plotOutput('tsne',width='100%',height='100%'))
##                               )## ,
##                   ## plotlyOutput('tsne',width='100%',height='100%')
##                 ),
    
##     ## absolutePanel(id = 'tsne', fixed=TRUE, top=0, left=0, right='auto',bottom='auto',width='100%',height='100%',plotlyOutput('tsne',width='100%',height='100%')),

##     absolutePanel(id = "controls", class = "panel panel-default", fixed=TRUE,
##                   draggable=FALSE, top=60, left = 20, right="auto", bottom = "auto",
##                   width=330, height="auto",

##                   h2("Color Controls"),

##                   actionButton('colorButton','Redraw Plot'),
    
##                   ## selectInput('genevsgroup',
##                   ##              label = 'Choose method of coloring',
##                   ##              choices = list('No Coloring' = 1,'Metadata' = 2),
##                   ##              selected = 1),

              
##                   ## radioButtons('genevsgroup',
##                   ##              label = 'Choose method of coloring',
##                   ##              choices = list('No Coloring' = 1,'Metadata' = 2),
##                   ##              selected = 1),

##                   selectInput('whichTSNE',
##                               label = 'What t-SNE to plot?',
##                               choices = tsne.choices,
##                               multiple=FALSE),
                           

##                   selectInput('colorfactors',
##                               label = 'What to color by?',
##                               choices = meta.choices,
##                               multiple=FALSE,
##                               ),

##                   shinyTree('tree',theme='proton')

##                   ## fileInput('euclid_input',
##                   ##           label='Euclidian coloring',
##                   ##           accept = c(
##                   ##               "text/tsv",
##                   ##               "text/tab-separated-values,text/plain",
##                   ##               ".tsv")
##                   ##           ),

##                   ## fileInput('spearman_input',
##                   ##           label='Spearman coloring',
##                   ##           accept = c(
##                   ##               'text/tsv',
##                   ##               'text/tab-separated-values,text/plain',
##                   ##               '.tsv')
##                   ##           ),

##                   ## fileInput('euclid_pca_input',
##                   ##           label='Euclid PCA coloring',
##                   ##           accept = c(
##                   ##               'text/tsv',
##                   ##               'text/tab-separated-values,text/plain',
##                   ##               '.tsv')
##                   ##           )


##                   ## checkboxGroupInput('colorfactors',
##                   ##                    label = 'What to color by?',
##                   ##                    choices = meta.choices,
##                   ##                    selected = 'project')

##                   )

##     )


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

