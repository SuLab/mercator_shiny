library(shiny)
library(plotly)
library(scatterD3)
library(shinyTree)
library(DT)

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

meta.choices <- readRDS('data/metadata_fields.RDS')
meta.choices <- c('No Coloring',meta.choices,'Mesh','Gene','KMeans','Louvain')

tsne.choices <- readRDS('data/tsne_pca_list_names.RDS')
## tsne.choices <- as.list(tsne.choices.vec)
names(tsne.choices) <- sapply(tsne.choices,function(x) {
    y <- strsplit(x,'[.]')[[1]]
    sprintf('Perplexity = %s, PCs = %s',y[1],y[2])
})
tsne.choices <- as.list(tsne.choices)


## names(tsne.choices) <- tsne.choices
## meta.choices <- meta.choices[c('project','tissue_general','tissue_detail','tumor_id','tumor_stage','tumor_sample_type','living','age','sex','ethnicity','extraction_kit','seq_platform','sources')]

fluidPage(
    titlePanel('Mercator'),
               fluidRow(
                   column(2,
                          h3("Color Controls"),                   
                          actionButton('colorButton','Redraw Plot'),
                          
                          selectInput('whichTSNE',
                                      label = 'What t-SNE to plot?',
                                      choices = tsne.choices,
                                      multiple=FALSE,
                                      selected=tsne.choices[10]),

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

                          ## selectInput('whichGene',
                          ##             label='Select a gene',
                          ##             choices = gene.choices,
                          ##             multiple=FALSE,
                          ##             selectize=FALSE,
                          ##             size=10),

                          selectInput('colorfactors',
                                      label = 'What to color by?',
                                      choices = meta.choices,
                                      multiple=FALSE
                                      ),

                          selectInput('kmeansChoice',
                                      label='Choose k for kmeans',
                                      choices=1:100,
                                      multiple=FALSE
                                      ),

                          selectInput('violinXFactors',
                                      label= 'Violin X-axis',
                                      choices = meta.choices,
                                      multiple=FALSE
                                      ),

                          shinyTree('tree',theme='proton'),

                          textInput('selectionName','Selection Name',value=''),

                          actionButton('saveSelection','Save Selection'),

                          DT::DTOutput('selectionList',width='17%')
                          ),
                   column(10,
                          tabsetPanel(
                              tabPanel('t-SNE',
                                       absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=105,left='15%',right='auto',bottom='auto',width='85%',height='95%',
                                                     plotlyOutput('tsne',width='100%',height='100%')
                                                     )
                                       ),
                              tabPanel('violin',
                                       absolutePanel(id='results',fixed=TRUE,draggable=FALSE,top=105,left='15%',right='auto',bottom='auto',width='85%',height='85%',
                                                     plotOutput('violin',width='100%',height='100%')
                                                     )

                                       )
                          )
                          )
               )
    
)

                          

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

