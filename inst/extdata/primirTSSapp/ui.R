library(shiny)

fluidPage(
  titlePanel("pri-miRNA TSS"),
  navbarPage(
      "",
      tabPanel(
        "Find pri-miRNA TSS",
        sidebarLayout(
          sidebarPanel(
            h3("Find pri-miRNA TSS"),
            br(),

            radioButtons("file1_test", "have histone peaks file",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "yes"),
            fileInput("file1", "choose histone peaks file",
                      multiple = TRUE,
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")),
            numericInput(inputId = "n1",
                         label = "merge threshold of histone peaks",
                         value = 250),

            radioButtons("file2_test", "have pol2 peaks files",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "no"),
            fileInput("file2", "choose pol2 peaks files",
                      multiple = TRUE,
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")),
            numericInput(inputId = "n2",
                         label = "merge threshold of pol2 peaks",
                         value = 250),

            radioButtons("file3_test", "have miRNA expression profile",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "no"),
            fileInput("file3", "choose expressed miRNAs",
                      multiple = TRUE,
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")),

            radioButtons("file4_test", "have DHS file",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "no"),
            fileInput("file4", "choose DHS files",
                      multiple = TRUE,
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")),
            radioButtons("DHS_check_result", "allmirdhs_byforce",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "yes"),

            radioButtons("file5_test", "have gene expression profile",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "no"),
            fileInput("file5", "choose expressed genes",
                      multiple = TRUE,
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")),
            radioButtons("gene_check_result", "allmirgene_byforce",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "yes"),

            numericInput(inputId = "flanking_num",
                         label = "eponine flanking_num",
                         value = 1000),

            numericInput(inputId = "threshold",
                         label = "eponine threshold",
                         value = 0.7),

            radioButtons("tf", "get TFs simultaneously",
                         c("YES" = "yes",
                           "NO" = "no"),
                         selected = "no"),
            numericInput(inputId = "tf_n",
                         label = "seek n bps upstream from TSSs",
                         value = 1000),
            numericInput(inputId = "min.score",
                         label = "threshold for scoring transcription factor binding sites",
                         value = 0.8),

            actionButton("wm_start", "Start the analysis"),

            downloadButton('tss_result.csv', 'download the result')


          ), # sidebarPanel

          mainPanel(

            h3("View the first six rows of the result"),

            tableOutput("wm_exprs_mx")
          ) # mainPanel

        ), # sidebarLayout

      br(),
      br()
    ),

    tabPanel(
      "Plot pri-miRNA TSS",
      sidebarLayout(
        sidebarPanel(
          h3("plot pri-miRNA TSS"),
          br(),

          radioButtons("file1.2_test", "have histone peaks file",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "yes"),
          fileInput("file1.2", "choose histone peaks files",
                    multiple = TRUE,
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
          numericInput(inputId = "n1.2",
                       label = "merge threshold of histone peaks",
                       value = 250),

          radioButtons("file2.2_test", "have pol2 peaks files",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "no"),
          fileInput("file2.2", "choose pol2 peaks files",
                    multiple = TRUE,
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
          numericInput(inputId = "n2.2",
                       label = "merge threshold of pol2 peaks",
                       value = 250),

          textInput("file3.2", "expressed miRNAs",
                    value = "hsa-mir-5697"),

          radioButtons("file4.2_test", "have DHS file",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "no"),
          fileInput("file4.2", "choose DHS files",
                    multiple = TRUE,
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
          radioButtons("DHS_check_result.2", "allmirdhs_byforce",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "yes"),

          radioButtons("file5.2_test", "have gene expression profile",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "no"),
          fileInput("file5.2", "choose expressed genes",
                    multiple = TRUE,
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
          radioButtons("gene_check_result.2", "allmirgene_byforce",
                       c("YES" = "yes",
                         "NO" = "no"),
                       selected = "yes"),

          numericInput(inputId = "flanking_num.2",
                       label = "eponine flanking_num",
                       value = 1000),

          numericInput(inputId = "threshold.2",
                       label = "eponine threshold",
                       value = 0.7),

          actionButton("wm_start.2", "Start to plot")

        ), # sidebarPanel

        mainPanel(

          h3("View the plot"),

          plotOutput("wm_plot", width = 1100, height = 700)
        ) # mainPanel

      ), # sidebarLayout

      br(),
      br()
    )


  )
) # fluidPage
