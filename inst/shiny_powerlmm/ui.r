
library(shiny)
library(shinydashboard)

dashboardPage(

     dashboardHeader(title = "powerlmm"),
     dashboardSidebar(
          sidebarMenu(
               menuItem("Study design", icon = icon("dashboard"), tabName = "dashboard"),
               menuItem("Power curves", icon = icon("line-chart"), tabName = "pcurve")
          )
     ),
     dashboardBody(
         tags$head(
             tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
         ),
          tabItems(
               tabItem("dashboard",
                    fluidRow(

                              box(title = "Power analysis for longitudinal linear mixed models with missing data",
                                  width = 12,
                                  solidHeader=TRUE,
                                        p("A flexible implementation for power calculations in longitudinal linear mixed models, with or without clustering, random slopes, missing data, and unbalanced designs"),

                                  box(title = "Basic setup", width = 12, status = "warning",
                                      column(6,
                                             radioButtons("multi_des", "Is 2 or 3-level design?:",
                                                          c("Yes, 3-level" = "3lvl",
                                                            "No, 2-level" = "2lvl"),
                                                          inline=TRUE),
                                             conditionalPanel("input.multi_des == '3lvl'",
                                                              radioButtons("partially_nested", "Nesting/clutering in both treatment arms?:",
                                                                           c("Both arms" = FALSE,
                                                                             "Only treatment (partially-nested)" = TRUE),
                                                                           inline=TRUE)),
                                             radioButtons("use_standardized", "Do you want to use standardized or unstandardized inputs?",
                                                          c("Standardized" = "yes",
                                                            "Unstandardized" = "no"),
                                                          inline=TRUE),
                                             actionButton("button1", "Update", icon("refresh"),
                                                          class = "btn btn-primary")
                                             ),
                                      column(6,
                                             imageOutput("img_study_design")

                                             )


                                  )
                                  )

                       ),
                       fluidRow(
                           box("Sample size", width=3,
                                        helpText("Setup sample size at each level"),
                                        numericInput("n1", "n1 (measurements)", 11),

                                        textInput("n2", "n2 (subjects per cluster)", 5),
                                        helpText("Input a comma separated list for varying cluster sizes"),
                                        conditionalPanel("input.multi_des == '3lvl'",
                                            numericInput("n3", "n3 (clusters)", 4)),
                                        numericInput("T_end", "T_end", 10),
                                        h2("Missing data"),
                               h4("Control"),
                                    fluidRow(

                                        column(6, numericInput("retention_cc", "Total dropout", 0.2)),
                                        column(6, numericInput("gamma_cc", "Shape", 1/2))),
                               h4("Treatment"),
                                    fluidRow(
                                       column(6, numericInput("retention_tx", "Total dropout", 0.3)),
                                       column(6, numericInput("gamma_tx", "Shape", 2))),
                                        actionButton("button2", "Update", icon("refresh"),
                                                     class = "btn btn-primary")
                              ),
                           box(h4("Number of subjects"),
                               column(6, tableOutput('table_n2'))
                              ),
                           div('data-spy'="affix", 'data-offset-top'="660", 'class'='col-sm-3',
                               box(width=12,
                                   valueBoxOutput("power_box", width=12),
                                   valueBoxOutput("total_n", width = 12)
                               )
                           ),


                           box(
                               tabBox(width=12,
                                      tabPanel("Plot",  plotOutput("plot_retention", height = 250)
                                      ),
                                      tabPanel("Table",
                                               tableOutput('table_retention')
                                      ))


                       )),
                        fluidRow(
                            div("class" = "h-section",
                                h2("Treatment effects")),
                           box("Fixed effects", width=3,
                                    numericInput("fixed.intercept", "Intercept", 0),
                                    numericInput("fixed.slope", "Slope control", 1),
                                    numericInput("cohend", "Treatment effect (Cohen's d)", 0.5),
                               actionButton("button3", "Update", icon("refresh"),
                                            class = "btn btn-primary")),

                           box(plotOutput("plot_ES", height = 350))
                           ),
                            fluidRow(
                                div("class" = "h-section",
                                    h2("Variance components")),
                                column(3,
                                       box("Random effects", width = NULL,
                                           h4("Level 1"),
                                           conditionalPanel("input.use_standardized == 'no'",
                                                            numericInput("sigma.error", "Error (within subjects)", 2.6)),
                                           h4("Level 2"),
                                           conditionalPanel("input.use_standardized == 'no'",
                                                            numericInput("sigma.person.intercept", "Intercept subjects", 2.9)),
                                           conditionalPanel("input.use_standardized == 'yes'",
                                                            sliderInput("icc_pre_persons", "Proportion (%) baseline variance at level 2", min = 0, max= 100, value =50)),
                                           conditionalPanel("input.use_standardized == 'no'",
                                                            numericInput("sigma.person.slope", "Slope subjects", 0)),
                                           sliderInput("cor.person", "Correlation intercept-slope", min = -1, max = 1, step = 0.01, value = -0.5),
                                           conditionalPanel("input.multi_des == '3lvl'",
                                                            h4("Level 3")),
                                           conditionalPanel("input.use_standardized == 'no' && input.multi_des == '3lvl'",
                                                            numericInput("sigma.cluster.intercept", "Intercept clusters", 0)),
                                           conditionalPanel("input.use_standardized == 'yes' && input.multi_des == '3lvl'",
                                                            sliderInput("icc_pre_clusters", "Proportion (%) baseline variance at level 3", min = 0, max = 100, value = 0)),
                                           conditionalPanel("input.use_standardized == 'yes' && input.multi_des == '3lvl'",
                                                            sliderInput("icc_slope", "Proportion (%) slope variance at level 3", min = 0, max = 100, value = 5)),
                                           conditionalPanel("input.use_standardized == 'yes'",
                                                            numericInput("var_ratio", "Variance ratio", 0.03)),
                                           conditionalPanel("input.use_standardized == 'no'  && input.multi_des == '3lvl'",
                                                            numericInput("sigma.cluster.slope", "Cluster slope", 0)),
                                           conditionalPanel("input.multi_des == '3lvl'",
                                                            sliderInput("cor.cluster", "Correlation intercept-slope", min = -1, max = 1, step = 0.01, value = 0)
                                           ),
                                           actionButton("button4", "Update", icon("refresh"),
                                                        class = "btn btn-primary"))),




                      column(6,
                             box(width = NULL,
                                 tabBox(width=12,
                                        tabPanel("Plot",  plotOutput("plot_sds", height = 250)
                                        ),
                                        tabPanel("Table",
                                                 h4("Standard deviations per time point"),
                                                 tableOutput('table_sds')
                                        ))
                             ),

                             box(width = NULL,
                                 h4("Percentage variance per level"),
                                  p("The proportation of variance at the cluster level is also the intraclass correlation, i.e. the correlation between two subjects belonging to the same cluster. With random slopes these variances will be a quadratic function of the time variable"),
                                  tabBox(width=12,
                                         tabPanel("Plot", plotOutput("plot_VPC", height = 250)
                                         ),
                                         tabPanel("Table",
                                                  tableOutput('table_VPC')
                                         ))

                             ),
                             box(width = NULL,
                                 h4("Correlation between time points"),
                                 tabBox(width=12,
                                        tabPanel("Plot",  plotOutput('plot_correlation_matrix')
                                        ),
                                        tabPanel("Table",
                                                 tableOutput('table_correlation_matrix')
                                        ))
                             ))









               ),
               fluidRow(
                   box(width=12,
                       verbatimTextOutput("debug"))
               )),
               tabItem("pcurve",
                       fluidRow(
                           box(title = "Setup power curves",status = "warning",
                               column(12,
                                      h3("Subjects"),
                                      splitLayout(
                                        numericInput("pcurve_n2_min", "Min", 2),
                                        numericInput("pcurve_n2_max", "Max", 30),
                                        numericInput("pcurve_n2_increment", "Increment", 3)
                                        ),
                                      conditionalPanel("input.multi_des == '3lvl'",
                                       h3("Clusters")),
                                      conditionalPanel("input.multi_des == '3lvl'",
                                      textInput("pcurve_n3", "Number of clusters", "2,4,8,12")),
                                      actionButton('pcurve_go_button', 'Draw plot', class = "pull-right")
                                      )
                                )

                       ),
                       fluidRow(

                            box(plotOutput("plot_power_curve", height = 400),
                                downloadButton("download_power_plot", "Download plot"))
                       ),
                       fluidRow(
                       box(tableOutput("table_power_curve"),
                           downloadButton("download_power_csv", "Download as CSV"))


                       )
               )
          )
     )
)
