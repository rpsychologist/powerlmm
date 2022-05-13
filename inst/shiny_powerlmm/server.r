library(shiny)
library(shinydashboard)
library(tidyr)
library(viridis)
library(ggplot2)
library(ggsci)



get_sample_size <- function(paras) {

    n2 <- powerlmm:::get_n2.plcp(paras)$treatment
    n3 <- powerlmm:::get_n3.plcp(paras)
    if(length(n2) == 1) {
        n2 <- rep(n2, n3$treatment)
    }
    n2_tx <- n2
    if(paras$partially_nested) {
        n2_cc <- c(sum(n2), rep(NA, length(n2) - 1))
    } else n2_cc <- n2

    data.frame(cluster = c(1:n3$treatment, "total"),
               control = c(n2_cc, sum(n2_cc, na.rm=TRUE)),
               treatment = c(n2_tx, sum(n2_tx)))
}

shinyServer(function(input, output, session) {

    ## UI
    output$ui_slope_lvl2 <- renderUI({
        if(input$use_standardized) {
            NULL
        } else {
            numericInput("sigma.person.slope", "Slope subjects", 0)
        }
    })
    output$ui_slope_lvl3 <- renderUI({
        if(input$use_standardized) {
            numericInput("icc_slope", "% slope variance at level 3", 0.05)
        } else {
            numericInput("sigma.cluster.slope", "Cluster subjects", 0.14)
        }
    })
    output$ui_var_ratio <- renderUI({
        if(input$use_standardized) {
            numericInput("var_ratio", "Variance ratio", 0.1)
        } else {
            NULL
        }
    })
    observe({
        if(input$multi_des == "2lvl") {
            updateNumericInput(session, "sigma.cluster.intercept", value = 0)
            updateNumericInput(session, "sigma.cluster.slope", value = 0)
            updateNumericInput(session, "icc_slope", value = 0)
            updateTextInput(session, "n2", label = "n2 (subjects per treatment)", value = 25)
            updateNumericInput(session, "n3", value = 1)
        } else {
            updateTextInput(session, "n2", label = "n2 (subjects per cluster)", value = 5)
            updateNumericInput(session, "n3", value = 4)
        }

        #if(length(as.numeric(unlist(strsplit(input$n2,",")))) > 1) {
        #    updateNumericInput(session, "n3", label = "n3 (clusters, unbalanced)", value = length(unlist(paras$n2)))
        #}
        p$run <- FALSE
        })

    observe({
        n2 <- as.numeric(unlist(strsplit(input$n2,",")))
        if(length(n2) > 1) {
            updateNumericInput(session, "n3", label = "n3 (clusters, unbalanced)", value = length(n2))
        }
    })


    p <- reactiveValues(paras = NULL, run = FALSE)
  observeEvent(c(input$button1,
                 input$button2,
                 input$button3,
                 input$button4), {
      n2 <- as.numeric(unlist(strsplit(input$n2,",")))
      if(length(n2) > 1) {
          n2 <- unequal_clusters(n2)
      }
      p_tx <- dropout_weibull(input$retention_tx, input$gamma_tx)
      p_cc <- dropout_weibull(input$retention_cc, input$gamma_cc)
      if(input$retention_tx == 0 & input$retention_cc == 0) {
          retention <- 0
      } else {
          retention <- per_treatment(p_cc, p_tx)
      }

      if(input$use_standardized == "yes") {
          validate(
              need((input$icc_pre_persons + input$icc_pre_clusters) < 100,
                   "Input '% baseline variance at level 2' and '% baseline variance at level 3' can't sum to over 100 %",
                   "")
          )
          paras <- study_parameters(n1 = input$n1,
                                    n2 = n2,
                                    n3 = input$n3,
                                    T_end = input$T_end,
                                    fixed_intercept = input$fixed.intercept,
                                    fixed_slope = -input$fixed.slope/input$T_end,
                                    icc_pre_subject = input$icc_pre_persons/100,
                                    icc_pre_cluster = input$icc_pre_clusters/100,
                                    icc_slope = input$icc_slope/100,
                                    var_ratio = input$var_ratio,
                                    cor_subject = input$cor.person,
                                    cor_cluster = input$cor.cluster,
                                    cor_within = 0,
                                    cohend = -input$cohend,
                                    partially_nested = as.logical(input$partially_nested),
                                    dropout = retention)
          p$paras <- paras
          updateNumericInput(session, "sigma.person.intercept", value = p$paras$sigma_subject_intercept)
          updateNumericInput(session, "sigma.person.slope", value = p$paras$sigma_subject_slope)
          updateNumericInput(session, "sigma.cluster.intercept", value = p$paras$sigma_cluster_intercept)
          updateNumericInput(session, "sigma.cluster.slope", value = p$paras$sigma_cluster_slope)
          updateNumericInput(session, "sigma.error", value = p$paras$sigma_error)


      } else {

          paras <- study_parameters(n1 = input$n1,
                                    n2 = n2,
                                    n3 = input$n3,
                                    T_end = input$T_end,
                                    fixed_intercept = input$fixed.intercept,
                                    fixed_slope =  -input$fixed.slope/input$T_end,
                                    sigma_subject_intercept = input$sigma.person.intercept,
                                    sigma_cluster_intercept = input$sigma.cluster.intercept,
                                    sigma_cluster_slope = input$sigma.cluster.slope,
                                    sigma_subject_slope = input$sigma.person.slope,
                                    sigma_error = input$sigma.error,
                                    cor_subject = input$cor.person,
                                    cor_cluster = input$cor.cluster,
                                    cor_within = 0,
                                    effect_size = -input$cohend,
                                    partially_nested = as.logical(input$partially_nested),
                                    dropout = retention)
          p$paras <- paras
          updateNumericInput(session, "icc_pre_persons", value = get_ICC_pre_subjects(p$paras)*100)
          updateNumericInput(session, "icc_pre_clusters", value = get_ICC_pre_clusters(p$paras)*100)
          updateNumericInput(session, "icc_slope", value = get_ICC_slope(v1 = p$paras$sigma_cluster_slope,
                                                                         u1 = p$paras$sigma_person_slope)*100)
          updateNumericInput(session, "var_ratio", value = get_var_ratio(v1 = p$paras$sigma_cluster_slope,
                                                                         u1 = p$paras$sigma_subject_slope,
                                                                         error = p$paras$sigma_error))
      }
      p$run <- TRUE
  })



    output$plot_ES <- renderPlot({
        if(p$run == FALSE) return()
        plot(p$paras, type = "effect") + theme_minimal()

    })
    output$plot_retention <- renderPlot({
        if(p$run == FALSE) return()
        plot(p$paras, type = "dropout") + theme_minimal()

    })
    output$table_retention <- renderTable({
        if(p$run == FALSE) return()
        get_dropout(p$paras)
    })
    output$plot_sds <- renderPlot({
        if(p$run == FALSE) return()
        p <- get_sds(p$paras)
        p <- plot(p)

        p + theme_minimal()

    })
    output$table_sds <- renderTable({
        if(p$run == FALSE) return()
          x <-  get_sds(p$paras)

          x[,c("time","SD_with_random_slopes")]

    })

    output$plot_VPC <- renderPlot({
        if(p$run == FALSE) return()
        p <- get_VPC(p$paras)
        p <- plot(p)

        p + theme_minimal()

    })
    output$table_VPC <- renderTable({
        if(p$run == FALSE) return()
        get_VPC(p$paras)



    })
    power_table <- eventReactive(input$pcurve_go_button, {
        n2_min <- input$pcurve_n2_min
        n2_max <- input$pcurve_n2_max
        n2_increment <- input$pcurve_n2_increment

        if(input$multi_des == "2lvl") {
            n3_list <- 1
        } else {
            n3_list <- as.numeric(unlist(strsplit(input$pcurve_n3,",")))
        }

        n2 <- seq(n2_min, n2_max, by = n2_increment)
        n3 <- n3_list

        progress <- shiny::Progress$new(max = length(n2) * length(n3))
        progress$set(message = "Computing power", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        updateProgress <- function(value = NULL, detail = NULL) {
            if (is.null(value)) {
                value <- progress$getValue()
                value <- value + (progress$getMax() - value) / 5
            }
            progress$set(value = value)
        }



        p <- p$paras
        ptm <- proc.time()

        res <- get_power_table(p,
                               n2 = n2,
                               n3 = n3, updateProgress = updateProgress)
        print(proc.time() - ptm)
        res
    })
    power_plot <- eventReactive(input$pcurve_go_button, {
        plot(power_table())
    })

    output$plot_power_curve <- renderPlot({
        input$pcurve_go_button

        power_plot()

    })
    output$table_power_curve <- renderTable({
        res <- power_table()
        res <- tidyr::spread(res, dropout, power)
        if(input$multi_des == "2lvl") {
            res <- dplyr::select(res, -n3)
        }
        res
        }
    )

    output$download_power_csv <- downloadHandler(
        filename = "power_curve.csv",
        content = function(file) {
            res <- power_table()
            res <- tidyr::spread(res, dropout, power)
            write.csv(res, file)
        }
    )
    output$download_power_plot <- downloadHandler(
        filename = "power_curve.png",
        content = function(file) {
            ggsave(file, plot = power_plot(), device = "png")
        }
    )

    output$table_VPC <- renderTable({
        if(p$run == FALSE) return()
        get_VPC(p$paras)
        })

    output$table_correlation_matrix <- renderTable({
        if(p$run == FALSE) return()
        res <- get_correlation_matrix(p$paras)
        rnames <-  rownames(res)
        res <- as.data.frame(res)
        #res <- cbind(' ' = rnames, res)
        res
    }, rownames = TRUE
    )

    output$plot_correlation_matrix <- renderPlot( {
        if(p$run == FALSE) return()
        res <- get_correlation_matrix(p$paras)
        res <- reshape2::melt(res)
        ggplot(res, aes(Var1, Var2, color = value, fill = value)) + geom_tile() +
            geom_text(aes(label = round(value,2)), hjust = "center", color = "black") +
            scale_fill_viridis() + scale_color_viridis() +
            labs(color = "correlation", fill = "correlation",
                 x = "Time", y = "Time",
                 title = "Subject-level correlation matrix") +
            theme_minimal()
    })


    output$table_n2 <- renderTable({
        if(p$run == FALSE) return()
        get_sample_size(p$paras)
    }, digits = 1)

    output$total_n <- renderValueBox({
        if(p$run == FALSE) {
            tmp <- "-"
        } else {
            tmp <- powerlmm:::get_tot_n.plcp(p$paras)$total
        }

        valueBox(
            tmp, "total sample size", icon = icon("list"),
            color = "purple"
        )

    })


    output$power_box <- renderValueBox({
        if(p$run == FALSE) {
            tmp <- "-"
        } else {
         tmp <- paste0(round(get_power(p$paras)$power,2)*100, "%")
        }
        valueBox(
            tmp, "Power", icon = icon("flash"),
            color = "green"
        )
    })


    output$img_study_design <- renderImage({
        if(input$multi_des == "3lvl") {
            if(!as.logical(input$partially_nested)) {
                filename <- normalizePath(file.path('www/img/three-level.png'))
                return(list(
                    src = filename,
                    contentType = "image/png",
                    alt = "Power calculations for three-level linear mixed models with missing data"
                ))
            } else {
                filename <- normalizePath(file.path('www/img/partially-nested.png'))
                return(list(
                    src = filename,
                    contentType = "image/png",
                    alt = "Power calculations for partially-nested three-level linear mixed models with missing data"
                ))
            }

        } else {
            filename <- normalizePath(file.path('www/img/two-level.png'))
            return(list(
                src = filename,
                contentType = "image/png",
                alt = "Power calculations for two-level linear mixed models with missing data"
            ))
        }

    }, deleteFile = FALSE)


    output$debug <- renderPrint({
       print(p$paras)
    }
    )



})
