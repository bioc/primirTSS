library(shiny)
library(readr)

source(system.file("extdata", "primirTSSapp", "extra_fun", "shiny_merge.R", package = "primirTSS"))
options(shiny.maxRequestSize=60*1024^3)

function(input, output) {

  wm_eSet <- eventReactive(input$wm_start, {

	withProgress(
	  message = 'Calculation in progress', detail = 'This may take a while...', value = 0, {

	  incProgress(0.3)
	    if (input$file1_test == "no") {
	      histone_bed <- "NO"
	    } else {
	      histone_bed <- read_csv(input$file1$datapath, guess_max = 100000000)
	    }

	    if (input$file2_test == "no") {
	      pol2_bed <- "NO"
	    } else {
	      pol2_bed <- readr::read_csv(input$file2$datapath, guess_max = 100000000)
	    }

	    if (histone_bed == "NO" && pol2_bed == "NO") {
	      stop("you should input as least one of histone peaks file or pol2 peaks file")
	    } else if (histone_bed == "NO" && pol2_bed != "NO") {
	      bed_merged <- env_shiny$peak_merge2(pol2_bed, n = input$n2)
	    } else if (histone_bed != "NO" && pol2_bed == "NO") {
	      bed_merged <- env_shiny$peak_merge2(histone_bed, n = input$n1)
	    } else if (histone_bed != "NO" && pol2_bed != "NO") {
	      bed_merged_h3 <- env_shiny$peak_merge2(histone_bed, n = input$n1)
	      bed_merged_pol2 <- env_shiny$peak_merge2(pol2_bed, n = input$n2)
	      bed_merged <- peak_join(bed_merged_h3, bed_merged_pol2)
	    }

	    if (input$file3_test == "no") {
	      expressed_mir <- "all"
	    } else {
	      expressed_mir <- readLines(input$file3$datapath)
	    }

	    if (input$file4_test == "yes") {
	      ignore_DHS_check <- FALSE
	      DHS <- readr::read_csv(input$file4$datapath, guess_max = 100000000)
	    } else {
	      ignore_DHS_check <- TRUE
	    }

	    if (input$DHS_check_result == "no") {
	      allmirdhs_byforce <- FALSE
	    } else {
	      allmirdhs_byforce <- TRUE
	    }

	    if (input$file5_test == "no") {
	      expressed_gene <- "all"
	    } else {
	      expressed_gene <- readLines(input$file5$datapath)
	    }

	    if (input$gene_check_result == "no") {
	      allmirgene_byforce <- FALSE
	    } else {
	      allmirgene_byforce <- TRUE
	    }

	  incProgress(0.6)

	  tss_result <- find_tss(bed_merged, expressed_mir = expressed_mir,
	                         flanking_num = input$flanking_num,
	                         threshold = input$threshold,
	                         ignore_DHS_check = ignore_DHS_check,
	                         DHS = DHS,
	                         allmirdhs_byforce = allmirdhs_byforce,
	                         expressed_gene = expressed_gene,
	                         allmirgene_byforce = allmirgene_byforce)

	  incProgress(0.8)

	  if (input$tf == "yes") {
	    a <- mir_tf(tss_result$tss_df, input$tf_n, input$min.score)
	  }

	})
    if (input$tf == "yes") {
      result <- dplyr::left_join(tss_result$tss_df, a, by = "mir_name")
    } else {
      result <- tss_result$tss_df
    }

    result
  })

  output$wm_exprs_mx <- renderTable({
    head(wm_eSet())
  },
  rownames = FALSE,
  digits = 0)

  output$tss_result.csv <- downloadHandler(
    filename = 'tss_result.csv',
    content = function(file) {
      readr::write_csv(wm_eSet(), file)
    }
  )


  wm_eSet_2 <- eventReactive(input$wm_start.2, {

    withProgress(
      message = 'Calculation in progress', detail = 'This may take a while...', value = 0, {

        incProgress(0.3)
        if (input$file1.2_test == "no") {
          histone_bed <- "NO"
        } else {
          histone_bed <- read_csv(input$file1.2$datapath, guess_max = 100000000)
        }

        if (input$file2.2_test == "no") {
          pol2_bed <- "NO"
        } else {
          pol2_bed <- readr::read_csv(input$file2.2$datapath, guess_max = 100000000)
        }

        if (histone_bed == "NO" && pol2_bed == "NO") {
          stop("you should input as least one of histone peaks file or pol2 peaks file")
        } else if (histone_bed == "NO" && pol2_bed != "NO") {
          bed_merged <- env_shiny$peak_merge2(pol2_bed, n = input$n2.2)
        } else if (histone_bed != "NO" && pol2_bed == "NO") {
          bed_merged <- env_shiny$peak_merge2(histone_bed, n = input$n1.2)
        } else if (histone_bed != "NO" && pol2_bed != "NO") {
          bed_merged_h3 <- env_shiny$peak_merge2(histone_bed, n = input$n1.2)
          bed_merged_pol2 <- env_shiny$peak_merge2(pol2_bed, n = input$n2.2)
          bed_merged <- peak_join(bed_merged_h3, bed_merged_pol2)
        }


        if (input$file4.2_test == "yes") {
          ignore_DHS_check <- FALSE
          DHS <- readr::read_csv(input$file4.2$datapath, guess_max = 100000000)
        } else {
          ignore_DHS_check <- TRUE
        }

        if (input$DHS_check_result.2 == "no") {
          allmirdhs_byforce <- FALSE
        } else {
          allmirdhs_byforce <- TRUE
        }

        if (input$file5.2_test == "no") {
          expressed_gene <- "all"
        } else {
          expressed_gene <- readLines(input$file5.2$datapath)
        }

        if (input$gene_check_result.2 == "no") {
          allmirgene_byforce <- FALSE
        } else {
          allmirgene_byforce <- TRUE
        }

        incProgress(0.6)
        track <- plot_primiRNA(input$file3.2, bed_merged,
                                     flanking_num = input$flanking_num.2,
                                     threshold = input$threshold.2,
                                     ignore_DHS_check,
                                     DHS, allmirdhs_byforce,
                                     expressed_gene = expressed_gene,
                                     allmirgene_byforce)
      })

    track
  })

  output$wm_plot <- renderPlot({
    wm_eSet_2()

  })

}







