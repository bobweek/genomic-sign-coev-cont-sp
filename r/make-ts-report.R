Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/quarto/bin")
rmarkdown::render('r/time-series-report.Rmd',output_file='time-series/time-series-report.html')
