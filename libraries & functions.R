libraries = c('ftsa', 'fda.usc', 'caret', 'ROCR', 'pROC', 'fda', 'RColorBrewer')

lapply(libraries, function(x) if (!(x %in"% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quetly = TRUE, character.only = TRUE)