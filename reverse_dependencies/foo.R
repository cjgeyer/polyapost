# clean
foo <- list.files()
bar <- grep("^foo\\.R", foo, invert = TRUE, value = TRUE)
baz <- paste("rm -r", paste(bar, collapse = " "))
system(baz)
list.files()
# my package
system("R CMD build ../package/polyapost")
# doit
options(repos = "https://cloud.r-project.org/")
library(tools)
out <- check_packages_in_dir(".", reverse = "all")
warnings()
summary(out)
summarize_check_packages_in_dir_results(".")
summarize_check_packages_in_dir_timings(".")
