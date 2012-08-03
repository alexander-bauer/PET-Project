require(knitr)
knit("FinalReport.Rmd")
system("pandoc -s -S -i -t slidy --mathjax -o index.html FinalReport.md")
